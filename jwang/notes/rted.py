import gurobipy as gb
import pandas as pd
import numpy as np
from andes.interop.pandapower import make_GSF
import logging
logger = logging.getLogger(__name__)


class rted:
    """
    Base class of RTED model, where DCOPF with SFR is modeled by Gurobipy.

    Parameters
    ----------
    ssp : pandapower instance

    Attributes
    ----------
    - `rted`: the DCOPF with SFR model
    - `dcopf`: the standard DCOPF
    - `{device}data`: DataFrame of {device}, {device} includes gen, line, load
    - `{device}dict`: Dict of {device}, {device} includes gen, line, load

    Notes
    -----
      - In ``rted``, there are two types of generators, where the type determines the generator limtis equation
      - All units are in p.u., where base_mva is the same of the input ``ssp.sn_mva``
      - Offline generators will be limited to zero, uncontrollable generators will be limited to ``p0``
      - All geenrators are set as Type I by default, which can be altered by ``def_typeII()`` after instantiated
      - Type I generator capacity limit: p_sch - prd >= pmin, p_sch + pru <= pmax
      - Type II generator capacity limit: prd <= prd_max, pru <= pru_max
    """

    def __init__(self, ssp):
        """
        """
        # --- store original info of ssp ---
        # --- def. data ---
        # --- 0. gsf ---
        busdata = ssp.bus[['name']]
        busdata = self.set_idx(busdata, 'bus')
        gsf = make_GSF(ssp)
        self.gsfdata = pd.DataFrame(gsf)
        self.gsfdata = self.set_idx(self.gsfdata, 'line')
        # self.gsfdata.drop(['name'], axis=1, inplace=True)

        # --- 1. gen ---
        ssp_gen = \
            pd.merge(left=ssp.poly_cost.rename(columns={'element': 'index'}),
                     right=ssp.gen.reset_index(),
                     on='index', how='right')
        ssp_gen.fillna(0, inplace=True)
        self.gendata = pd.DataFrame()
        # --- a. controllable gen cost, quadratic ---
        self.gendata['gen_pp'] = ssp_gen.name
        self.gendata['cp1'] = ssp_gen.cp1_eur_per_mw.tolist()
        self.gendata['cp2'] = ssp_gen.cp2_eur_per_mw2.tolist()
        self.gendata['cp0'] = ssp_gen.cp0_eur.tolist()
        # --- b. controllable gen bus ---
        self.gendata['bus'] = ssp_gen.bus.tolist()
        # --- c. alter gen limit, 0 for offline gen, constant for uncontrolled gen ---
        self.gendata['ctrl'] = ssp_gen.controllable.astype(int)
        ssp_gen['unctrl'] = 1 - ssp_gen.controllable
        self.gendata['p0'] = ssp.gen.p_mw / ssp.sn_mva
        pmax = ssp_gen.in_service * (ssp_gen.controllable * ssp_gen.max_p_mw + ssp_gen.unctrl * ssp_gen.p_mw)
        pmin = ssp_gen.in_service * (ssp_gen.controllable * ssp_gen.min_p_mw + ssp_gen.unctrl * ssp_gen.p_mw)
        self.gendata['pmax'] = list(pmax / ssp.sn_mva)
        self.gendata['pmin'] = list(pmin / ssp.sn_mva)
        # --- d. gen type ---
        # set all gen as type I by default
        # Type I ramp: pmin =< p_sch +- pr <= pmax
        # Type II ramp: prdmax =< pr <= prumax
        self.gendata['type'] = 'I'

        self.ps_defined = False
        # SFR cost
        self.gendata['c_ru'] = 0
        self.gendata['c_rd'] = 0
        # ramp limits
        self.gendata['rampu'] = 10
        self.gendata['rampd'] = 10

        self.var_defined = False
        # previous setpoints
        self.gendata['p_pre'] = 0
        # SFR limits
        self.gendata['pru_max'] = 0
        self.gendata['prd_max'] = 0

        # --- set index ---
        self.gendata = self.set_idx(self.gendata, 'gen')
        # convert to a dict
        self.gendict = self.gendata.T.to_dict()

        # --- 1.b Generation Shift Factor ---
        gsfT = self.gsfdata.T
        gsfT['bus'] = gsfT.index
        self.gen2data = self.gendata.merge(gsfT, on='bus', how='left')
        self.gen2data.index = self.gendata.index
        self.gen2dict = self.gen2data.T.to_dict()

        # --- 0. line ---
        ssp_line = pd.merge(ssp.bus[['vn_kv']].reset_index(),
                            ssp.line[['from_bus', 'max_i_ka']].rename(columns={'from_bus': 'index'}),
                            on='index', how='right')
        ssp_line['lim'] = ssp_line.max_i_ka * ssp_line.vn_kv / ssp.sn_mva
        ssp_line.sort_index(inplace=True)
        self.linedata = ssp_line[['lim']]
        # add power surplus, where the controlled gen is removed
        sup = pd.DataFrame()
        sup['bus'] = ssp.bus.index
        sup = sup.merge(ssp.gen[['bus', 'p_mw', 'controllable']],
                        'left', 'bus').fillna(0).rename(columns={'p_mw': 'gen'})
        sup = sup.merge(ssp.load[['bus', 'p_mw']], 'left', 'bus').fillna(0).rename(columns={'p_mw': 'load'})
        sup['uc'] = 1 - sup.controllable
        sup['net'] = sup.uc * (sup.gen-sup.load) / ssp.sn_mva
        sup2 = sup.groupby('bus').sum().reset_index(drop=True)
        sup2['bus'] = sup2.index
        sup2 = sup[['bus', 'net']].groupby('bus').sum().reset_index(drop=True)
        sup2['bus'] = sup2.index
        self.linedata['sup'] = np.matmul(gsf, sup2.net.values)
        # --- set index ---
        self.linedata = self.set_idx(self.linedata, 'line')
        self.linedict = self.linedata.T.to_dict()

        # --- 3. load: forcated load, fluc. up and down ---
        self.ptotal = ssp.load.p_mw.sum() / ssp.sn_mva
        self.dpd_defined = False
        self.dpd_u = 0.2
        self.dpd_d = 0.2

    def set_idx(self, df, name='device'):
        """
        Set the index of input DataFrame as string with format
        f'{name}{index}'.

        Parameters
        ----------
        df: DataFrame
            input DataFrame
        """
        dfout = df.copy()
        dfout['idx'] = list(range(1, dfout.shape[0]+1))
        dfout['name_c'] = [name] * dfout.shape[0]
        dfout['name'] = dfout['name_c'] + dfout['idx'].astype(str)
        dfout.set_index(dfout['name'], inplace=True)
        dfout.drop(['idx', 'name_c', 'name'], axis=1, inplace=True)
        return dfout

    def set_list(self, inlist, dfname, attr):
        """Input the list into a datafram column, with len check"""
        df = getattr(self, dfname)
        if len(inlist) > 0:
            if len(inlist) != len(df):
                logger.error(f'Input list does not match {dfname} length!')
            df[f'{attr}'] = inlist
        return True

    def def_type(self, gen=[]):
        """
        Define which generators are type II.

        Parameters
        ----------
        gen : list
            index of generators (in attribute gendata)
        """
        self.gendata.type.loc[gen] = 'II'
        self.gendict = self.gendata.T.to_dict()
        return True

    def def_var(self, ptotal=10, dpd_u=0.2, dpd_d=0.2, p_pre=[], gen=[], pru_max=[], prd_max=[]):
        """
        Define model parameters.

        Parameters
        ----------
        ptotal : number
            Total load (p.u.)
        dpd_u : number
            SFR Up requirement (p.u.)
        dpd_d : number
            SFR Down requirement (p.u.)
        p_pre : list
            Previous generator setpoints (p.u.)
        gen : list
            index of generators, only used for Type II generators
        pru_max : list
            Max SFR Up capacity for Type II gen (p.u.)
        prd_max : list
            Max SFR Down capacity for Type II gen (p.u.)
        """
        # load
        self.ptotal = ptotal
        # --- SFR requirements ---
        self.dpd_u = dpd_u
        self.dpd_d = dpd_d

        # --- gen data ---
        self.set_list(p_pre, 'gendata', 'p_pre')
        # --- Type II generators---
        self.gendata['pru_max'] = 0
        self.gendata['prd_max'] = 0
        if gen:
            pass
        else:
            gendict_II = dict()
            for (new_key, new_value) in self.gendict.items():
                if new_value['type'] == 'II':
                    gendict_II[new_key] = new_value
            gen = gendict_II.keys()
        # --- length check ---
        if len(gen) != len(pru_max):
            logger.warning('pru_max length does not match Type II generator numbers!')
        if len(gen) != len(prd_max):
            logger.warning('prd_max length does not match Type II generator numbers!')

        if len(pru_max) > 0:
            pass
        else:
            pru_max = [0] * len(gen)

        if len(prd_max) > 0:
            pass
        else:
            prd_max = [0] * len(gen)

        self.gendata['pru_max'].loc[gen] = pru_max
        self.gendata['pru_max'].loc[gen] = prd_max
        self.gendict = self.gendata.T.to_dict()

        self.var_defined = True
        return self.var_defined

    def def_ps(self, c_ru=[], c_rd=[], rampu=[], rampd=[]):
        """
        Define power system parameters.

        Parameters
        ----------
        c_ru : list
            SFR Up cost
        c_rd : list
            SFR Down cost
        rampu : list
            Ramp Up limit
        rampd : list
            Ramp Down limit
        """
        # --- gen SFR cost, ramp---
        self.set_list(c_ru, 'gendata', 'c_ru')
        self.set_list(c_rd, 'gendata', 'c_rd')
        self.set_list(rampu, 'gendata', 'rampu')
        self.set_list(rampd, 'gendata', 'rampd')
        self.gendict = self.gendata.T.to_dict()

        self.ps_defined = True
        return self.ps_defined

    def _build_vars(self, mdl, gendict):
        GEN = gendict.keys()
        # --- 1. gen: p_sch, pr_gu, pr_gd ---
        self.p_sch = mdl.addVars(GEN, name='p_sch', vtype=gb.GRB.CONTINUOUS, obj=0,
                                 ub=self.gendata.pmax.tolist(), lb=self.gendata.pmin.tolist())
        self.p_ru = mdl.addVars(GEN, name='p_ru', vtype=gb.GRB.CONTINUOUS, obj=0)
        self.p_rd = mdl.addVars(GEN, name='p_rd', vtype=gb.GRB.CONTINUOUS, obj=0)
        return mdl

    def _build_obj(self, mdl, gendict):
        GEN = gendict.keys()

        cost_pg = sum(self.p_sch[gen] * gendict[gen]['cp1'] * gendict[gen]['ctrl']
                      + self.p_sch[gen] * self.p_sch[gen] * gendict[gen]['cp2'] * gendict[gen]['ctrl']
                      + gendict[gen]['cp0'] * gendict[gen]['ctrl']
                      for gen in GEN)
        cost_pru = sum(self.p_ru[gen] * gendict[gen]['c_ru'] for gen in GEN)
        cost_prd = sum(self.p_rd[gen] * gendict[gen]['c_rd'] for gen in GEN)

        self.obj = mdl.setObjective(expr=cost_pg + cost_pru + cost_prd, sense=gb.GRB.MINIMIZE)
        return mdl

    def _build_cons(self, mdl, gendict, gen2dict, linedict, dpd_u, dpd_d, ptotal):

        # --- filter Type II gen ---
        gendict_I = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 'I':
                gendict_I[new_key] = new_value
        gendict_II = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 'II':
                gendict_II[new_key] = new_value

        GEN = gendict.keys()
        GENI = gendict_I.keys()
        GENII = gendict_II.keys()
        LINE = linedict.keys()

        # --- 4. power balance ---
        p_sum = sum(self.p_sch[gen] for gen in GEN)
        mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- 1.a Type I GEN capacity limits ---
        mdl.addConstrs((self.p_sch[gen] + self.p_ru[gen] <= gendict[gen]['pmax'] for gen in GENI),
                       name='Pg_max')
        mdl.addConstrs((self.p_sch[gen] - self.p_rd[gen] >= gendict[gen]['pmin'] for gen in GENI),
                       name='Pg_min')
        # --- 1.b Type II Gen capacity and SFR limits---
        mdl.addConstrs((self.p_sch[gen] <= gendict[gen]['pmax'] for gen in GENII),
                       name='Pg_max')
        mdl.addConstrs((self.p_sch[gen] >= gendict[gen]['pmin'] for gen in GENII),
                       name='Pg_min')
        mdl.addConstrs((self.p_ru[gen] <= gendict[gen]['pru_max'] for gen in GENII),
                       name='PRU_max')
        mdl.addConstrs((self.p_rd[gen] <= gendict[gen]['prd_max'] for gen in GENII),
                       name='PRD_max')

        # --- 2. SFR requirements ---
        # --- a) RegUp --
        mdl.addConstr(sum(self.p_ru[gen] for gen in GEN) == dpd_u, name='RegUp')
        # --- b) RegDn --
        mdl.addConstr(sum(self.p_rd[gen] for gen in GEN) == dpd_d, name='RegDn')

        # --- 5. line limits ---
        for line in LINE:
            lhs1 = sum(self.p_sch[gen] * gen2dict[gen][line] for gen in GEN)
            mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['lim'], name=f'{line}_U')
            mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['lim'], name=f'{line}_D')

        # --- 6. ramp limits ---
        mdl.addConstrs((self.p_sch[gen] - gendict[gen]['p_pre'] <= gendict[gen]['rampu']
                       for gen in GEN), name='RampU')
        mdl.addConstrs((gendict[gen]['p_pre'] - self.p_sch[gen] <= gendict[gen]['rampd']
                       for gen in GEN), name='RampD')

        return mdl

    def build_rted(self):
        """Build RTED model (DCOPF with SFR)"""

        if not self.var_defined:
            logger.warning('PS vars ``dpd_u``, ``dpd_d``, ``p_pre``, ``pru_max``, ``prd_max`` are not given.'
                           'Default values will be used. They can be defined by ``def_var``:')
        if not self.ps_defined:
            logger.warning('PS params ``c_ru``, ``c_rd``, ``rampu``, ``rampd`` are not given.'
                           'Default values will be used. They can be defined by ``def_ps``')

        self.rted = gb.Model('DCOPF+SFR')

        # --- A. decision vars ---
        self.rted = self._build_vars(self.rted, self.gendict)

        # --- B. obj. ---
        self.rted = self._build_obj(self.rted, self.gendict)

        # --- C. constraints ---
        self.rted = self._build_cons(self.rted, self.gendict, self.gen2dict, self.linedict,
                                     self.dpd_u, self.dpd_d, self.ptotal)
        return True

    def build_sdc(self):
        """Build standard DCOPF"""
        self.dcopf = gb.Model('StdDCOPF')

        GEN = self.gendict.keys()

        # --- A. decision vars ---
        self.dcopf = self._build_vars(self.dcopf, self.gendict)

        # --- B. obj. ---
        self.dcopf = self._build_obj(self.dcopf, self.gendict)

        # --- C. constraints ---
        # set the SFR requirements as 0
        self.dcopf = self._build_cons(self.dcopf, self.gendict, self.gen2dict, self.linedict, 0, 0, self.ptotal)

        # --- limit SFR vars to 0---
        self.dcopf.addConstrs((self.p_ru[gen] == 0 for gen in GEN), name='PRU_0')
        self.dcopf.addConstrs((self.p_rd[gen] == 0 for gen in GEN), name='PRD_0')
        return True

    def get_res(self, model='rted'):
        """
        Get optimization resutlts, can be used after effectively solving the optimization problem.

        Parameters
        ----------
        model : str

            'rted' or 'dcopf'

        Returns
        ----------
        DataFrame

            The output DataFrame contains setpoints ``p_sch``,
            SFR Up/Dn capacity ``pru`` and ``prd``,
            and AGC Up/Dn participation factor: ``bu`` and ``bd``
        """
        # --- gather data --
        pru = []
        prd = []
        p_sch = []
        for gen in self.gendict.keys():
            p_sch.append(self.p_sch[gen].X)
            pru.append(self.p_ru[gen].X)
            prd.append(self.p_rd[gen].X)
        # --- build output table ---
        dcres = pd.DataFrame()
        dcres['gen_pp'] = self.gendata.gen_pp
        dcres['p_sch'] = p_sch
        dcres['pru'] = pru
        dcres['prd'] = prd
        dcres['bu'] = dcres['pru'] / dcres['pru'].sum()
        dcres['bd'] = dcres['prd'] / dcres['prd'].sum()
        dcres.fillna(0, inplace=True)

        # --- cost ---
        mdl = getattr(self, model)
        total_cost= mdl.getObjective().getValue()

        sfr_cost = np.sum(dcres['pru'] * self.gendata.c_ru + dcres['prd'] * self.gendata.c_rd)
        logger.warn(
            f'{model} cost (p.u.): Total={np.round(total_cost, 3)}, '
            + f'GEN={np.round(total_cost-sfr_cost, 3)}, SFR={np.round(sfr_cost, 3)}')
        return dcres
