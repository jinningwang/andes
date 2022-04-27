from andes.interop.pandapower import to_pandapower
from andes.interop.pandapower import make_GSF, build_group_table
import gurobipy as gb
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)


class system:
    """
    Base class of jams system.

    Parameters
    ----------
    name : str
        Name of the system.

    Attributes
    ----------
    bus: pandas.DataFrame
        Bus data.
    gen: pandas.DataFrame
        Generator data.
    line: pandas.DataFrame
        Line data.
    load: pandas.DataFrame
        Load data.
    gen_gsf: pandas.DataFrame
        Generator shift factor of gen bus data.
    cost: pandas.DataFrame
        Cost data.
    """

    def __init__(self, name='system'):
        self.name = name

    def update_dict(self, model=None):
        """
        Update model DataFrame into model dict.

        Parameters
        ----------
        model : list
            list of models that need to be updated.
            If None is given, update all models.
        """
        # --- validity check ---
        if not hasattr(self, 'cost'):
            self._default_cost()

        # --- build dict ---
        if not model:
            mdl_list = ['bus', 'gen', 'line', 'gen_gsf', 'cost']
        else:
            mdl_list = model
        for mdl in mdl_list:
            mdl_df = getattr(self, mdl).copy()
            mdl_df.set_index(['idx'], inplace=True)
            setattr(self, mdl+'dict', mdl_df.T.to_dict())

    def from_andes(self, ssa):
        """
        Create jams system from ANDES system.

        Parameters
        ----------
        ssa : andes.system.system
            ANDES system.

        Notes
        -----
        All generators are set as controllable.
        """
        # --- base mva ---
        self.mva = ssa.config.mva

        # --- bus ---
        bus_cols = ['idx', 'u', 'name', 'Vn', 'vmax', 'vmin', 'v0', 'a0', 'area', 'zone', 'owner']
        self.bus = ssa.Bus.as_df()[bus_cols]
        self.bus.sort_values('idx', inplace=True)

        # --- generator ---
        stg_cols = ['idx', 'u', 'name', 'Sn', 'Vn', 'bus', 'p0',
                    'pmax', 'pmin', 'v0']
        self.gen = build_group_table(ssa, 'StaticGen', stg_cols).reset_index(drop=True)
        self.gen['ctrl'] = 1
        # TODO: later on, merge 'ramp5', 'ramp10', 'ramp30'
        self.gen['ramp5'] = 100
        self.gen['ramp10'] = 200
        self.gen['ramp30'] = 600
        # --- later on ---
        # self.gen['ramp5'] = self.gen['ramp5'] / self.mva
        # self.gen['ramp10'] = self.gen['ramp10'] / self.mva
        # self.gen['ramp30'] = self.gen['ramp30'] / self.mva
        # if self.gen['ramp5'].max() == 0:
        #     self.gen['ramp5'] = 100
        #     self.gen['ramp10'] = 100
        #     self.gen['ramp30'] = 100

        # --- load ---
        pq_cols = ['idx', 'u', 'name', 'bus', 'Vn', 'p0', 'q0',
                   'vmax', 'vmin', 'owner']
        self.load = ssa.PQ.as_df()[pq_cols]
        self.load.sort_values(by='idx', inplace=True)

        # --- line ---
        line_cols = ['idx', 'u', 'name', 'bus1', 'bus2', 'Sn', 'fn', 'Vn1', 'Vn2',
                     'trans', 'tap', 'phi', 'rate_a', 'rate_b', 'rate_c']
        ssa_line = ssa.Line.as_df()
        self.line = ssa_line[line_cols][ssa_line['trans'] == 0].reset_index(drop=True)
        self.load.sort_values(by='idx', inplace=True)
        if self.line['rate_a'].max() == 0:
            self.line['rate_a'] = 2000
            self.line['rate_b'] = 2000
            self.line['rate_c'] = 2000
        self.line['rate_a'] = self.line['rate_a'] / self.mva
        self.line['rate_b'] = self.line['rate_b'] / self.mva
        self.line['rate_c'] = self.line['rate_c'] / self.mva

        # --- GSF ---
        ssp = to_pandapower(ssa)
        gsf_matrix = make_GSF(ssp)
        self.gsf_matrix = gsf_matrix
        gsfdata = pd.DataFrame(gsf_matrix)
        gsfdata.columns = self.bus['idx']
        gsfdata['line'] = self.line.idx
        gsfdata.set_index('line', inplace=True)
        gsfT = gsfdata.T
        gsfT['bus'] = self.bus['idx']
        self.gen_gsf = self.gen[['idx', 'name', 'bus']].merge(gsfT, on='bus', how='left')
        self.gen_gsf.sort_values(by='idx', inplace=True)

        # add power surplus, where the controlled gen is removed
        sup = pd.DataFrame()
        sup['bus'] = self.bus['idx']
        sup = sup.merge(self.load[['bus', 'p0']],
                        on='bus', how='left').fillna(0).rename(columns={'p0': 'load'})
        sup['net'] = (-1 * sup.load)
        sup2 = sup[['bus', 'net']].groupby('bus').sum()
        self.line['sup'] = np.matmul(gsf_matrix, sup2.net.values)

        # --- update dict ---
        self.update_dict(model=None)

    def _default_cost(self):
        """
        Default cost data: c1=1, all other are 0.
        """
        self.cost = pd.DataFrame()
        self.cost['idx'] = self.gen['idx']
        self.cost['c2'] = 0
        self.cost['c1'] = 1
        self.cost['c0'] = 0
        self.cost['cr'] = 0


class dcopf(system):
    """
    DCOPF class.

    Parameters
    ----------
    name : str
        Name of the system.

    Attributes
    ----------
    bus: pandas.DataFrame
        Bus data.
    gen: pandas.DataFrame
        Generator data.
    line: pandas.DataFrame
        Line data.
    load: pandas.DataFrame
        Load data.
    gen_gsf: pandas.DataFrame
        Generator shift factor of gen bus data.
    cost: pandas.DataFrame
        Cost data.
    """

    def __init__(self, name='dcopf'):
        super().__init__(name)
        self.mdl = gb.Model(name)
        self.res_cost = 0

    def build(self):
        self.update_dict()
        # --- build DCOPF model ---
        self.mdl = self._build_vars(self.mdl)
        self.mdl = self._build_obj(self.mdl)
        self.mdl = self._build_cons(self.mdl)
        logger.warning('Successfully build DCOPF model.')

    def _build_vars(self, mdl):
        GEN = self.gendict.keys()
        # --- uncontrollable generators limit to p0 ---
        gencp = self.gen.copy()
        gencp['pmax'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        gencp['pmin'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        # --- offline geenrators limit to 0 ---
        gencp['pmax'][gencp.u == 0] = 0
        gencp['pmin'][gencp.u == 0] = 0
        # --- gen: pg ---
        self.pg = mdl.addVars(GEN, name='pg', vtype=gb.GRB.CONTINUOUS, obj=0,
                              ub=gencp.pmax.tolist(), lb=gencp.pmin.tolist())
        return mdl

    def _build_obj(self, mdl):
        GEN = self.gendict.keys()
        gendict = self.gendict
        costdict = self.costdict
        # --- minimize generation cost ---
        cost_pg = sum(self.pg[gen] * costdict[gen]['c1']
                      + self.pg[gen] * self.pg[gen] * costdict[gen]['c2']
                      + costdict[gen]['c0'] * gendict[gen]['u']  # online status
                      for gen in GEN)
        self.obj = mdl.setObjective(expr=cost_pg, sense=gb.GRB.MINIMIZE)
        return mdl

    def _build_cons(self, mdl):
        ptotal = self.load.p0.sum()

        gendict = self.gendict
        linedict = self.linedict
        gen_gsfdict = self.gen_gsfdict

        GEN = gendict.keys()
        LINE = linedict.keys()

        # --- power balance ---
        p_sum = sum(self.pg[gen] for gen in GEN)
        mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- line limits ---
        for line in LINE:
            lhs1 = sum(self.pg[gen] * gen_gsfdict[gen][line] for gen in GEN)
            mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')
        return mdl

    def get_res(self):
        """
        Get resutlts, can be used after mdl.optimize().

        Returns
        -------
        DataFrame
            The output DataFrame contains setpoints ``pg``
        """
        self.build()
        self.mdl.optimize()
        # --- check if mdl is sovled ---
        if not hasattr(self.pg[self.gen.idx[0]], 'X'):
            logger.warning('DCOPF has no valid resutls!')
            pg = [0] * self.gen.shape[0]
        else:
            logger.warning('Successfully solve DCOPF.')
            # --- gather data --
            pg = []
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
            # --- cost ---
            self.res_cost = self.mdl.getObjective().getValue()
            logger.info(f'Total cost={np.round(self.res_cost, 3)}')
        # --- build output table ---
        dcres = pd.DataFrame()
        dcres['gen'] = self.gen['idx']
        dcres['pg'] = pg
        dcres.fillna(0, inplace=True)
        return dcres


class rted(dcopf):
    """
    RTED class.

    Parameters
    ----------
    name : str
        Name of the system.

    Attributes
    ----------
    bus: pandas.DataFrame
        Bus data.
    gen: pandas.DataFrame
        Generator data.
    line: pandas.DataFrame
        Line data.
    load: pandas.DataFrame
        Load data.
    gen_gsf: pandas.DataFrame
        Generator shift factor of gen bus data.
    cost: pandas.DataFrame
        Cost data.
    """

    def __init__(self, name='rted'):
        super().__init__(name)
        # self.build()

    def from_andes(self, ssa):
        super().from_andes(ssa)
        self.gen['p_pre'] = 0
        self.gen['band'] = self.gen['pmax'] - self.gen['pmin']

    def to_dcopf(self):
        """
        Convert to DCOPF model.

        Returns
        -------
        dcopf: dcopf
            The output dcopf class.
        """
        ssd = dcopf()
        ssd.bus = self.bus
        ssd.gen = self.gen
        ssd.load = self.load
        ssd.line = self.line
        ssd.gen_gsf = self.gen_gsf
        ssd.cost = self.cost
        return ssd

    def set_p_pre(self):
        """
        Get ``p_pre`` from DCOPF.
        """
        ssd = self.to_dcopf()
        self.gen['p_pre'] = ssd.get_res()['pg']
        logger.warning('Successfully set p_pre from DCOPF results.')

    def def_sfr(self, du, dd):
        """
        Define the SFR requirements.

        Parameters
        ----------
        du : float
            SFR Up requirement.
        dd: float
            SFR Down requirement.
        """
        self.du = du
        self.dd = dd

    def data_check(self):
        """
        Check data consistency.
        """
        if not hasattr(self.cost, 'cru'):
            self.cost['cru'] = 0
            logger.warning('No RegUp cost data (``cru`` in ``cost``), set to 0.')
        if not hasattr(self.cost, 'crd'):
            self.cost['crd'] = 0
            logger.warning('No RegDn cost data(``crd`` in ``cost``), set to 0.')
        if not hasattr(self, 'du'):
            self.du = 0
            logger.warning('No RegUp requirement data (``du``), set to 0.')
        if not hasattr(self, 'dd'):
            self.dd = 0
            logger.warning('No RegDn requirement data (``dd``), set to 0.')
        if not hasattr(self.gen, 'ramp5'):
            self.gen['ramp5'] = 20
            logger.warning('No ramp AGC data (``ramp5`` in ``gen``), set to 20.')

    def build(self):
        """
        Build gurobipy model.

        Returns
        -------
        mdl : gurobipy.Model
            The gurobipy model.
        """
        self.data_check()

        # --- build RTED model ---
        self.update_dict()
        self.mdl = self._build_vars(self.mdl)
        self.mdl = self._build_obj(self.mdl)
        self.mdl = self._build_cons(self.mdl)
        logger.info('Successfully build RTED model.')

    def _build_vars(self, mdl):
        GEN = self.gendict.keys()
        # --- uncontrollable generators limit to p0 ---
        gencp = self.gen.copy()
        gencp['pmax'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        gencp['pmin'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        # --- offline geenrators limit to 0 ---
        gencp['pmax'][gencp.u == 0] = 0
        gencp['pmin'][gencp.u == 0] = 0
        # --- gen: pg ---
        self.pg = mdl.addVars(GEN, name='pg', vtype=gb.GRB.CONTINUOUS, obj=0,
                              ub=gencp.pmax.tolist(), lb=gencp.pmin.tolist())
        # --- RegUp, RegDn ---
        self.pru = mdl.addVars(GEN, name='pru', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=gencp.band.tolist(), lb=[0] * gencp.shape[0])
        self.prd = mdl.addVars(GEN, name='prd', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=gencp.band.tolist(), lb=[0] * gencp.shape[0])
        return mdl

    def _build_obj(self, mdl):
        GEN = self.gendict.keys()
        gendict = self.gendict
        costdict = self.costdict
        # --- minimize generation cost ---
        cost_pg = sum(self.pg[gen] * costdict[gen]['c1']
                      + self.pg[gen] * self.pg[gen] * costdict[gen]['c2']
                      + costdict[gen]['c0'] * gendict[gen]['u']  # online status
                      for gen in GEN)
        # --- RegUp, RegDn cost ---
        cost_ru = sum(self.pru[gen] * costdict[gen]['cru'] for gen in GEN)
        cost_rd = sum(self.pru[gen] * costdict[gen]['crd'] for gen in GEN)
        self.obj = mdl.setObjective(expr=cost_pg + cost_ru + cost_rd, sense=gb.GRB.MINIMIZE)
        return mdl

    def _build_cons(self, mdl):
        ptotal = self.load.p0.sum()

        gendict = self.gendict
        linedict = self.linedict
        gen_gsfdict = self.gen_gsfdict

        GEN = gendict.keys()
        LINE = linedict.keys()

        # --- power balance ---
        p_sum = sum(self.pg[gen] for gen in GEN)
        mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- line limits ---
        for line in LINE:
            lhs1 = sum(self.pg[gen] * gen_gsfdict[gen][line] for gen in GEN)
            mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')

        # --- GEN capacity ---
        mdl.addConstrs((self.pg[gen] + self.pru[gen] <= gendict[gen]['pmax'] for gen in GEN),
                       name='PG_max')
        mdl.addConstrs((self.pg[gen] - self.prd[gen] >= gendict[gen]['pmin'] for gen in GEN),
                       name='PG_mim')

        # --- SFR requirements ---
        # --- a) RegUp --
        mdl.addConstr(sum(self.pru[gen] for gen in GEN) == self.du, name='RegUp')
        # --- b) RegDn --
        mdl.addConstr(sum(self.prd[gen] for gen in GEN) == self.dd, name='RegDn')

        # --- ramp limits ---
        mdl.addConstrs((self.pg[gen] - gendict[gen]['p_pre'] <= gendict[gen]['ramp5']
                       for gen in GEN), name='RampU')
        mdl.addConstrs((self.pg[gen] - gendict[gen]['p_pre'] >= -1 * gendict[gen]['ramp5']
                       for gen in GEN), name='RampD')
        return mdl

    def get_res(self):
        """
        Get resutlts, can be used after mdl.optimize().

        Returns
        -------
        DataFrame
            The output DataFrame contains setpoints ``pg``

        """
        self.build()
        self.mdl.optimize()
        # --- check if mdl is sovled ---
        if not hasattr(self.pg[self.gen.idx[0]], 'X'):
            logger.warning('RTED has no valid resutls!')
            pg = [0] * self.gen.shape[0]
            pru = [0] * self.gen.shape[0]
            prd = [0] * self.gen.shape[0]
        else:
            logger.warning('Successfully solve RTED.')
            # --- gather data --
            pg = []
            pru = []
            prd = []
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
                pru.append(self.pru[gen].X)
                prd.append(self.prd[gen].X)
            # --- cost ---
            self.res_cost = self.mdl.getObjective().getValue()
            logger.info(f'Total cost={np.round(self.res_cost, 3)}')
        # --- build output table ---
        dcres = pd.DataFrame()
        dcres['gen'] = self.gen['idx']
        dcres['pg'] = pg
        dcres['pru'] = pru
        dcres['prd'] = prd
        dcres['bu'] = dcres['pru'] / dcres['pru'].sum()
        dcres['bd'] = dcres['prd'] / dcres['prd'].sum()
        dcres.fillna(0, inplace=True)
        return dcres


class rted2(rted):
    """
    RTED2 class, where type II generator is supported.

    Parameters
    ----------
    name : str
        Name of the system.

    Attributes
    ----------
    bus: pandas.DataFrame
        Bus data.
    gen: pandas.DataFrame
        Generator data.
    line: pandas.DataFrame
        Line data.
    load: pandas.DataFrame
        Load data.
    gen_gsf: pandas.DataFrame
        Generator shift factor of gen bus data.
    cost: pandas.DataFrame
        Cost data.
    """

    def __init__(self, name='rted'):
        super().__init__(name)
        # self.build()

    def data_check(self):
        super().data_check()
        if not hasattr(self.gen, 'type'):
            self.gen['type'] = '1'
        if not hasattr(self.gen, 'prumax'):
            self.gen['prumax'] = 0
        if not hasattr(self.gen, 'prdmax'):
            self.gen['prdmax'] = 0

    def def_type2(self, gen_idx, prumax, prdmax):
        """
        Define type 2 generator.

        Parameters
        ----------
        gen_idx : str
            Generator index that will be set to type 2.
        """
        self.gen['type'] = 1
        self.gen['prumax'] = 0
        self.gen['prdmax'] = 0
        for idx, rmax, dmax in zip(gen_idx, prumax, prdmax):
            row = self.gen[self.gen['idx'] == idx].index[0]
            self.gen['type'].iloc[row] = 2
            self.gen['prumax'].iloc[row] = rmax
            self.gen['prdmax'].iloc[row] = dmax

    def build(self):
        self.data_check()

        # --- build RTED model ---
        self.update_dict()
        self.mdl = gb.Model(self.name)
        self.mdl = self._build_vars(self.mdl)
        self.mdl = self._build_obj(self.mdl)
        self.mdl = self._build_cons(self.mdl)
        logger.info('Successfully build RTED2 model.')

    def _build_cons(self, mdl):
        ptotal = self.load.p0.sum()

        gendict = self.gendict
        linedict = self.linedict
        gen_gsfdict = self.gen_gsfdict

        GEN = gendict.keys()
        LINE = linedict.keys()

        # --- power balance ---
        p_sum = sum(self.pg[gen] for gen in GEN)
        mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- line limits ---
        for line in LINE:
            lhs1 = sum(self.pg[gen] * gen_gsfdict[gen][line] for gen in GEN)
            mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')

        # --- GEN capacity ---
        # --- filter Type II gen ---
        gendict_I = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 1:
                gendict_I[new_key] = new_value
        gendict_II = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 2:
                gendict_II[new_key] = new_value
        GENI = gendict_I.keys()
        GENII = gendict_II.keys()
        # --- a Type I GEN capacity limits ---
        mdl.addConstrs((self.pg[gen] + self.pru[gen] <= gendict[gen]['pmax'] for gen in GENI),
                       name='PG_max')
        mdl.addConstrs((self.pg[gen] - self.prd[gen] >= gendict[gen]['pmin'] for gen in GENI),
                       name='PG_mim')
        # --- b Type II Gen capacity and SFR limits---
        mdl.addConstrs((self.pg[gen] <= gendict[gen]['pmax'] for gen in GENII),
                       name='PG_max')
        mdl.addConstrs((self.pg[gen] >= gendict[gen]['pmin'] for gen in GENII),
                       name='PG_mim')
        mdl.addConstrs((self.pru[gen] <= gendict[gen]['prumax'] for gen in GENII),
                       name='PRU_max')
        mdl.addConstrs((self.prd[gen] <= gendict[gen]['prdmax'] for gen in GENII),
                       name='PRD_max')

        # --- SFR requirements ---
        # --- a) RegUp --
        mdl.addConstr(sum(self.pru[gen] for gen in GEN) == self.du, name='RegUp')
        # --- b) RegDn --
        mdl.addConstr(sum(self.prd[gen] for gen in GEN) == self.dd, name='RegDn')

        # --- ramp limits ---
        mdl.addConstrs((self.pg[gen] - gendict[gen]['p_pre'] <= gendict[gen]['ramp5']
                       for gen in GEN), name='RampU')
        mdl.addConstrs((self.pg[gen] - gendict[gen]['p_pre'] >= -1 * gendict[gen]['ramp5']
                       for gen in GEN), name='RampD')
        return mdl


class rted3(rted2):
    """
    RTED3 class, inherits from RTED2,
    where soft boud of SFR is applied.

    Punishment coeeficient for insufficient SFR and load
    are set to 1000 and 2000 by default.
    """
    def __init__(self, name='rted3'):
        super().__init__(name)
        self.ks = 1000
        self.kl = 2000

    def _build_obj(self, mdl):
        GEN = self.gendict.keys()
        gendict = self.gendict
        costdict = self.costdict
        # --- minimize generation cost ---
        cost_pg = sum(self.pg[gen] * costdict[gen]['c1']
                      + self.pg[gen] * self.pg[gen] * costdict[gen]['c2']
                      + costdict[gen]['c0'] * gendict[gen]['u']  # online status
                      for gen in GEN)
        # --- RegUp, RegDn cost ---
        cost_ru = sum(self.pru[gen] * costdict[gen]['cru'] for gen in GEN)
        cost_rd = sum(self.prd[gen] * costdict[gen]['crd'] for gen in GEN)
        soft_ru = self.du - sum(self.pru[gen] for gen in GEN)
        soft_rd = self.dd - sum(self.prd[gen] for gen in GEN)
        cost_sfr = cost_ru + cost_rd + self.ks * (soft_ru + soft_rd)
        # --- laod cost ---
        ptotal = self.load.p0.sum()
        pg = sum(self.pg[gen] for gen in GEN)
        cost_load = self.kl * (ptotal - pg)
        self.obj = mdl.setObjective(expr=cost_pg + cost_sfr + cost_load,
                                    sense=gb.GRB.MINIMIZE)
        return mdl

    def _build_cons(self, mdl):
        ptotal = self.load.p0.sum()

        gendict = self.gendict
        linedict = self.linedict
        gen_gsfdict = self.gen_gsfdict

        GEN = gendict.keys()
        LINE = linedict.keys()

        # --- SFR requirements ---
        # --- a) RegUp --
        mdl.addConstr(sum(self.pru[gen] for gen in GEN) <= self.du, name='RegUp')
        # --- b) RegDn --
        mdl.addConstr(sum(self.prd[gen] for gen in GEN) <= self.dd, name='RegDn')

        # --- power balance ---
        p_sum = sum(self.pg[gen] for gen in GEN)
        mdl.addConstr(p_sum <= ptotal, name='PowerBalance')

        # --- line limits ---
        for line in LINE:
            lhs1 = sum(self.pg[gen] * gen_gsfdict[gen][line] for gen in GEN)
            mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')

        # --- GEN capacity ---
        # --- filter Type II gen ---
        gendict_I = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 1:
                gendict_I[new_key] = new_value
        gendict_II = dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 2:
                gendict_II[new_key] = new_value
        GENI = gendict_I.keys()
        GENII = gendict_II.keys()
        # --- a Type I GEN capacity limits ---
        mdl.addConstrs((self.pg[gen] + self.pru[gen] <= gendict[gen]['pmax'] for gen in GENI),
                       name='PG_max')
        mdl.addConstrs((self.pg[gen] - self.prd[gen] >= gendict[gen]['pmin'] for gen in GENI),
                       name='PG_mim')
        # --- b Type II Gen capacity and SFR limits---
        mdl.addConstrs((self.pg[gen] <= gendict[gen]['pmax'] for gen in GENII),
                       name='PG_max')
        mdl.addConstrs((self.pg[gen] >= gendict[gen]['pmin'] for gen in GENII),
                       name='PG_mim')
        mdl.addConstrs((self.pru[gen] <= gendict[gen]['prumax'] for gen in GENII),
                       name='PRU_max')
        mdl.addConstrs((self.prd[gen] <= gendict[gen]['prdmax'] for gen in GENII),
                       name='PRD_max')

        # --- AGC ramp limits ---
        mdl.addConstrs((self.pg[gen] - gendict[gen]['p_pre'] <= gendict[gen]['ramp5']
                       for gen in GEN), name='RampU')
        mdl.addConstrs((gendict[gen]['p_pre'] - self.pg[gen] <= gendict[gen]['ramp5']
                       for gen in GEN), name='RampD')
        return mdl
