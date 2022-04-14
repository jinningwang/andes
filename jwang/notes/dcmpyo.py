import pyomo.environ as pyo
import pandas as pd
import numpy as np
from makeGSF import makeGSF_ppn
import logging
logger = logging.getLogger(__name__)


class dcm:
    """
    Base class for DCOPF with SFR model in Pyomo.

    Attributes:

    - `mdl`: the DCOPF with SFR model.

    Notes
    -----
      - Input `ssp.load` should be the expected load
      - Input `rampu` and `rampd` should be converted to 5-min-based
      - Attrbute `ptotal` is calculated from `ssp.load.p_mw` insidely and converted to p.u.
    """

    def __init__(self, ssp, len_ev):
        # --- store original info of ssp ---
        # --- 2. gen table ---
        self.gen_ssp = ssp.gen[['bus', 'p_mw']]

        # --- def. length ---
        self.len_cgen = ssp.gen[ssp.gen.controllable].shape[0]
        self.len_ev = len_ev

        # --- def. data ---
        # --- bus ---
        # --- load ---
        self.loaddata = ssp.load[['bus', 'p_mw']]
        self.loaddata['p'] = self.loaddata.p_mw / ssp.sn_mva
        self.loaddict = self.loaddata[['bus', 'p']].T.to_dict()

        # --- 0. gsf ---
        busdata = ssp.bus[['name']]
        busdata = self.set_idx(busdata, 'bus')
        gsf = makeGSF_ppn(ssp)
        self.gsfdata = pd.DataFrame(gsf)
        self.gsfdata = self.set_idx(self.gsfdata, 'line')
        # self.gsfdata.drop(['name'], axis=1, inplace=True)

        # --- 1. gen ---
        ssp_gen = \
            pd.merge(left=ssp.poly_cost.rename(columns={'element': 'index'}),
                     right=ssp.gen.reset_index(),
                     on='index', how='left')
        self.gendata = pd.DataFrame()
        # --- a. controllable gen cost, linear ---
        self.gendata['cg'] = ssp_gen.cp1_eur_per_mw[ssp_gen.controllable].tolist()
        # --- b. controllable gen bus ---
        self.gendata['bus'] = ssp_gen.bus[ssp_gen.controllable].tolist()
        # --- c. controllable gen limit ---
        self.gendata['pmax'] = list(ssp_gen.max_p_mw[ssp_gen.controllable] / ssp.sn_mva)
        self.gendata['pmin'] = list(ssp_gen.min_p_mw[ssp_gen.controllable] / ssp.sn_mva)

        self.ps_defined = False
        # TODO: SFR cost
        self.gendata['cg_ru'] = [0] * self.len_cgen
        self.gendata['cg_rd'] = [0] * self.len_cgen
        # TODO: ramp
        self.gendata['rampu'] = [10] * self.len_cgen
        self.gendata['rampd'] = [10] * self.len_cgen
        # --- set index ---
        self.gendata = self.set_idx(self.gendata, 'cgen')
        # convert to a dict
        self.gendict = self.gendata.T.to_dict()

        # --- 1. sup gen ---
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
        self.linedata['sup'] = np.matmul(gsf, sup.net.values)
        # --- set index ---
        self.linedata = self.set_idx(self.linedata, 'line')
        self.linedict = self.linedata.T.to_dict()

        # --- 2. ev ---
        self.evdata = pd.DataFrame()
        self.evdata['ce_ru'] = [0] * self.len_ev
        self.evdata['ce_rd'] = [0] * self.len_ev

        # TODO: replace pr_eu, pr_ed with computed value
        self.evpr_defined = False
        self.evdata['pr_eu_max'] = [0] * self.len_ev
        self.evdata['pr_ed_max'] = [0] * self.len_ev
        # --- set index ---
        self.evdata = self.set_idx(self.evdata, 'ev')
        self.evdict = self.evdata.T.to_dict()

        # --- 3. orcasted load up and down ---
        self.ptotal = ssp.load.p_mw.sum() / ssp.sn_mva
        self.dpd_defined = False
        self.dpd_u = 0.2
        self.dpd_d = 0.2

        # --- build model ---
        # self.build_crt()

    def set_idx(self, df, name='device'):
        """Set index of the DataFrame."""
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

    def defvar(self, dpd_u=0.2, dpd_d=0.2, pr_eu_max=[], pr_ed_max=[]):
        """Define model parameters (p.u.): dpd_u, dpd_d, pr_eu_max, pr_ed_max"""
        # --- load, flucation ---
        self.dpd_u = dpd_u
        self.dpd_d = dpd_d

        # --- ev SFR capacity ---
        self.set_list(pr_eu_max, 'evdata', 'pr_eu_max')
        self.set_list(pr_ed_max, 'evdata', 'pr_ed_max')
        self.evdict = self.evdata.T.to_dict()

        self.dpd_defined = True
        self.evpr_defined = True
        return self.dpd_defined & self.evpr_defined

    def defps(self, cg_ru=[], cg_rd=[], rampu=[], rampd=[], ce_ru=[], ce_rd=[]):
        """Define power system parameters: cg_ru, cg_rd, rampu, rampd, ce_ru, ce_rd"""
        # --- gen SFR cost, ramp---
        self.set_list(cg_ru, 'gendata', 'cg_ru')
        self.set_list(cg_rd, 'gendata', 'cg_rd')
        self.set_list(rampu, 'gendata', 'rampu')
        self.set_list(rampd, 'gendata', 'rampd')
        self.gendict = self.gendata.T.to_dict()

        # --- ev SFR cost---
        self.set_list(ce_ru, 'evdata', 'ce_ru')
        self.set_list(ce_rd, 'evdata', 'ce_rd')
        self.evdict = self.evdata.T.to_dict()

        self.ps_defined = True
        return self.ps_defined

    def _build_vars_dc(self, mdl, GEN):
        # --- 1. gen: p_sch  ---
        mdl.p_sch = pyo.Var(GEN, domain=pyo.Reals)
        return mdl

    def _build_obj_dc(self, mdl, gendict, GEN):
        cost_pg = sum(mdl.p_sch[gen] * gendict[gen]['cg'] for gen in GEN)
        mdl.obj = pyo.Objective(expr=cost_pg, sense=pyo.minimize)
        return mdl

    def _build_cons_dc(self, mdl, gen2dict, linedict, GEN, LINE, ptotal):
        mdl.cons = pyo.ConstraintList()

        # --- 4. power balance ---
        p_sum = sum(mdl.p_sch[gen] for gen in GEN)
        #   + sum(self.mdl.dp_g[gen] for gen in GEN) \
        #   + sum(self.mdl.dp_e[ev] for ev in EV)
        # TODO: do we need consider the load flucation here? I don't think so
        mdl.cons.add(expr=p_sum == ptotal)

        # --- 5. line limits ---
        for line in LINE:
            lhs1 = sum(mdl.p_sch[gen] * gen2dict[gen][line] for gen in GEN)
            mdl.cons.add(expr=lhs1+linedict[line]['sup'] <= linedict[line]['lim'])

        return mdl

    def _build_vars(self, mdl, GEN, EV):
        # --- 1. gen: p_sch, pr_gu, pr_gd, dp_g, b_gu, b_gd ---
        mdl.p_sch = pyo.Var(GEN, domain=pyo.Reals)
        mdl.pr_gu = pyo.Var(GEN, domain=pyo.NonNegativeReals)
        mdl.pr_gd = pyo.Var(GEN, domain=pyo.NonNegativeReals)
        # self.mdl.dp_g = pyo.Var(GEN, domain=pyo.Reals)
        mdl.b_gu = pyo.Var(GEN, domain=pyo.Reals, bounds=(0, 1))
        mdl.b_gd = pyo.Var(GEN, domain=pyo.Reals, bounds=(0, 1))
        # --- 2. ev: dp_e, b_eu, b_ed ---
        mdl.pr_eu = pyo.Var(EV, domain=pyo.NonNegativeReals)
        mdl.pr_ed = pyo.Var(EV, domain=pyo.NonNegativeReals)
        # self.mdl.dp_e = pyo.Var(EV, domain=pyo.Reals)
        mdl.b_eu = pyo.Var(EV, domain=pyo.Reals, bounds=(0, 1))
        mdl.b_ed = pyo.Var(EV, domain=pyo.Reals, bounds=(0, 1))
        return mdl

    def _build_obj(self, mdl, gendict, evdict, GEN, EV):
        cost_pg = sum(mdl.p_sch[gen] * gendict[gen]['cg'] for gen in GEN)
        cost_pru = sum(mdl.pr_gu[gen] * gendict[gen]['cg_ru'] for gen in GEN)
        cost_prd = sum(mdl.pr_gd[gen] * gendict[gen]['cg_rd'] for gen in GEN)
        # TODO: We don't have to include EV in the obj.?
        cost_eru = sum(mdl.pr_eu[ev] * evdict[ev]['ce_ru'] for ev in EV)
        cost_erd = sum(mdl.pr_ed[ev] * evdict[ev]['ce_rd'] for ev in EV)
        mdl.obj = pyo.Objective(expr=cost_pg + cost_pru + cost_prd + cost_eru + cost_erd, sense=pyo.minimize)
        return mdl

    def _build_cons(self, mdl, gendict, gen2dict, evdict, linedict, GEN, EV, LINE, dpd_u, dpd_d, ptotal):
        mdl.cons = pyo.ConstraintList()
        # --- 1. SFR capacity ---
        for gen in GEN:
            mdl.cons.add(mdl.p_sch[gen] + mdl.pr_gu[gen] <= gendict[gen]['pmax'])
            mdl.cons.add(mdl.p_sch[gen] - mdl.pr_gd[gen] >= gendict[gen]['pmin'])
        for ev in EV:
            mdl.cons.add(mdl.pr_eu[ev] <= evdict[ev]['pr_eu_max'])
            mdl.cons.add(mdl.pr_ed[ev] <= evdict[ev]['pr_ed_max'])

        # --- 2. SFR requirements ---
        # --- a) RegUp --
        dp_gu_sum = sum(mdl.b_gu[gen] * mdl.pr_gu[gen] for gen in GEN)
        dp_eu_sum = sum(mdl.b_eu[ev] * mdl.pr_eu[ev] for ev in EV)
        mdl.cons.add(expr=dp_gu_sum + dp_eu_sum >= dpd_u)
        # --- b) RegDn --
        dp_gd_sum = sum(mdl.b_gd[gen] * mdl.pr_gd[gen] for gen in GEN)
        dp_ed_sum = sum(mdl.b_ed[ev] * mdl.pr_ed[ev] for ev in EV)
        mdl.cons.add(expr=dp_gd_sum + dp_ed_sum >= dpd_d)

        # --- 3. beta equality ---
        # --- a) RegUp --
        b_gu_sum = sum(mdl.b_gu[gen] for gen in GEN)
        b_eu_sum = sum(mdl.b_eu[ev] for ev in EV)
        mdl.cons.add(expr=b_gu_sum + b_eu_sum == 1)
        # --- b) RegDn --
        b_gd_sum = sum(mdl.b_gd[gen] for gen in GEN)
        b_ed_sum = sum(mdl.b_ed[ev] for ev in EV)
        mdl.cons.add(expr=b_gd_sum + b_ed_sum == 1)

        # --- 4. power balance ---
        p_sum = sum(mdl.p_sch[gen] for gen in GEN)
        #   + sum(self.mdl.dp_g[gen] for gen in GEN) \
        #   + sum(self.mdl.dp_e[ev] for ev in EV)
        # TODO: do we need consider the load flucation here? I don't think so
        mdl.cons.add(expr=p_sum == ptotal)

        # --- 5. line limits ---
        for line in LINE:
            lhs1 = sum(mdl.p_sch[gen] * gen2dict[gen][line] for gen in GEN)
            mdl.cons.add(expr=lhs1+linedict[line]['sup'] <= linedict[line]['lim'])

        # TODO: ramp limits
        # --- 6. ramp limits ---

        return mdl

    def build(self):
        """Build the DCOPF with SFR model"""

        if not (self.dpd_defined & self.evpr_defined):
            logger.warning('Computed vars `dpd_u`, `dpd_d`, `pr_eu`, `pr_ed` are not fully given.'
                           'Default values are used.' + '\n'
                           'They can be defined by defvar(dpd_u=0.2, dpd_d=0.2, pr_eu=[], pr_ed=[])')
        if not self.ps_defined:
            logger.warning('PS vars `cg_ru`, `cg_rd`, `rampu`, `rampd`, `cg_eu`, `cg_ed` are not given.'
                           'Default values are used.' + '\n'
                           'They can be defined by defps(cg_ru=[], cg_rd=[], rampu=[], rampd=[], cg_eu=[],' + '\n'
                           'cg_ed=[])')

        mdl = pyo.ConcreteModel()

        GEN = self.gendict.keys()
        EV = self.evdict.keys()
        LINE = self.linedict.keys()

        # --- A. decision vars ---
        mdl = self._build_vars(mdl, GEN, EV)

        # --- B. obj. ---
        mdl = self._build_obj(mdl, self.gendict, self.evdict, GEN, EV)

        # --- C. constraints ---
        mdl = self._build_cons(mdl, self.gendict, self.gen2dict,
                               self.evdict, self.linedict, GEN, EV, LINE,
                               self.dpd_u, self.dpd_d, self.ptotal)
        self.mdl = mdl
        return True

    def builddc(self):
        """"Build the DCOPF without SFR model"""

        mdldc = pyo.ConcreteModel()
        GEN = self.gendict.keys()
        LINE = self.linedict.keys()

        # --- A. decision vars ---
        mdldc = self._build_vars_dc(mdldc, GEN)
        # --- B. obj. ---
        mdldc = self._build_obj_dc(mdldc, self.gendict, GEN)
        # --- C. constraints ---
        mdldc = self._build_cons_dc(mdldc, self.gen2dict, self.linedict, GEN, LINE, self.ptotal)

        self.mdldc = mdldc
        return True
