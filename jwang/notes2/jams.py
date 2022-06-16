from andes.interop.pandapower import to_pandapower, add_gencost
from andes.interop.pandapower import make_GSF, build_group_table
import gurobipy as gb
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)


class dcbase:
    """
    Base class of DC optimal power flow.

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
        """
        Base class of DC optimal power flow.

        Parameters
        ----------
        name : str
            Name of the system.
        """
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
        # --- build dict ---
        if not model:
            mdl_list = ['bus', 'gen', 'line', 'gen_gsf', 'cost']
        else:
            mdl_list = model
        for mdl in mdl_list:
            mdl_df = getattr(self, mdl)
            if mdl == 'line':
                sup = pd.DataFrame()
                sup['bus'] = self.bus['idx']
                sup = sup.merge(self.load[['bus', 'p0', 'sf']],
                                on='bus', how='left').fillna(0).rename(columns={'p0': 'load'})
                sup['net'] = (-1 * sup.load * sup.sf)
                sup2 = sup[['bus', 'net']].groupby('bus').sum()
                mdl_df['sup'] = np.matmul(self.gsf_matrix, sup2.net.values)
            mdl_df.index = mdl_df.idx
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
        All generators are set as controllable by default.
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
        # TODO: later on, merge 'ramp_5', 'ramp_10', 'ramp_30'
        self.gen['ramp_5'] = 100
        self.gen['ramp_10'] = 100
        self.gen['ramp_30'] = 100
        self.gen['sf'] = 1  # scaling factor
        # --- later on ---
        # self.gen['ramp_5'] = self.gen['ramp_5'] / self.mva
        # self.gen['ramp_10'] = self.gen['ramp_10'] / self.mva
        # self.gen['ramp_30'] = self.gen['ramp_30'] / self.mva
        # if self.gen['ramp_5'].max() == 0:
        #     self.gen['ramp_5'] = 100
        #     self.gen['ramp_10'] = 100
        #     self.gen['ramp_30'] = 100

        # --- load ---
        pq_cols = ['idx', 'u', 'name', 'bus', 'Vn', 'p0', 'q0',
                   'vmax', 'vmin', 'owner']
        self.load = ssa.PQ.as_df()[pq_cols].reset_index(drop=True)
        self.load['sf'] = 1  # scaling factor
        self.load.sort_index(inplace=True)

        # --- line ---
        line_cols = ['idx', 'u', 'name', 'bus1', 'bus2', 'Sn', 'fn', 'Vn1', 'Vn2',
                     'trans', 'tap', 'phi', 'rate_a', 'rate_b', 'rate_c']
        ssa_line = ssa.Line.as_df()
        self.line = ssa_line[line_cols][ssa_line['trans'] == 0].reset_index(drop=True)
        if self.line['rate_a'].max() == 0:
            self.line['rate_a'] = 2000 / self.mva
            self.line['rate_b'] = 2000 / self.mva
            self.line['rate_c'] = 2000 / self.mva

        # --- GSF ---
        ssp = to_pandapower(ssa)
        gsf_matrix = make_GSF(ssp)  # TODO: remove dependency on pandapower
        self.gsf_matrix = gsf_matrix
        gsfdata = pd.DataFrame(data=self.gsf_matrix,
                               columns=self.bus['idx'].values,
                               index=self.line['idx'].values)
        gsfT = gsfdata.T
        gsfT['bus'] = gsfT.index
        # gsfT['bus'] = self.bus['idx']
        self.gen_gsf = self.gen[['idx', 'name', 'bus']].merge(gsfT, on='bus', how='left')
        self.gen_gsf.index = self.gen_gsf['idx']

        # add power surplus, where the controlled gen is removed
        sup = pd.DataFrame()
        sup['bus'] = self.bus['idx']
        sup = sup.merge(self.load[['bus', 'p0']],
                        on='bus', how='left').fillna(0).rename(columns={'p0': 'load'})
        sup['net'] = (-1 * sup.load)
        sup2 = sup[['bus', 'net']].groupby('bus').sum()
        self.line['sup'] = np.matmul(gsf_matrix, sup2.net.values)

        # --- update dict ---
        self.data_check(skip_cost=True)
        self.update_dict(model=['bus', 'gen', 'line', 'load', 'gen_gsf'])

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

    def _data_complete(self, obj, attr, df=False):
        """
        Fill missing data.

        Parameters
        ----------
        obj : str
            Object name.
        attr : str
            Attribute name.
        df : bool, optional
            if obj is a pandas.DataFrame, set ``df=True``
        """
        # TODO: finish this
        if not hasattr(obj, attr):
            if df:
                df = getattr(obj, attr)

    def data_check(self, skip_cost=False, info=True):
        """
        Check data consistency:

        1. If scaling factors of gen and load are valid.
        1. If gen upper and lower limits are valid.
        1. If cost data exists, when set ``skip_cost=True``.

        Parameters
        ----------
        skip_cost : bool
            True to skip cost data check
        """
        warning_list = []
        absent_list = []
        if not self.gen.sf.between(left=0, right=1, inclusive='both').all():
            warning_list.append(['sf of gen'])
        if not self.load.sf.between(left=0, right=1, inclusive='both').all():
            warning_list.append(['sf of load'])
        equal = self.gen.pmin <= self.gen.pmax
        if not equal.all():
            warning_list.append(['pmax and pmin of gen'])

        if not skip_cost:
            if not hasattr(self, 'cost'):
                self._default_cost()
                absent_list.append('cost')
        if info:
            if len(warning_list) > 1:
                logger.warning(f'Suspected data: {warning_list}')
            if len(absent_list) > 1:
                logger.warning(f'Missing data: {absent_list}')
        return [warning_list, absent_list]


class dcopf(dcbase):
    """
    DC optimal power flow.

    In the modeling, power generation cost is minimized, s.t.:
    power balance, line limits, and generator active power limits.

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

    def __init__(self, name='DCOPF', OutputFlag=1):
        """
        DC Optimal Power Flow class.

        Parameters
        ----------
        name : str
            Name of the system.
        OutputFlag : int
            Output flag, 0 or 1.
        """
        super().__init__(name)
        self.OutputFlag = OutputFlag
        self.mdl = gb.Model(self.name)
        self.mdl.setParam('OutputFlag', OutputFlag)
        self.mdl_l_FLAG = False  # indicate if gb.model is loaded
        self.mdl_m_FLAG = False  # indicate if gb.model is modified

    def from_andes(self, ssa, no_build=False):
        super().from_andes(ssa)
        if not no_build: self.build(info=False)

    def build(self, info=True):
        """
        Build DCOPF mpdel as attribute ``mdl``, will call `update_dict()` first.

        After build, following attributes will be available:
        `mdl`: model
        `pg`: Vars, power generation; named as ``pg``, indexed by `gen.idx`
        `obj`: Obj function, power generation cost;
        `pb`: Constr, power balance; named as ``pb``
        `llu`: Constrs, line limit up; named as ``llu``, indexed by `line.idx`
        `lld`: Constrs, line limit down; named as ``lld``, indexed by `line.idx`
        """
        self.mdl = gb.Model(self.name)
        self.mdl.setParam('OutputFlag', self.OutputFlag)
        self.data_check(skip_cost=False, info=info)
        self.update_dict()
        # --- build DCOPF model ---
        self.build_vars()
        self.build_obj()
        self.build_cons()
        self.mdl_l_FLAG = True
        if info:
            logger.warning(f'{self.name} GB model is loaded.')

    def build_vars(self):
        """
        Build gurobi vars as attribtue ``pg``.
        """
        GEN = self.gendict.keys()
        # --- uncontrollable generators limit to p0 ---
        gencp = self.gen.copy()
        gencp['pmax'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        gencp['pmin'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        # --- offline geenrators limit to 0 ---
        gencp['pmax'][gencp.u == 0] = 0
        gencp['pmin'][gencp.u == 0] = 0
        # --- gen: pg ---
        self.pg = self.mdl.addVars(GEN, name='pg', vtype=gb.GRB.CONTINUOUS, obj=0,
                                   ub=gencp.pmax.tolist(), lb=gencp.pmin.tolist())
        return self.pg

    def build_obj(self):
        """
        Build gurobi objective as attribtue ``obj``.
        """
        GEN = self.gendict.keys()
        # --- minimize generation cost ---
        cost_pg = sum(self.pg[gen] * self.costdict[gen]['c1']
                      + self.pg[gen] * self.pg[gen] * self.costdict[gen]['c2']
                      + self.costdict[gen]['c0'] * self.gendict[gen]['u']  # online status
                      for gen in GEN)
        self.obj = self.mdl.setObjective(expr=cost_pg, sense=gb.GRB.MINIMIZE)
        return [self.obj, cost_pg]

    def build_cons(self):
        """
        Build gurobi constraints as attribtues ``pb``, ``llu```, and ``lld``.
        """
        ptotal = self.load.p0 * self.load.sf
        ptotal = np.sum(ptotal)

        GEN = self.gendict.keys()
        LINE = self.linedict.keys()

        # --- power balance ---
        p_sum = sum(self.pg[gen] * self.gendict[gen]['sf'] for gen in GEN)
        self.pb = self.mdl.addConstr(p_sum == ptotal, name='pb')

        # --- line limits ---
        lhs = {}
        for line in LINE:
            lhs[line] = sum(self.pg[gen] * self.gen_gsfdict[gen][line] for gen in GEN)

        self.llu = self.mdl.addConstrs((lhs[line]+self.linedict[line]['sup'] <= self.linedict[line]['rate_a']
                                        for line in LINE),
                                       name='llu')
        self.lld = self.mdl.addConstrs((lhs[line]+self.linedict[line]['sup'] >= -self.linedict[line]['rate_a']
                                        for line in LINE),
                                       name='lld')
        return self.mdl

    def solve(self, info=True, no_build=False):
        """
        Build and solve the model, will call "build()" first.

        Parameters
        ----------
        info: bool
            If True, will print out the solving information.
        no_build: bool
            If True, will not call ``build()`` and use existing ``mdl``.

        Returns
        -------
        res: DataFrame
            The output DataFrame contains setpoints ``pg``
        """

        pg = [0] * self.gen.shape[0]
        if not no_build:
            self.build(info=True)
        self.mdl.optimize()
        if self.mdl.Status == gb.GRB.OPTIMAL:
            if info:
                logger.warning(f'{self.name} is solved.')
            pg = []
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
            # --- cost ---
            total_cost = self.mdl.getObjective().getValue()
            logger.warning(f'{self.name}: total cost={np.round(total_cost, 3)}')
        else:
            if info:
                logger.warning(f'{self.name} solved to {self.mdl.Status}, please check.')
            pg = [0] * self.gen.shape[0]
        # --- build output table ---
        self.res = pd.DataFrame()
        self.res['gen'] = self.gen['idx'].values
        self.res['pg'] = pg
        self.res.fillna(0, inplace=True)
        return self.res

    def diagnostic(self):
        """
        Diagnostic of the model.
        """
        logger.warning(f'{self.name} diagnostic process:')

        self.data_check(skip_cost=False, info=True)

        if self.mdl_m_FLAG:
            self.mdl.solve()
            self.mdl_m_FLAG = False

        if self.mdl.Status == gb.GRB.INFEASIBLE:
            logger.warning(f'{self.name} is infeasible.')

        # --- load analysis ---
        pg = self.gen.pmax.sum()
        pl = self.load.p0.sum()
        ol = (pl-pg)/pg
        if pg < pl:
            logger.warning(f'Overload by {np.round(ol*100, 2)}%.')
            sf = np.floor(pg / pl * 100) / 100
            self.load.sf = sf
            self.update_dict()
            self.build(info=False)
            self.solve(info=False)
            self.mdl_m_FLAG = True
            if self.mdl.Status == gb.GRB.OPTIMAL:
                logger.warning(f'Scale load to {sf*100}% CAN address it.')
            else:
                logger.warning(f'Scale load to {sf*100}% CANNOT address it.')

            # reset load
            self.load.sf = 1
            self.update_dict()
            self.build(info=False)
            self.mdl.optimize()
            self.mdl_m_FLAG = False

        # --- IIS ---
        self.mdl.computeIIS()
        c = self.mdl.getConstrs()
        susp_c_idx = [i for i, e in enumerate(self.mdl.IISConstr) if e != 0]
        susp_c = []
        for i in susp_c_idx:
            susp_c.append(c[i])

        vars = self.mdl.getVars()
        susp_vl_idx = [i for i, e in enumerate(self.mdl.IISLB) if e != 0]
        susp_vl = []
        for i in susp_vl_idx:
            susp_vl.append(vars[i])

        susp_vu_idx = [i for i, e in enumerate(self.mdl.IISUB) if e != 0]
        susp_vu = []
        for i in susp_vu_idx:
            susp_vu.append(vars[i])

        if len(susp_c) > 0:
            logger.warning(f'{self.name} IISConstrs: {susp_c}')
        if len(susp_vl) > 0:
            logger.warning(f'{self.name} IISLB: {susp_vl}')
        if len(susp_vu) > 0:
            logger.warning(f'{self.name} IISUB: {susp_vu}')
        return True


class rted(dcopf):
    """
    Real-time economic dispatch (RTED) using DCOPF.

    In the modeling, cost of generation and SFR are minimized, s.t.:
    power balance, line limits, generator active power limits,
    and ramping limits.

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

    def __init__(self, name='RTED', OutputFlag=1):
        """
        Real-time economic dispatch (RTED) using DCOPF.

        Parameters
        ----------
        name : str
            Name of the system.
        """
        super().__init__(name=name, OutputFlag=OutputFlag)

    def from_andes(self, ssa, no_build=False):
        super().from_andes(ssa, no_build=True)
        self.gen['p_pre'] = 0
        if not no_build: self.build(info=False)

    def def_sfr(self, sfrur, sfrdr):
        """
        Define the SFR requirements as attribtues ``sfrur``, ``sfrdr``.

        Parameters
        ----------
        sfru : float
            SFR Up requirement.
        sfrd: float
            SFR Down requirement.
        """
        self.sfrur = sfrur
        self.sfrdr = sfrdr

    def data_check(self, skip_cost=False, info=True):
        """
        Check data consistency:

        1. If scaling factors of gen and load are valid.
        1. If gen upper and lower limits are valid.
        1. If cost data exists, when set ``skip_cost=True``.
        1. If cost data has cru, crd
        1. If SFR requirements data ``sfrur``, ``sfrdr`` exist.
        1. If gen data has ``ramp_5``.

        Parameters
        ----------
        skip_cost : bool
            True to skip cost data check
        """
        [warning_list, absent_list] = super().data_check(skip_cost=skip_cost, info=False)
        if not skip_cost:
            if not hasattr(self.cost, 'cru'):
                self.cost['cru'] = 0
                absent_list.append('cost.cru')
            if not hasattr(self.cost, 'crd'):
                self.cost['crd'] = 0
                absent_list.append('cost.crd')
            if not hasattr(self, 'sfrur'):
                self.sfrur = 0
                absent_list.append('cost.sfrur')
            if not hasattr(self, 'sfrdr'):
                self.sfrdr = 0
                absent_list.append('sfrdr')
            if not hasattr(self.gen, 'ramp_5'):
                self.gen['ramp_5'] = 100
                absent_list.append('ramp_5')
        if info:
            if len(warning_list) > 1:
                logger.warning(f'Suspected data: {warning_list}')
            if len(absent_list) > 1:
                logger.warning(f'Missing data: {absent_list}')
        return [warning_list, absent_list]

    def build(self, info=True):
        """
        Build RTED model as the attribute ``mdl``, will call `update_dict()` first.

        After build, following attributes will be available:
        `mdl`: model
        `pg`: Vars, power generation; named as ``pg``, indexed by `gen.idx`
        `pru`: Vars, RegUp power; named as ``pru``, indexed by `gen.idx`
        `prd``: Vars, RegUp power; named as ``prd``, indexed by `gen.idx`
        `obj`: Obj function, power generation cost;
        `pb`: Constr, power balance; named as ``pb``
        `llu`: Constrs, line limit up; named as ``llu``, indexed by `line.idx`
        `lld`: Constrs, line limit down; named as ``lld``, indexed by `line.idx`
        `pgmax`: Constrs, generator limit up; named as ``pgmax``, indexed by `gen.idx`
        `pgmin`: Constrs, generator limit down; named as ``pgmin``, indexed by `gen.idx`
        `sfru`: Constr, SFR up; named as ``sfru``
        `sfrd`: Constr, SFR down; named as ``sfrd``
        `rampu`: Constrs, ramping limit up; named as ``rampu``, indexed by `gen.idx`
        `rampd`: Constrs, ramping limit down; named as ``rampd``, indexed by `gen.idx`
        """
        super().build(info=info)

    def build_vars(self):
        """
        Build gurobi vars as attribtues ``pg``, ``pru``, and ``prd``.
        """
        super().build_vars()
        GEN = self.gendict.keys()
        # --- RegUp, RegDn ---
        self.pru = self.mdl.addVars(GEN, name='pru', vtype=gb.GRB.CONTINUOUS, obj=0)
        self.prd = self.mdl.addVars(GEN, name='prd', vtype=gb.GRB.CONTINUOUS, obj=0)
        return [self.pg, self.pru, self.prd]

    def build_obj(self):
        """
        Build gurobi objective as attribtue ``obj``.
        """
        GEN = self.gendict.keys()
        [_, cost_pg] = super().build_obj()
        # --- RegUp, RegDn cost ---
        cost_ru = sum(self.pru[gen] * self.costdict[gen]['cru'] for gen in GEN)
        cost_rd = sum(self.pru[gen] * self.costdict[gen]['crd'] for gen in GEN)
        self.obj = self.mdl.setObjective(expr=cost_pg + cost_ru + cost_rd,
                                         sense=gb.GRB.MINIMIZE)
        return [self.obj, cost_pg, cost_ru, cost_rd]

    def build_cons(self):
        """
        Build gurobi constraints as attribtues ``pb``, ``llu``, ``lld``, ``pgmax``,
        ``pgmin``, ``sfru``, ``sfrd``, ``rampu``, and ``rampd``.
        """
        super().build_cons()
        GEN = self.gendict.keys()

        # --- GEN capacity ---
        self.pgmax = self.mdl.addConstrs((self.pg[gen] + self.pru[gen] <= self.gendict[gen]['pmax']
                                          for gen in GEN),
                                         name='pgmax')
        self.pgmin = self.mdl.addConstrs((self.pg[gen] - self.prd[gen] >= self.gendict[gen]['pmin']
                                          for gen in GEN),
                                         name='pgmin')

        # --- SFR requirements ---
        # --- a) RegUp --
        self.sfru = self.mdl.addConstr(sum(self.pru[gen] for gen in GEN) == self.sfrur, name='sfru')
        # --- b) RegDn --
        self.sfrd = self.mdl.addConstr(sum(self.prd[gen] for gen in GEN) == self.sfrdr, name='sfrd')

        # --- ramp limits ---
        self.rampu = self.mdl.addConstrs((self.pg[gen] - self.gendict[gen]['p_pre'] <= self.gendict[gen]['ramp_5']
                                          for gen in GEN), name='rampu')
        self.rampd = self.mdl.addConstrs((self.gendict[gen]['p_pre'] - self.pg[gen] <= self.gendict[gen]['ramp_5']
                                          for gen in GEN), name='rampd')
        return self.mdl

    def solve(self, info=True, no_build=False,
              disable_sfr=False, disable_ramp=False):
        """
        Build and solve the model, will call "build()" first.

        Parameters
        ----------
        info: bool
            If True, will print out the solving information.
        no_build: bool
            If True, will not call ``build()`` and use existing ``mdl``.

        Returns
        -------
        res: DataFrame
            The output DataFrame contains setpoints ``pg``,
            RegUp power ``pru``, RegDn power ``prd``,
            RegUp factor ``bu``, and RegDn factor ``bd``.
        """

        pg = [0] * self.gen.shape[0]
        pru = [0] * self.gen.shape[0]
        prd = [0] * self.gen.shape[0]
        if not no_build:
            self.build(info=True)
        if disable_ramp:
            disable_list = [self.rampu, self.rampd]
            rl = ['rampu', 'rampd']
        if disable_sfr:
            disable_list = [self.pgmax, self.pgmin, self.sfru, self.sfrd,
                            self.rampu, self.rampd]
            rl = ['pgmax', 'pgmin', 'sfru', 'sfrd', 'rampu', 'rampd']
        if disable_sfr or disable_ramp:
            for c in disable_list:
                self.mdl.remove(c)
            logger.warning(f'{self.name} removed Constrs: {rl}')
        self.mdl.optimize()
        if self.mdl.Status == gb.GRB.OPTIMAL:
            if info:
                logger.warning(f'{self.name} is solved.')
            pg = []
            pru = []
            prd = []
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
                pru.append(self.pru[gen].X)
                prd.append(self.prd[gen].X)
            # --- cost ---
            total_cost = self.mdl.getObjective().getValue()
            logger.warning(f'{self.name}: total cost={np.round(total_cost, 3)}')
        else:
            if info:
                logger.warning('Optimization ended with status %d' % self.mdl.Status)
            pg = [0] * self.gen.shape[0]
            pru = [0] * self.gen.shape[0]
            prd = [0] * self.gen.shape[0]
        # --- build output table ---
        self.res = pd.DataFrame()
        self.res['gen'] = self.gen['idx'].values
        self.res['pg'] = pg
        self.res['pru'] = pru
        self.res['prd'] = prd
        self.res['bu'] = self.res['pru'] / self.res['pru'].sum()
        self.res['bd'] = self.res['prd'] / self.res['prd'].sum()
        self.res.fillna(0, inplace=True)
        return self.res


class rted2(rted):
    """
    Real-time economic dispatch (RTED) using DCOPF,
    where type2 generator is supported.

    In the modeling, cost of generation and SFR are minimized, s.t.:
    power balance, line limits, generator active power limits,
    and ramping limits.

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

    def __init__(self, name='RTED2', OutputFlag=1):
        """
        Real-time economic dispatch (RTED) using DCOPF,
        where type2 generator is supported.

        Parameters
        ----------
        name : str
            Name of the system.
        """
        super().__init__(name=name, OutputFlag=OutputFlag)

    def from_andes(self, ssa):
        super().from_andes(ssa)
        self.gen['type'] = '1'
        self.gen['prumax'] = 0
        self.gen['prdmax'] = 0
        self.update_dict(model=['gen'])

    def data_check(self, skip_cost=False, info=True):
        """
        Check data consistency:

        1. If scaling factors of gen and load are valid.
        1. If gen upper and lower limits are valid.
        1. If cost data exists, when set ``skip_cost=True``.
        1. If cost data has cru, crd
        1. If SFR requirements data ``sfrur``, ``sfrdr`` exist.
        1. If gen data has ``ramp_5``.
        1. If gen data has ``type``, ``prumax``, ``prdmax``.

        Parameters
        ----------
        skip_cost : bool
            True to skip cost data check
        """
        [warning_list, absent_list] = super().data_check(skip_cost=skip_cost, info=False)
        if not skip_cost:
            if not hasattr(self.gen, 'type'):
                self.gen['type'] = '1'
                absent_list.append('gen.type')
            if not hasattr(self.gen, 'prumax'):
                self.gen['prumax'] = 0
                absent_list.append('gen.prumax')
            if not hasattr(self.gen, 'prdmax'):
                self.gen['prdmax'] = 0
                absent_list.append('gen.prdmax')
        if info:
            if len(warning_list) > 1:
                logger.warning(f'Suspected data: {warning_list}')
            if len(absent_list) > 1:
                logger.warning(f'Missing data: {absent_list}')
        return [warning_list, absent_list]

    def def_type2(self, gen_idx, prumax, prdmax):
        """
        Define type 2 generator.

        Length of ``gen_idx`` must be equal to length of ``prumax`` and ``prdmax``.

        Parameters
        ----------
        gen_idx : list of str
            Generator index that will be set to type 2.
        """
        for n in range(len(gen_idx)):
            row = self.gen[self.gen['idx'] == gen_idx[n]].index[0]
            self.gen['type'].loc[row] = 2
            self.gen['prumax'].loc[row] = prumax[n]
            self.gen['prdmax'].loc[row] = prdmax[n]

    def build(self, info=True):
        """
        Build RTED model as the attribute ``mdl``, will call `update_dict()` first.

        After build, following attributes will be available:
        `mdl`: model
        `pg`: Vars, power generation; named as ``pg``, indexed by `gen.idx`
        `pru`: Vars, RegUp power; named as ``pru``, indexed by `gen.idx`
        `prd``: Vars, RegUp power; named as ``prd``, indexed by `gen.idx`
        `obj`: Obj function, power generation cost;
        `pb`: Constr, power balance; named as ``pb``
        `llu`: Constrs, line limit up; named as ``llu``, indexed by `line.idx`
        `lld`: Constrs, line limit down; named as ``lld``, indexed by `line.idx`
        `pgmax`: Constrs, generator limit up; named as ``pgmax``, indexed by `gen.idx`
        `pgmin`: Constrs, generator limit down; named as ``pgmin``, indexed by `gen.idx`
        `prumax`: Constrs, SFR upper limit; named as ``prumax``, indexed by `gen.idx` of type2 generator
        `prdmax`: Constrs, SFR lower limit; named as ``prdmax``, indexed by `gen.idx` of type2 generator
        `sfru`: Constr, SFR up; named as ``sfru``
        `sfrd`: Constr, SFR down; named as ``sfrd``
        `rampu`: Constrs, ramping limit up; named as ``rampu``, indexed by `gen.idx`
        `rampd`: Constrs, ramping limit down; named as ``rampd``, indexed by `gen.idx`
        """
        super().build(info=info)

    def build_cons(self):
        super().build_cons()
        self.mdl.remove(self.pgmax)
        self.mdl.remove(self.pgmin)

        # --- GEN capacity ---
        # --- filter Type II gen ---
        gendict_I = dict()
        gendict_II = dict()
        for (new_key, new_value) in self.gendict.items():
            if int(new_value['type']) == 1:
                gendict_I[new_key] = new_value
            elif int(new_value['type']) == 2:
                gendict_II[new_key] = new_value
        GENI = gendict_I.keys()
        GENII = gendict_II.keys()
        # --- a Type I GEN capacity limits ---
        pgmaxI = self.mdl.addConstrs((self.pg[gen] + self.pru[gen] <= self.gendict[gen]['pmax'] for gen in GENI),
                                     name='pgmax')
        pgminI = self.mdl.addConstrs((self.pg[gen] - self.prd[gen] >= self.gendict[gen]['pmin'] for gen in GENI),
                                     name='pgmin')
        # --- b Type II Gen capacity and SFR limits---
        pgmaxII = self.mdl.addConstrs((self.pg[gen] <= self.gendict[gen]['pmax'] for gen in GENII),
                                      name='pgmax')
        pgminII = self.mdl.addConstrs((self.pg[gen] >= self.gendict[gen]['pmin'] for gen in GENII),
                                      name='pgmin')
        self.prumax = self.mdl.addConstrs((self.pru[gen] <= self.gendict[gen]['prumax'] for gen in GENII),
                                          name='prumax')
        self.prdmax = self.mdl.addConstrs((self.prd[gen] >= -1 * self.gendict[gen]['prdmax'] for gen in GENII),
                                          name='prdmax')
        self.pgmax = pgmaxI | pgmaxII
        self.pgmin = pgminI | pgminII
        return self.mdl
