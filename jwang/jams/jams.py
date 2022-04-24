from andes.interop.pandapower import to_pandapower, add_gencost
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
        Generator-Shunt data.
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
            mdl_df = getattr(self, mdl)
            mdl_df.index = mdl_df.index
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
        # TODO: later on, merge 'ramp_agc', 'ramp_10', 'ramp_30'
        self.gen['ramp_agc'] = 100
        self.gen['ramp_10'] = 100
        self.gen['ramp_30'] = 100
        # --- later on ---
        # self.gen['ramp_agc'] = self.gen['ramp_agc'] / self.mva
        # self.gen['ramp_10'] = self.gen['ramp_10'] / self.mva
        # self.gen['ramp_30'] = self.gen['ramp_30'] / self.mva
        # if self.gen['ramp_agc'].max() == 0:
        #     self.gen['ramp_agc'] = 100
        #     self.gen['ramp_10'] = 100
        #     self.gen['ramp_30'] = 100

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
        gsfdata['line'] = self.line['idx']
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
        self.update_dict()

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
    """

    def __init__(self, name='dcopf'):
        super().__init__(name)
        # self.build()

    def from_andes(self, ssa):
        super().from_andes(ssa)

    def build(self):
        self.update_dict()
        self.mdl = gb.Model(self.name)

        # --- build GB model ---
        self.mdl = self._build_vars(self.mdl)
        self.mdl = self._build_obj(self.mdl)
        self.mdl = self._build_cons(self.mdl)
        logger.info('Successfully build DCOPF model.')

    def _build_vars(self, mdl):
        GEN = self.gendict.keys()
        # --- uncontrollable generators limit to p0 ---
        gencp = self.gen.copy()
        gencp['pmax'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        gencp['pmin'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        # --- offline geenrators limit to 0 ---
        gencp['pmax'][gencp.u == 0] = 0
        gencp['pmin'][gencp.u == 0] = 0
        # --- gen: p_sch ---
        self.p_sch = mdl.addVars(GEN, name='p_sch', vtype=gb.GRB.CONTINUOUS, obj=0,
                                 ub=gencp.pmax.tolist(), lb=gencp.pmin.tolist())
        return mdl

    def _build_obj(self, mdl):
        GEN = self.gendict.keys()
        gendict = self.gendict
        costdict = self.costdict
        # --- minimize generation cost ---
        cost_pg = sum(self.p_sch[gen] * costdict[gen]['c1']
                      + self.p_sch[gen] * self.p_sch[gen] * costdict[gen]['c2']
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
        p_sum = sum(self.p_sch[gen] for gen in GEN)
        mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- line limits ---
        for line in LINE:
            lhs1 = sum(self.p_sch[gen] * gen_gsfdict[gen][line] for gen in GEN)
            # mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            # mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')
        return mdl

    def get_res(self):
        """
        Get resutlts, can be used after mdl.optimize().

        Returns
        -------
        DataFrame
            The output DataFrame contains setpoints ``p_sch``
        """

        # --- check if mdl is sovled ---
        if not hasattr(self.p_sch[self.gen.idx[0]], 'X'):
            logger.warning('DCOPF has no valid resutls!')
            p_sch = [0] * self.gen.shape[0]
        else:
            # --- gather data --
            p_sch = []
            for gen in self.gendict.keys():
                p_sch.append(self.p_sch[gen].X)
            # --- cost ---
            total_cost = self.mdl.getObjective().getValue()
            logger.info(f'Total cost={np.round(total_cost, 3)}')
        # --- build output table ---
        dcres = pd.DataFrame()
        dcres['gen'] = self.gen['idx']
        dcres['p_sch'] = p_sch
        dcres.fillna(0, inplace=True)
        return dcres
