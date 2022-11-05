"""EV Aggregator based on State Space Modeling"""

import itertools
from tqdm import tqdm
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import OrderedDict
# import swifter
# from swifter import set_defaults
# set_defaults(
#     npartitions=None,
#     dask_threshold=1,
#     scheduler='proceselfs',
#     progress_bar=False,
#     progress_bar_desc=None,
#     allow_dask_on_strings=False,
#     force_parallel=True,
# )

import logging
logger = logging.getLogger(__name__)


# --- Functions ---
def update_xl(inl_input):
    '''
    Update online staus records, x series.

    Parameters:
    inl_input: dict
        Single EV information, ['sx', 'xl', 'u', 'u0', 'c2', 'ts', 'c0', 'bd']
    '''
    [sx, xl0, u, u0, c2, ts, c0, bd] = inl_input.copy()
    xl = xl0.copy()
    case1 = (u*u0 == 1) & (c2 == c0)  # --- cont. online & same c ---
    case2 = (u*u0 == 1) & (c2 != c0) & (bd == 1)  # --- cont. online & change c by boundary ---
    case3 = (u*u0 == 1) & (c2 != c0) & (bd == 0)  # --- cont. online & change c by control ---
    case4 = (1-u0)*u == 1  # offline -> online

    if case1 | case2:
        if len(xl[0]) == 0:
            xl[0] = [[c2]]
            xl[1] = [[sx]]
            xl[2] = [[ts]]
        else:
            xl[0][-1].append(c2)
            xl[1][-1].append(sx)
            xl[2][-1].append(ts)
    elif case3 | case4:
        if len(xl[0]) == 0:
            xl[0] = [[c2]]
            xl[1] = [[sx]]
            xl[2] = [[ts]]
        else:
            xl[0].append([c2])
            xl[1].append([sx])
            xl[2].append([ts])
    return xl


def safe_div(x, y):
    '''
    Safe division, return 0 if y is 0.

    Parameters
    ----------
    x: float
        numerator.
    y: float
        denominator.
    '''
    if y == 0:
        return 0
    else:
        return x/y


# --- EV State Space Model---
class ev_ssm():
    '''
    Class of EV State Space Model.

    EV parameters:
    Pc, Pd, nc, nd, Q follow uniform distribution.
    soci, socd, ts, tf follows normal distribution.

    The EV parameters are defined as two types: parameters following
    uniform distribution are stored in ``ev_ufparam``, while the parameters
    following normal distribution are stored in ``ev_nfparam``.

    Attributes
    ----------
    xtab: pandas.DataFrame
        EV states table, only online EVs are counted.
    ne: int
        Number of online EVs
    N: int
        Number of total EVs
    Np: int
        SSM update cycle.
    ev: pandas.DataFrame
        EV dataset
        u: online status
        u0: previous online status
        sx: SOC interval
        xl: list of lists, [[states], [sx], [t]]
    ts: float
        current time (in 24H)

    Notes
    -----
    ev_ufparam:
        Pl: rated charging/discharging power (kW) lower bound
        Pu: rated charging/discharging power (kW) upper bound
        nl: charging/discharging efficiency lower bound
        nu: charging/discharging efficiency upper bound
        Ql: Battery capacity (kWh) lower bound
        Qu Battery capacity (kWh) upper bound
        socl: Minimum SoC value
        socu: Maximum SoC value
    ev_nfparam:
        soci: initial SoC
        socd: demanded SoC
        ts1: start charging time, [24H]
        ts2: start charging time, [24H]
        tf1: finish charging time, [24H]
        tf2: finish charging time, [24H]
    '''

    def report(self, is_report=True):
        """
        Report EVA info.

        Parameters
        ----------
        is_report: bool
            If True, report EVA info.
        """
        # --- EV summary info ---
        self.Q = self.ev.Q.sum()/1e3
        self.ev['wQ'] = self.ev['Q'] * self.ev['u']
        self.wQ = self.ev['wQ'].sum()/1e3
        self.ev['wsoc'] = self.ev['soc'] * self.ev['Q'] * self.ev['u']
        self.wsoc = self.ev['wsoc'].sum() / self.wQ / 1e3
        cond_ol = self.ev['u'] > 0
        cond_c = self.ev['c'] > 0
        cond_d = self.ev['c'] < 0
        cond_nlc = self.ev['lc'] < 1
        self.data['Pcc'] = -1 * self.ev[cond_ol & cond_c]['Pc'].sum()/1e3
        self.data['Pdc'] = self.ev[cond_ol & cond_d]['Pd'].sum()/1e3
        self.data['Ptc'] = self.data['Pdc'] + self.data['Pcc']

        self.data['Pcc2'] = -1 * self.ev[cond_ol & cond_c & cond_nlc]['Pc'].sum()/1e3
        self.data['Pdc2'] = self.ev[cond_ol & cond_d & cond_nlc]['Pd'].sum()/1e3
        self.data['Ptc2'] = self.data['Pdc2'] + self.data['Pcc2']
        self.ev.drop(columns=['wQ', 'wsoc'], inplace=True)

        cid = self.ev['c'][self.ev['u'] == 1].value_counts().to_dict()
        msg_c = 'Ctrl: '
        for key in cid.keys():
            msg_c += f"{key}={cid[key]}; "

        if is_report:
            # --- report info ---
            msg_time = f"{self.config['name']}: ts={np.round(self.data['ts'], 4)}[H], {self.config['N']} EVs, Total Q={self.Q.round(2)} MWh\n"
            msg_soc = f"Online {self.data['ne']}, Q={self.wQ.round(2)} MWh, SoC={self.wsoc.round(4)}\n"
            msg_p = f"Power(MW): Pt={self.data['Ptc'].round(4)}, Pc={self.data['Pcc'].round(4)}, Pd={self.data['Pdc'].round(4).round(4)}\n"
            logger.warning(msg_time + msg_soc + msg_p + msg_c)

    def g_ts(self, ts):
        """
        Update time and time series.

        Parameters
        ----------
        ts: float
            current time (in 24H).
        """
        # NOTE: tolerance is 0.1s, and large than 24 will be reduced with 24.
        if abs(ts - self.data['ts']) < self.config['t_tol']:
            logger.warning(f"{self.config['name']}: {ts}[H] is too close to current time={self.data['ts']} [H]")
        else:
            ts = ts if ts < 24 else ts-24
        return ts

    def __init__(self, ts=0, N=10000, step=1, tp=100, t_agc=4,
                 lr=0.1, lp=100, seed=None, name="EVA",
                 n_pref=1, is_report=True,
                 tt_mean=0.5, tt_var=0.02, tt_lb=0, tt_ub=1,
                 ict=False, ecc=False, agc=False, t_dur=1):
        """
        Notes
        -----
        1. For efficiency, the EVs that are out of the time range
        [``ts``, ``ts+t_dur``] will be droped.

        1. ``n_pref`` is the number of pereference levels, lower level is
        more likely to response AGC. EVs are evenly distributed on each
        level.
        [NOT SOLVED] This feature will result in control error if ``n_pref`` > 1.

        Parameters
        ----------
        ts: float
            Current time in hour, [0, 24].
        N: int
            Number of EVs
        tp: int
            SSM update period (second).
        t_agc: int
            AGC update period (second).
        lr: float
            learning rate of updating SSM A
        lp: int
            SSM update length, data length for update A.
        step: int
            Step size in seconds.
        seed: int
            Random seed. ``None`` for random.
        name: str
            EVA name.
        n_pref: int
            Number of pereference level, n_pref >= 1.
        is_report: bool
            Paselfd to ``report()``, True to report EVA info.
        ict: bool
            True to turn on increasing charging time control.
        ecc:
            True to turn on error correction.
        """
        # --- 1. init ---
        # TODO: improve documenation
        # --- sl/su: soc charging/dicharging boundary ---
        self.config = OrderedDict(ict=ict, ecc=ecc, agc=agc,
                                  N=N, Ns=20, step=step, tp=tp, n_pref=n_pref,
                                  seed=seed, t_dur=t_dur, t_agc=t_agc,
                                  Np=int(tp/step), lr=lr, lp=lp,
                                  t_tol=0.1/3600, Pdbd=1e-4, er_tol=1, iter_tol=100,
                                  socl=0.005, socu=0.995, name=name)
        self.data = OrderedDict(ts=ts, Pi=0, Pr=0, Prc=0, Per=0, Pet=0)
        # Pdbd: deadband of Power
        self.build(tt_mean=tt_mean, tt_var=tt_var,
                   tt_lb=tt_lb, tt_ub=tt_ub)
        self.report(is_report=is_report)
        d2 = dict(ne=self.ev['u'].sum())
        self.data = {**self.data, **d2}

        # --- output: estimated FRC ---
        self.prumax = 0
        self.prdmax = 0
        self.g_u()
        [self.data['Pt'], _, _, _, _] = self.ep(ne=self.data['nec'])
        self.g_frc()

        # time series data
        self.tsd = pd.DataFrame(data=self.data, index=[0])
        # --- adaptive soc coeff ---
        sdu = self.ev['socd'].max()
        sdl = self.ev['socd'].min()
        Qave = self.ev['Q'].mean() / 1e3  # MWh
        ncave = self.ev['nc'].mean()
        self.Th = (sdu - sdl) * Qave / (2 * self.Pave * ncave)

    def build(self, tt_mean, tt_var, tt_lb, tt_ub):
        """
        Build the ev DataFrame.

        Returns
        -------
        ev: pandas.DataFrame
            EV dataset
        """
        # --- 1. generate EV dataset ---
        self.ev = pd.DataFrame()

        # --- 1a. uniform distribution parameters range ---
        ev_ufparam = dict(Pl=5.0, Pu=7.0,
                          nl=0.88, nu=0.95,
                          Ql=20.0, Qu=30.0,
                          socl=0, socu=1)

        cols = ['Pc', 'Pd', 'nc', 'nd', 'Q']
        cols_bound = {'Pc':   ['Pl', 'Pu'],
                      'Pd':   ['Pl', 'Pu'],
                      'nc':   ['nl', 'nu'],
                      'nd':   ['nl', 'nu'],
                      'Q':    ['Ql', 'Qu']}
        np.random.seed(self.config['seed'])
        for col in cols:
            idxl = cols_bound[col][0]
            idxh = cols_bound[col][1]
            self.ev[col] = np.random.uniform(
                low=ev_ufparam[idxl],
                high=ev_ufparam[idxh],
                size=self.config['N'])
        # Assumption: (1) Pd == Pc; (2) nd == nc;
        self.ev['Pd'] = self.ev['Pc']
        self.ev['nd'] = self.ev['nc']
        # soc interval
        unit = ev_ufparam['socu']/self.config['Ns']
        self.soc_intv = {}
        decimal = 4
        for i in range(self.config['Ns']):
            intv_single = [np.around(i*unit, decimal), np.around((i+1)*unit, decimal)]
            self.soc_intv[i] = intv_single
        #  --- 1b. normal distribution parameters range ---
        ev_pdf_name = ['soci', 'socd', 'ts1', 'ts2', 'tf1', 'tf2', 'tt']
        ev_pdf_data = {'mean': [0.3,    0.8,    -6.5,  17.5,   8.9,  32.9, tt_mean],
                       'var': [0.05,   0.03,   3.4,   3.4,   3.4,  3.4, tt_var],
                       'lb': [0.2,    0.7,    0.0,   5.5,    0.0,  20.9, tt_lb],
                       'ub': [0.4,    0.9,    5.5,   24.0,   20.9, 24.0, tt_ub],
                       'info':  ['initial SoC', 'demanded SoC',
                                 'start charging time 1', 'start charging time 2',
                                 'finish charging time 1', 'finish charging time 2',
                                 'tolerance of increased charging time']}
        ev_pdf = pd.DataFrame(data=ev_pdf_data, index=ev_pdf_name).transpose()
        for col in ev_pdf_name:
            self.ev[col] = stats.truncnorm(
                (ev_pdf[col]['lb'] - ev_pdf[col]['mean']) / ev_pdf[col]['var'],
                (ev_pdf[col]['ub'] - ev_pdf[col]['mean']) / ev_pdf[col]['var'],
                loc=ev_pdf[col]['mean'], scale=ev_pdf[col]['var']).rvs(self.config['N'],
                                                                       random_state=self.config['seed'])

        # ---1c. other data: ts1, ts2, tf1, tf2 ---
        et = self.ev.copy()
        r1 = 0.5  # ratio of t1
        tp1 = self.ev[['ts1', 'tf1']].sample(n=int(et.shape[0]*r1), random_state=self.config['seed'])
        tp2 = self.ev[['ts2', 'tf2']].sample(n=int(et.shape[0]*(1-r1)), random_state=self.config['seed'])
        tp = pd.concat([tp1, tp2], axis=0).reset_index(drop=True).fillna(0)
        tp['ts'] = tp['ts1'] + tp['ts2']
        tp['tf'] = tp['tf1'] + tp['tf2']
        check = tp['ts'] > tp.tf
        row_idx = tp[check].index
        mid = tp.tf.iloc[row_idx].values
        tp.tf.iloc[row_idx] = tp['ts'].iloc[row_idx]
        tp['ts'].iloc[row_idx] = mid
        check = tp['ts'] > tp.tf
        self.ev['ts'] = tp['ts']
        self.ev['tf'] = tp['tf']
        self.ev['u'] = 1

        # Initialize delta power
        self.ev['dP'] = 0

        self.states = list(itertools.product([1, 0, -1], self.soc_intv.keys()))
        self.states_str = [str(s[1])+'S'+str(s[0]) for s in self.states]

        # --- update soc interval and online status ---
        # TODO: do we need to consider the AGC participation?
        # soc is initialized considering random behavior
        # required time to charge to socd
        tr = (self.ev["socd"] - self.ev["soci"]) * self.ev["Q"] / self.ev["Pc"] / self.ev["nc"]
        # stay time
        tc = self.data['ts'] - self.ev["ts"]
        tc[tc < 0] = 0  # reset negative time to 0
        # charge
        self.ev['soc'] = self.ev["soci"]+tc*self.ev["Pc"]*self.ev["nc"]/self.ev["Q"]
        # ratio of stay/required time
        kt = tc/tr
        kt[kt < 1] = 1
        mask = kt[kt > 1].index
        # higher than required charging time, log scale higher than socd
        socp = self.ev["socd"]+np.log(kt)*(1-self.ev["socd"])
        self.ev.loc[mask, 'soc'] = socp.iloc[mask]
        # reset out of range EV soc
        self.ev.loc[self.ev['soc'] >= 1, 'soc'] = 1
        self.ev.loc[self.ev['soc'] <= 0, 'soc'] = 0
        mask_bd1 = self.ev['soc'].values <= self.config['socl']
        mask_bd2 = self.ev['soc'].values >= self.config['socu']
        mask_bd = mask_bd1 | mask_bd2
        self.ev['bd'] = 0
        self.ev.loc[mask_bd, 'bd'] = 1  # indicator if in range

        self.ev['lc'] = 0  # low charge, 0 for regular, 1 for low charge
        # --- ev online status: u0 as u ---
        self.ev['u0'] = 0
        self.g_u()
        self.ev['u0'] = self.ev['u0']

        # Initialize control signal, randomly assign c/dc
        self.ev['c'] = 1
        self.ev['c2'] = 0
        self.g_c(is_test=True, Pi=0)
        # `IS` for demanded charging SoC EVs
        self.ev['c'] = self.ev[['soc', 'c', 'socd']].apply(
            lambda x: 0 if (x[0] >= x[2]) & (x[1] == 1) else x[1], axis=1)
        self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)

        self.g_x()
        self.r_state()

        # initialize x series
        self.ev['xl'] = [[[], [], []]] * self.config['N']
        self.data['ne'] = self.ev['u'].sum()  # number of online EVs

        # pereference
        # Note: pereference is chosen from a list [0, 1, 2, n_pref]
        # the probability of each level is evenly distributed
        # lower one is preferred to be firstly chosen to response AGC
        self.pref = list(range(self.config['n_pref']))
        self.rho = [1 / self.config['n_pref']] * self.config['n_pref']
        self.ev['pref'] = np.random.choice(self.pref, p=self.rho,
                                           size=self.config['N'])

        ev_cols = ['u', 'u0',  'soc', 'bd', 'c', 'c2', 'c0', 'sx', 'dP', 'xl',
                   'soci', 'socd', 'Pc', 'Pd', 'nc', 'nd', 'Q', 'ts', 'tf', 'tt',
                   'pref', 'lc']
        self.ev = self.ev[ev_cols]
        self.ev[['ts', 'tf', 'tt']] = self.ev[['ts', 'tf', 'tt']].astype(float)
        self.ev['agc'] = 0  # `agc` is indicator of participation of AGC
        self.ev['mod'] = 0  # `mod` is indicator of modified of over or under charge

        # load history data of ina
        if self.config['ict']:
            ina = np.genfromtxt('ev_ina_ict.csv', delimiter=',')
        else:
            ina = np.genfromtxt('ev_ina.csv', delimiter=',')

        # initialization of number of actions; truncated normal distribution;
        sx0 = np.ceil(self.ev['soc'] / (1 / self.config['Ns'])) - 1
        # size of each sx
        sx0d = sx0.value_counts().sort_index()
        for i in sx0d.index:
            i = int(i)
            a, b = ina[i, 2], ina[i, 3]
            if a == b:
                b = a + 0.01
            pdf = stats.norm(loc=0, scale=ina[i, 1])
            res = pdf.rvs(sx0d[float(i)], random_state=self.config['seed']).round(0)
            mask = self.ev[sx0 == i].index
            self.ev.loc[mask, 'na'] = ina[i, 0] * (self.ev['tf'].iloc[mask] - self.data['ts']) + res
        mask = self.ev[(self.ev['soc'] < self.ev['socd']) & (self.ev['na'] < 0)].index
        na0 = 1000 * (self.ev['socd'] - self.ev['soc'])
        self.ev.loc[mask, 'na'] = na0.iloc[mask]
        # DEBUG: scale up soc [0.6, 0.7] na0
        mask = self.ev[(self.ev['soc'] < 0.7) & (self.ev['soc'] > 0.5)].index
        self.ev.loc[mask, 'na'] = na0.iloc[mask] * 10
        # for fully charged EVs, reset their na to 0
        mask = self.ev[(self.ev['soc'] >= self.ev['socd'])].index
        self.ev.loc[mask, 'na'] = 0
        self.ev['na'] = self.ev['na'].astype(int)
        self.ev['nama'] = 0  # `nama` is the number of action mileage

        # initialization of max number of actions;
        pcn = self.ev['Pc'] * self.ev['nc']
        self.ev['nam'] = ((self.ev['tf'].mean() - self.ev['ts'].mean() + self.ev['tt']) * pcn
                          - self.ev['socd'] * self.ev['Q']) / (pcn * self.config['t_agc'] / 3600)
        self.ev['nam'] = self.ev['nam'].astype(int)

        if self.config['ict']:
            mask = self.ev[(self.ev['na'] >= self.ev['nam'])].index
            self.ev.loc[mask, 'na'] = self.ev['nam'].iloc[mask]
            self.ev.loc[mask, 'lc'] = 1
        cond2 = self.ev['soc'] <= self.config['socl']
        self.ev.loc[cond2, 'lc'] = 1
        self.g_u()
        # self.data['nec'] = self.ev['u'].sum() - self.ev["lc"].sum()  # numebr of online EVs with lc == 0

        self.g_BCD()
        self.n_step = 1
        # --- drop EVs that are not in time range ---
        set0 = set(self.ev.index)
        set1 = set(self.ev[(self.ev['tf'] <= self.data['ts'])].index)
        set2 = set(self.ev[(self.ev['ts'] >= self.data['ts'] + self.config['t_dur'])].index)
        set_in = set0 - set1 - set2
        self.ev = self.ev.iloc[list(set_in)].reset_index(drop=True)
        self.ev['soc0'] = self.ev['soc']
        self.ev['na0'] = self.ev['na'].astype(int)
        self.g_x()
        self.r_state()
        return True

    def g_u(self):
        """
        Update online status of EV at given time.
        """
        # --- online status ---
        self.ev['u0'] = self.ev['u'].values.astype(int)
        self.ev['u'] = (self.ev['ts'] <= self.data['ts']) & (self.ev['tf'] >= self.data['ts'])
        self.ev['u'] = self.ev['u'].astype(int)
        self.data['ne'] = self.ev['u'].sum()  # number of online EVs
        # number of online and non low-charge EVs
        self.data['nec'] = self.ev[(self.ev['lc'] == 0) & (self.ev['u'] == 1)].shape[0]
        return True

    def g_x(self):
        """
        Update EV x and SSM x.

        In the output table, columns stand for soc interval, rows stand for charging status.

        0 for charging, 1 for idle, 2 for discharging.

        Returns
        -------
        xtab: pd.DataFrame
            Table of SSM status of *online* EVs
        rtab: pd.DataFrame
            Table of SSM status of *online* EVs, by percentage
        """
        # --- find single EV sx ---
        self.ev['sx'] = np.ceil(self.ev['soc'] / (1 / self.config['Ns'])) - 1
        self.ev['sx'] = self.ev['sx'].astype(int)

        # --- update c2 ---
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})

        # --- output tab ---
        # TODO: consider low charged path?
        cond_u = self.ev['u'] == 1
        cond_nlc = self.ev['lc'] == 0
        states = pd.Series([(-1, -1)] * self.ev.shape[0])
        mask = states[cond_u & cond_nlc].index
        states.loc[mask] = pd.Series(list(zip(self.ev['c2'], self.ev['sx']))).loc[mask]
        res = dict(states.value_counts())
        self.xtab = pd.DataFrame(columns=range(self.config['Ns']), index=[0, 1, 2], data=0)
        for key in res.keys():
            if key[1] > -1:
                self.xtab.loc[key[0], key[1]] = res[key]
        self.xtab.fillna(0, inplace=True)
        all_ev = self.xtab.sum().sum()
        # TODO: use estimated ne rather than actual ne?
        self.rtab = self.xtab.div(all_ev)
        return True

    def save_A(self, csv):
        """
        Save A matrix to csv file.

        Parametes
        ---------
        csv: str
            csv file name.
        """
        As = pd.DataFrame(data=self.A)
        As.to_csv(csv, index=False)
        logger.warning(f"{self.config['name']}: Save A to %s." % csv)
        return True

    def load_A(self, csv):
        """
        Load A matrix from csv files.

        Warning: load A from csv files will overwrite existing A!

        Parametes
        ---------
        csv: str
            csv file name.
        """
        self.A = pd.read_csv(csv).values
        logger.warning(f"{self.config['name']}: Load A from %s." % csv)

    def plot_agc(self, figsize=(6, 3), style='default', tu='s', **kwargs):
        """
        Plot the AGC results.

        Parameters
        ----------
        figsize: tuple
            Figure size.
        style: str
            Plt style.
        tu: str
            Time unit, 'h' for hour, 's' for second.
        """
        plt.style.use(style)
        fig_agc, ax_agc = plt.subplots(figsize=figsize, **kwargs)
        self.tsd['ts2'] = 3600 * (self.tsd['ts'] - self.tsd['ts'].iloc[0])
        if tu == 's':
            x = 'ts2'
            xlabel = 'Time [s]'
        elif tu == 'h':
            x = 'ts'
            xlabel = 'Time [H]'
        self.tsd.plot(x=x, y=['Pi', 'Prc'],
                      label=['Control', 'Response'],
                      ax=ax_agc)

        ax_agc.set_xlabel(xlabel)
        ax_agc.set_ylabel('Power (MW)')
        ax_agc.set_title('AGC response')
        ax_agc.set_xlim(self.tsd[x].iloc[0], self.tsd[x].iloc[-1])
        ax_agc.grid()
        ax_agc.legend()
        self.tsd.drop(columns=['ts2'], inplace=True)
        return fig_agc, ax_agc

    def plot(self, figsize=(6, 3), style='default', tu='s',
             bbox_to_anchor=(0.85, 0.2), loc='lower right', **kwargs):
        """
        Plot the results.

        Parameters
        ----------
        figsize: tuple
            Figure size.
        style: str
            Plt style.
        tu: str
            Time unit, 'h' for hour, 's' for second.
        """
        plt.style.use(style)
        fig_ev, ax_ev = plt.subplots(figsize=figsize, **kwargs)
        self.tsd['ts2'] = 3600 * (self.tsd['ts'] - self.tsd['ts'].iloc[0])
        if tu == 's':
            x = 'ts2'
            xlabel = 'Time [s]'
        elif tu == 'h':
            x = 'ts'
            xlabel = 'Time [H]'
        self.tsd.plot(x=x, y=['Ptc', 'Pcc', 'Pdc'],
                      label=['Total', 'Charging', 'Discharging'],
                      ax=ax_ev, legend=False)
        ax2 = ax_ev.twinx()
        self.tsd.plot(x=x, y='ne', label='Online EVs',
                      color='orange', ax=ax2,
                      legend=False)
        ax2.set_ylabel('Online EVs')
        ax_ev.set_xlabel(xlabel)
        ax_ev.set_ylabel('Power (MW)')
        ax_ev.set_title(f"{self.config['name']}")
        ax_ev.set_xlim(self.tsd[x].iloc[0], self.tsd[x].iloc[-1])
        ax_ev.grid()
        fig_ev.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, **kwargs)
        self.tsd.drop(columns=['ts2'], inplace=True)
        return fig_ev, ax_ev

    def test(self, tf=9.05):
        """
        Build the matrix A with test ev data.

        During the test, the EV behaviors and controls will be uniformed.
        After the test, the EV will be reset to initial states.

        Parameters
        ----------
        ts: float
            start time [H]
        tf: float
            end time [H]
        """
        t_step = self.config['step'] / 3600
        ev_cols = ['u', 'u0', 'soc', 'sx', 'bd',
                   'ts', 'tf', 'c', 'c0', 'dP', 'c2']
        ev_copy = self.ev[ev_cols].copy()

        self.ev['tf'] = 23.9
        self.ev['ts'] = 0.1
        self.ev['soc'] = np.random.uniform(low=0.001, high=0.999, size=self.config['N'])
        self.g_u()

        # initialize control signal: randomly assign C/I/D
        self.ev['c'] = np.random.choice([-1, 0, 1], p=[0.33, 0.33, 0.34], size=self.config['N'])
        # revise
        # offline -> I
        self.ev['c'] = self.ev[['u', 'c']].apply(lambda x: x[1] if x[0] != 0 else 0, axis=1)
        # sl, DC -> C
        self.ev['c'] = self.ev[['u', 'c', 'soc']].apply(
            lambda x: 1 if (x[2] < self.config['socl']) & (x[1] == -1) else x[1], axis=1)
        # demanded SoC, C -> I
        # self.ev['c'] = self.ev[['u', 'c', 'soc', 'socd']].apply(
        #     lambda x: 1 if (x[2] >= x[3]) & (x[1] == 1) else x[1], axis=1)

        # record control signal
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        self.ev['c0'] = self.ev['c2']
        self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)

        # --- update x ---
        self.g_x()
        self.g_xl()

        for t in tqdm(np.arange(self.data['ts']+t_step, tf, t_step), desc=f"{self.config['name']} MCS"):
            self.data['ts'] = self.g_ts(t)
            self.g_u()  # update online status
            self.g_c(is_test=True)  # update control signal
            # --- update soc interval and online status ---
            # charging/discharging power, kW
            self.ev['dP'] = self.ev[['Pc', 'Pd', 'nc', 'nd', 'c', 'u']].apply(
                lambda x: x[0]*x[2]*x[5]*x[4] if x[4] >= 0 else x[1]*x[3]*x[5]*x[4], axis=1)
            # --- update and modify SoC ---
            self.ev['soc'] = self.ev.soc + t_step * self.ev['dP'] / self.ev['Q']
            masku = self.ev['soc'].values >= self.config['socu']
            maskl = self.ev['soc'].values <= self.config['socl']
            self.ev.loc[masku, 'soc'] = self.config['socu']
            self.ev.loc[maskl, 'soc'] = self.config['socl']
            self.ev['soc'] = self.ev.soc + t_step * self.ev['dP'] / self.ev['Q']
            # --- boundary ---
            self.ev['bd'] = self.ev[['soc', 'socd']].apply(
                lambda x: 1 if (x[0] <= self.config['socl']) | (x[0] >= x[1]) else 0)

            # record control signal
            self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
            self.ev['c0'] = self.ev['c2']
            # format
            self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)
            # update control signal
            # offline -> I
            self.ev['c'] = self.ev[['u', 'c']].apply(lambda x: x[1] if x[0] != 0 else 0, axis=1)
            # sl, DC -> C
            self.ev['c'] = self.ev[['u', 'c', 'soc']].apply(
                lambda x: 1 if (x[2] < self.config['socl']) & (x[1] == -1) else x[1], axis=1)
            # demanded SoC, C -> I
            self.ev['c'] = self.ev[['u', 'c', 'soc', 'socd']].apply(
                lambda x: 1 if (x[2] >= x[3]) & (x[1] == 1) else x[1], axis=1)

            # --- update x ---
            self.g_x()
            self.g_xl()

        # build A matrix
        self.g_A()

        # reset EV data
        self.ev[ev_cols] = ev_copy[ev_cols]
        # self.reset(ts0, clean_xl=False)

    def run(self, tf=10, Pi=0,
            is_updateA=False, is_rstate=False,
            is_test=False, disable=False):
        """
        Run the ev aggregator from ``ts`` to ``tf``.

        Parameters
        ----------
        tf: int
            running time [H]
        Pi: float
            AGC input signal (MW)
        is_update: bool
            `True` for updating SSM A during the simulation, False for not.
            Set `is_g_SSM=True` will record EV status.
        is_record: bool
            `True` for recording the EV status in a series, False for not recording.
            Set `is_record=False` can speed up the simulation.
        is_test: bool
            `g_c` control mode
        disable: bool
            tqdm progress bar
        """
        t_step = self.config['step'] / 3600  # t_step is in hours
        if tf - self.data['ts'] < 1e-5:
            logger.warning(f"{self.config['name']}: end time {tf}[H] is too close to start time {self.data['ts']}[H]," +
                           "simulation will not start.")
        else:
            for t in tqdm(np.arange(self.data['ts'], tf + 0.1/3600, t_step),
                          desc=f"{self.config['name']} MCS",
                          disable=disable):
                if abs(t - self.data['ts']) < self.config['t_tol']:
                    continue
                # --- update SSM A ---
                Pi_input = Pi
                if self.n_step % self.config['Np'] == 0:
                    if is_updateA:
                        self.g_A(is_update=True)
                    if is_rstate:
                        self.r_state()
                        self.data['Per'] = self.data['Pi'] - self.data['Prc']  # error of AGC response
                        self.data['Pet'] = self.data['Pt'] - self.data['Ptc']  # error of total output power
                    self.report(is_report=False)
                self.data['ts'] = self.g_ts(t)
                self.g_u()  # update online status
                # TODO: add warning when Pi is 0
                Pi_input += self.config['ecc'] * (self.data['Pet'])  # DEBUG: remove term ` + self.data['Pet']`
                self.g_c(Pi=Pi_input, is_test=is_test)
                lc0 = self.ev['lc'].copy()
                self.g_c(Pi=Pi_input, is_test=is_test)  # update control signal
                # --- update soc interval and online status ---
                # charging/discharging power, kW
                self.ev['dP'] = 0
                masku = self.ev[self.ev['c'].astype(float) > 0].index.astype(int)
                maskl = self.ev[self.ev['c'].astype(float) < 0].index.astype(int)
                self.ev.loc[masku, 'dP'] = self.ev.loc[masku, 'Pc'] * self.ev.loc[masku, 'nc']
                self.ev.loc[maskl, 'dP'] = self.ev.loc[maskl, 'Pd'] * self.ev.loc[maskl, 'nd']
                self.ev['dP'] = self.ev['dP'] * self.ev['u']
                # --- update and modify outranged SoC ---
                self.ev['soc'] = self.ev['soc'] + t_step * self.ev['dP'] / self.ev['Q']
                masku = self.ev[self.ev['soc'].astype(float) >= 1].index
                maskl = self.ev[self.ev['soc'].astype(float) <= 0].index
                self.ev.loc[masku, 'soc'] = 1
                self.ev.loc[maskl, 'soc'] = 0

                # --- update x ---
                self.g_x()

                if is_updateA:
                    self.g_xl()
                self.report(is_report=False)  # record power
                # Actual AGC response: AGC switched power if not modified. (MW)
                self.data['Prc'] = np.sum(self.ev['agc'] * self.ev['Pc'] *
                                          (1 - self.ev['mod']) * (1 - lc0)) * 1e-3
                self.data['Pi'] = Pi
                self.data['Pt'] = self.Pt

                self.tsd = pd.concat([self.tsd, pd.DataFrame(data=self.data, index=[0])],
                                     ignore_index=True)
                self.tsd.iloc[-1, 1] = Pi  # col 'Pr'
                self.n_step += 1
            # self.ict()

    def g_A(self, is_update=False):
        """
        Compute A matrix: cols: x(k), row: x(k+1)
        The sum of col should be 1.

        Parameters
        ----------
        is_update: bool
            incremental update. `True` for updating SSM A, `False` for generating SSM A.
        """
        # --- gather results ---
        states = []
        ctrls = []
        for item in self.ev.xl.tolist():
            if len(item[0]) > 0:
                states.append(item[1][0])
                ctrls.append(item[0][0])
        data = []
        for x, y in zip(ctrls, states):
            data0 = []
            for c, s in zip(x, y):
                rx = c * self.config['Ns'] + s
                data0.append(rx)
            data.append(data0)

        A0 = np.zeros((3*self.config['Ns'], 3*self.config['Ns']))
        for d in data:
            for i in range(len(d)-1):
                A0[d[i+1], d[i]] += 1

        # TODO: Consider data length limit: self.config['lp']
        if is_update:
            n = int(self.config['lr'] / (1-self.config['lr']))
            A0 = self.A * len(data) * n * np.ones((60,))
            self.A = self.A + A0
            row_sum = self.A.sum(axis=0)
            self.A = self.A / row_sum
        else:
            row_sum = A0.sum(axis=0)
            self.A = A0 / row_sum

    def g_xl(self):
        """
        Update EV x series.
        """
        self.ev['tnow'] = self.data['ts']
        col = ['sx', 'xl', 'u', 'u0', 'c2', 'tnow', 'c0', 'bd']
        self.ev['xl'] = self.ev[col].apply(update_xl, axis=1)
        self.ev.drop(['tnow'], axis=1, inplace=True)
        return True

    def g_c(self, Pi=0, is_test=False):
        """
        Generate the charging signal.
        `is_test=True` is recommended for initially building SSM A.

        EV start to charge with rated power as soon as it plugs in.
        The process will not be interrupted until receive control signal
        or achieved demanmded SoC level.

        Parameters
        ----------
        is_test: bool
            `True` for test mode, `False` for normal mode.
            normal mode: ``CS`` for online EVs
            test mode: only revise control signal
        """
        # record last reformatted ctrl signal
        self.ev['c0'] = self.ev['c2'].copy()
        # `CS` for just arrived EVs
        mask = self.ev[(self.ev['u0'] == 0) & (self.ev['u'] == 1)].index
        self.ev.loc[mask, 'c'] = 1

        if is_test:
            pass
        else:
            if self.config['agc'] & (abs(Pi) > self.config['Pdbd']):  # deadband
                self.r_agc(Pi=Pi)
            else:
                # --- revise control ---
                # `CS` for low charged EVs, and set 'lc' to 1
                mask = self.ev[(self.ev['soc'] <= self.config['socl']) & (self.ev['u']) == 1].index
                self.ev.loc[mask, ['lc', 'c']] = 1
                self.ev['mod'] = 0
                self.ev.loc[mask, 'mod'] = 1
                # `IS` for full EVs
                mask = self.ev[(self.ev['soc'] >= self.ev['socd'])].index
                self.ev.loc[mask, 'c'] = 0
        # `IS` for offline EVs
        self.ev['c'] = self.ev['c'].astype(float) * self.ev['u'].astype(float)
        # reformatted control signal c2
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        # format
        self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)

    def g_res(self, x0, n=1):
        """
        Estimate EV_SSM status `x0`, stores in attribute pd.DataFrame ``etab``

        Parameters
        ----------
        x0: numpy.array, (60, )
            [CS, DS, IS], initial distribution of states.
        n: int
            number of steps.

        Returns
        -------
        xr: numpy.array, (60, )
            [CS, DS, IS], estimated distribution of states.
        etab: pd.DataFrame
            the estimated states of the next n steps.
        """
        if not hasattr(self, 'A'):
            raise ValueError('Matrix A is not defined')
        # --- out ---
        An = np.linalg.matrix_power(self.A, n=n)
        self.x0 = np.matmul(An, x0)
        self.etab = pd.DataFrame(data=self.x0.reshape(3, 20),
                                 columns=range(20),
                                 index=range(3))

    def reset(self):
        """
        NOT TEST YET.
        Reset the pd.DataFrame ``ev`` to the initial time.

        Parameters
        ----------
        tnow: int
            the current time [H].
        clean_xl: bool
            `True` to clean `ev['xl']`, `False` to keep it.
        """
        A = self.A.copy()
        # TODO: reset with __init__ function
        self.__init__(self, ts=0, N=10000, step=1, tp=100, t_agc=4,
                      lr=0.1, lp=100, seed=None, name='EVA',
                      n_pref=1, is_report=True,
                      tt_mean=0.5, tt_var=0.02, tt_lb=0, tt_ub=1,
                      ict=True, ecc=True, t_dur=1)

    def r_state(self):
        """
        Record the states distribution of online EVs `x0` from `rtab`.

        Returns
        -------
        x0: np.array
            states distribution (60, ) [CS, IS, DS]
        """
        if hasattr(self, 'rtab'):
            self.x0 = self.rtab.values.reshape(-1,)
        else:
            self.x0 = np.zeros(60,)
            logger.warning(f"{self.config['name']}: `rtab` is not available!")
        return self.x0

    def cp(self):
        """
        Analytical output power (MW).
        Total, charge, discharge, upper, lower

        Returns
        -------
        out: list of float
            [Pt, Pc, Pd, Pu, Pl] (MW)
        """
        masku = self.ev[(self.ev['u'] == 1) & (self.ev['c'].values > 0)].index
        maskl = self.ev[(self.ev['u'] == 1) & (self.ev['c'].values < 0)].index
        self.ev['Ps'] = 0
        self.ev.loc[masku, 'Ps'] = self.ev.loc[masku, 'Pc']
        self.ev.loc[maskl, 'Ps'] = -1 * self.ev.loc[maskl, 'Pd']

        Pt = - self.ev['Ps'].to_numpy().sum()
        Pc = - self.ev['Ps'][self.ev['Ps'].values > 0].sum()
        Pd = self.ev['Ps'][self.ev['Ps'].values < 0].sum()
        ue1_ridx = self.ev[self.ev['u'].values == 1].index
        Pl = - self.ev['Pc'].iloc[ue1_ridx].sum()
        Pu = self.ev['Pd'].iloc[ue1_ridx].sum()
        self.ev.drop(['Ps'], axis=1, inplace=True)

        out = [Pt, Pu, Pl, Pc, Pd]
        out = [x/1000 for x in out]  # kW to MW
        return out

    def ep(self, ne=None):
        """
        Estimate output power (MW).

        Returns
        -------
        out: list of float
            [Pt, Pu, Pl], total, upper, lower (MW)
        """
        # TODO: `ne` should be replaced with an recorded value
        # ne = self.data['ne']
        if not ne:
            ne = self.data['ne']
        self.Pt = np.matmul(ne * self.D, self.x0)[0]
        self.Pu = np.matmul(ne * self.Db, self.x0)[0]
        self.Pl = np.matmul(ne * self.Dd, self.x0)[0]
        self.Pa = ne * np.matmul(self.Da - self.D, self.x0)[0]
        self.Pb = ne * np.matmul(self.Db - self.Da, self.x0)[0]
        self.Pc = ne * np.matmul(self.Dc - self.D, self.x0)[0]
        self.Pd = ne * np.matmul(self.Dd - self.Dc, self.x0)[0]
        out = [self.Pt, self.Pu, self.Pl, self.Pa, self.Pc]
        return out

    def g_frc(self, nea=None):
        """
        Estimate frequency regulation capacity (FRC)[MW].

        Parameters
        ----------
        nea: int
            adjusted EV numbers, multiplier of the FRC,
            'nea' is set to `self.data['ne']` if not given.
            FRC = FRC0 * nea / ne

        Returns
        -------
        out: list of float
            [prumax, prdmax] (MW)
        """
        self.report(is_report=False)
        [_, _, _, _, _] = self.ep(ne=self.data['nec'])
        RU = self.Pu - self.Pt
        RD = self.Pt - self.Pl
        # --- demanded-SOC limited FRC is not that reasonable ---
        # TODO: 0.8 is the estimated soc demanded, may need revision
        # self.prumax = max(min((self.wsoc - 0.8) * self.wQ / T, RU), 0)
        # self.prdmax = min((1 - self.wsoc) * self.wQ / T, RD)
        # --- SSM FRC is not that reasonable ---
        if not nea:
            nea = self.data['nec']
        self.prumax = RU * nea / self.data['ne'] * self.data['nec'] / self.data['ne']
        self.prdmax = RD * nea / self.data['ne'] * self.data['nec'] / self.data['ne']
        return [self.prumax, self.prdmax]

    def g_BCD(self):
        """
        Build SSM B, C, D matrix.
        Different from the paper, the D matrix does not contain `ne`.

        The output of EVA (MW) `Pt=ne*D*x0`

        The output of mode {m} ``P{m}=ne*D{m}*x0``, where `m` in [a, b, c, d]
        """
        # --- B ---
        B1 = -1 * np.eye(self.config['Ns'])
        B2 = np.eye(self.config['Ns'])
        B3 = np.zeros((self.config['Ns'], self.config['Ns']))
        self.B = np.vstack((B1, B2, B3))

        # --- C ---
        C1 = np.zeros((self.config['Ns'], self.config['Ns']))
        C2 = -1 * np.eye(self.config['Ns'])
        C3 = np.eye(self.config['Ns'])
        self.C = np.vstack((C1, C2, C3))

        # --- D ---
        kde = stats.gaussian_kde(self.ev.Pc)
        Pave = 0  # P average
        step = 0.01
        # Note: Pd is assumed to be equal to Pc
        for Pl in np.arange(self.ev.Pc.min(), self.ev.Pc.max(), step):
            Pave += (kde.integrate_box(Pl, Pl+step)) * (Pl + 0.005 * step)
        Pave = Pave / 1000  # kW to MW
        self.Pave = Pave
        D1 = -1 * np.ones((1, self.config['Ns']))
        D2 = np.zeros((1, self.config['Ns']))
        D3 = np.ones((1, self.config['Ns']))
        self.D = Pave * np.hstack((D1, D2, D3))

        # --- D a,b,c,d ---
        D1 = np.zeros((1, self.config['Ns']))
        D2 = np.zeros((1, self.config['Ns']))
        D3 = np.ones((1, self.config['Ns']))
        self.Da = Pave * np.hstack((D1, D2, D3))

        D1 = np.ones((1, self.config['Ns']))
        D2 = np.ones((1, self.config['Ns']))
        D3 = np.ones((1, self.config['Ns']))
        self.Db = Pave * np.hstack((D1, D2, D3))
        self.Db[0, self.config['Ns']] = 0  # low charged EV don't DC

        D1 = -1 * np.ones((1, self.config['Ns']))
        D2 = np.zeros((1, self.config['Ns']))
        D3 = np.zeros((1, self.config['Ns']))
        self.Dc = Pave * np.hstack((D1, D2, D3))

        D1 = -1 * np.ones((1, self.config['Ns']))
        D2 = -1 * np.ones((1, self.config['Ns']))
        D3 = -1 * np.ones((1, self.config['Ns']))
        self.Dd = Pave * np.hstack((D1, D2, D3))
        self.Dd[0, 2*self.config['Ns']-1] = 0  # high charged EV don't C

        return True

    def r_agc(self, Pi):
        """
        Alter `ev.ctrl` based on `us` `vs` from `g_agc` to response AGC; update `x0`

        Parameters
        ----------
        Pi: float
            Power input (MW)

        Returns
        -------
        u: numpy.ndarray
            SSM vector `u`, (Ns,); mode (a), (d): CS - IS
        v: numpy.ndarray
            SSM vector `v`, (Ns,); mode (b), (c): IS - DS
        us: numpy.ndarray
            SSM vector `us`, (Ns+1,); probability mode (a), (d): CS - IS
        vs: numpy.ndarray
            SSM vector `vs`, (Ns+1,); probability mode (b), (c): IS - DS
        usp: numpy.ndarray
            (Ns+1, n_pref), probability mode (a), (d) with pereference: CS - IS
        vsp: numpy.ndarray
            (Ns+1, n_pref), probability mode (b), (c) with pereference: IS - DS
        """
        if Pi >= 0:
            Pi_cap = min(Pi, self.prumax)
        elif Pi < 0:
            Pi_cap = max(Pi, -1*self.prdmax)
        # initialize output
        u = np.zeros(self.config['Ns'])
        v = np.zeros(self.config['Ns'])
        us = np.zeros(self.config['Ns']+1)
        vs = np.zeros(self.config['Ns']+1)
        usp = np.repeat(us.reshape(self.config['Ns']+1, 1), self.config['n_pref'], axis=1)
        vsp = np.repeat(vs.reshape(self.config['Ns']+1, 1), self.config['n_pref'], axis=1)
        # corrected control
        iter = 0
        Pr = self.data['Pr']
        error = Pi_cap - Pr
        while (abs(error) >= self.config['er_tol']) & (iter < self.config['iter_tol']):
            error = Pi_cap - Pr
            error0 = error
            Pt0 = self.Pt
            u, v, us, vs = self.g_agc(error)
            # pereference signal
            usp = np.repeat(us.reshape(self.config['Ns']+1, 1), self.config['n_pref'], axis=1)
            vsp = np.repeat(vs.reshape(self.config['Ns']+1, 1), self.config['n_pref'], axis=1)
            for i_sx in range(self.config['Ns']):
                pu = us[i_sx]
                pv = vs[i_sx]
                if (us[-1] == 1) & (vs[-1] == 1):
                    for i in range(len(self.rho)):
                        su = np.sum([m*n for m, n in zip(self.rho[0:i], usp[i_sx, 0:i])])
                        sv = np.sum([m*n for m, n in zip(self.rho[0:i], vsp[i_sx, 0:i])])
                        usp[i_sx, i] = max(min((pu-su)/self.rho[i], 1), 0)
                        vsp[i_sx, i] = max(min((pv-sv)/self.rho[i], 1), 0)
                if (us[-1] == -1) & (vs[-1] == -1):
                    for i in range(len(self.rho)-1, -1, -1):
                        su = np.sum([m*n for m, n in zip(self.rho[i+1:], usp[i_sx, i+1:])])
                        sv = np.sum([m*n for m, n in zip(self.rho[i+1:], vsp[i_sx, i+1:])])
                        usp[i_sx, i] = max(min((pu-su)/self.rho[i], 1), 0)
                        vsp[i_sx, i] = max(min((pv-sv)/self.rho[i], 1), 0)

            # EV switch control states according to aggregator signal
            # --- old version ---
            # self.ev.c = self.ev[['u', 'c', 'sx', 'pref']].apply(lambda x: r_agc_sev(x, us, vs, usp, vsp), axis=1)
            # --- new version ---
            # self.ev['c'] = self.ev['c'] * self.ev['u']  # offline
            self.ev['pv'] = vsp[self.ev['sx'], self.ev['pref']]  # prob for vs
            self.ev['pu'] = usp[self.ev['sx'], self.ev['pref']]  # prob for us

            cond_ol = self.ev['u'] == 1  # online EV
            cond_nlc = self.ev['lc'] == 0  # non-lc EV

            self.ev['p'] = 1.1  # default p
            mask_p = self.ev[cond_ol & cond_nlc].index
            self.ev.loc[mask_p, 'p'] = np.random.uniform(low=0, high=1, size=len(mask_p))
            cond_pv = self.ev['p'] <= self.ev['pv']  # prob of EV to vs
            cond_pu = self.ev['p'] <= self.ev['pu']  # prob of EV to us
            cond_d = self.ev['c'] == -1  # EV in DS
            cond_i = self.ev['c'] == 0  # EV in IS
            cond_c = self.ev['c'] == 1  # EV in CS
            if us[-1] == 1:  # positive signal, I to D, C to I
                maskvs = self.ev[cond_ol & cond_nlc & cond_pv & cond_i].index
                maskus = self.ev[cond_ol & cond_nlc & cond_pu & cond_c].index
                self.ev.loc[maskvs, 'c'] = -1
                self.ev.loc[maskus, 'c'] = 0
            elif us[-1] == -1:  # negative signal, I to C, D to I
                maskvs = self.ev[cond_ol & cond_nlc & cond_pv & cond_d].index
                maskus = self.ev[cond_ol & cond_nlc & cond_pu & cond_i].index
                self.ev.loc[maskvs, 'c'] = 0
                self.ev.loc[maskus, 'c'] = 1
                mask = self.ev[((self.ev['c'] != 0)) & (self.ev['lc'] == 1)].index

            self.g_x()
            # --- record output ---
            # TODO: modification of random traveling behavior
            self.x0 = self.x0 + np.matmul(self.B, u) + np.matmul(self.C, v)
            self.g_u()
            self.ep(ne=self.data['nec'])
            # dx0 = np.matmul(self.B, u) + np.matmul(self.C, v)
            # self.data['Pr'] = np.matmul(self.D, dx0)[0]
            # self.data['Pr'] += np.matmul(self.data['nec'] * self.D,
            #                             np.matmul(self.B, u) + np.matmul(self.C, v))[0]

            Pr += self.Pt - Pt0
            error = Pi_cap - Pr
            iter += 1
            if abs(error0 - error) < self.config['er_tol']:  # tolerante of control error
                break

        self.data['Pr'] = Pr
        # --- deliver agc signal ---
        # `IS` for offline EVs
        self.ev['c'] = self.ev['c'].astype(float) * self.ev['u'].astype(float)
        # update agc
        self.ev['agc'] = 0
        cond_ol = self.ev['u'] == 1
        cond_nlc = self.ev['lc'] == 0
        cond_socl = self.ev['soc'] < self.ev['socd']
        cond_socu = self.ev['soc'] >= self.ev['socd']
        # online and non lc, not full
        mask_nf = self.ev[cond_ol & cond_socl & cond_nlc].index
        self.ev.loc[mask_nf, 'agc'] = 1 - self.ev.loc[mask_nf, 'c']
        # online and non lc, full
        mask_f = self.ev[cond_ol & cond_socu & cond_nlc].index
        self.ev.loc[mask_f, 'agc'] = 0 - self.ev.loc[mask_f, 'c']
        # update counter
        self.ev['na'] += self.ev['agc']
        self.ev['nama'] += self.ev['agc'].abs()
        # --- revise control ---
        # `CS` for low charged EVs, and set 'lc' to 1
        mask = self.ev[(self.ev['soc'] <= self.config['socl']) & (self.ev['u']) == 1].index
        self.ev.loc[mask, ['lc', 'c']] = 1
        self.ev['mod'] = 0
        self.ev.loc[mask, 'mod'] = 1
        # `CS` for lc not full EVs
        mask = self.ev[(self.ev['soc'] < self.ev['socd']) & (self.ev['lc'] == 1)].index
        self.ev.loc[mask, 'c'] = 1
        # `IS` for lc and full EVs
        mask = self.ev[(self.ev['soc'] >= self.ev['socd']) & (self.ev['lc'] == 1)].index
        self.ev.loc[mask, 'c'] = 0
        # update lc according to na when ICT is on
        if self.config['ict']:
            lc0 = self.ev['lc'].values.copy().astype(int)
            self.ev['lc'] = (self.ev['na'] >= self.ev['nam'])
            self.ev['lc'] = self.ev['lc'].values | lc0  # once lc, never AGC
            self.ev['lc'] = self.ev['lc'].astype(int)

        # --- record output ---
        self.uv = [u, v, us, vs, usp, vsp]
        # TODO: ps array
        self.usp = np.repeat(us[0:self.config['Ns']].reshape(self.config['Ns'], 1),
                             self.config['n_pref'], axis=1)
        self.vsp = np.repeat(vs[0:self.config['Ns']].reshape(self.config['Ns'], 1),
                             self.config['n_pref'], axis=1)
        return u, v, us, vs, usp, vsp, Pi_cap, error

    def g_agc(self, Pi):
        """
        Generate control signal `u` `v` `us` `vs` as response to AGC input.

        v1.0, C to I, high to low, I to D;

        v2.0, C to I, I to D, high to low;

        Parameters
        ----------
        Pi: float
            Power input (MW)

        Returns
        -------
        u: numpy.NDArray
            vector `u`, (Ns,); mode (a), (d): CS - IS
        v: numpy.NDArray
            vector `v`, (Ns,); mode (b), (c): IS - DS
        us: numpy.NDArray
            vector `us`, (Ns+1,); probability mode (a), (d): CS - IS
        vs: numpy.NDArray
            vector `vs`, (Ns+1,); probability mode (b), (c): IS - DS
        """
        x = self.x0.copy()

        u = np.zeros(self.config['Ns'])
        v = np.zeros(self.config['Ns'])
        us = np.zeros(self.config['Ns']+1)
        vs = np.zeros(self.config['Ns']+1)

        if Pi >= 0:  # RegUp
            # --- step I ---
            ru = min(Pi, self.Pa) / (self.Pave * self.data['nec'])  # total RegUp power
            u = np.zeros(self.config['Ns'])   # C->I
            v = np.zeros(self.config['Ns'])   # I->D
            for j in range(self.config['Ns']):
                a = ru - np.sum(x[j+1: self.config['Ns']]) - np.sum(x[j+1+self.config['Ns']: 2*self.config['Ns']])
                u[j] = min(max(a, 0), x[j])
                v[j] = min(max(a-x[j], 0), x[j+self.config['Ns']])
            # --- step II ---
            for j in range(self.config['Ns']):
                us[j] = min(safe_div(u[j], x[j]), 1)
                vs[j] = min(safe_div(v[j], x[j+self.config['Ns']]+u[j]), 1)   # TODO: switch two times?
            # --- step III ---
            us[-1] = 1
            vs[-1] = 1

        elif Pi < 0:  # RegDn
            # --- step I ---
            rv = max(Pi, self.Pc) / (self.Pave * self.data['nec'])
            v = np.zeros(self.config['Ns'])
            for j in range(self.config['Ns']):
                a = rv + np.sum(x[2*self.config['Ns']: 2*self.config['Ns']+j])
                v[j] = max(min(a, 0), -1 * x[j+2*self.config['Ns']])

            ru = min(Pi-self.Pc, 0) / (self.Pave * self.data['nec'])
            u = np.zeros(self.config['Ns'])
            for j in range(self.config['Ns']):
                a = ru - np.sum(v[0:j]) - np.sum(u[0:j-1])
                u[j] = max(min(a, 0), -1 * x[j+self.config['Ns']])
            # --- step II ---
            for j in range(self.config['Ns']):
                vs[j] = min(safe_div(-1*v[j], x[j+2*self.config['Ns']]), 1)
                us[j] = min(safe_div(-1*u[j], x[j+self.config['Ns']]-v[j]), 1)
            # --- step III ---
            us[-1] = -1
            vs[-1] = -1

        return u, v, us, vs

    def sac(self, x):
        """
        SoC adaptive coefficient.

        Parameters
        ----------
        x: float
            SoC, [0, 1]

        Returns
        -------
        np.array
            [fcs(x), fis(x), fds(x)]
        """
        sdu = self.ev['socd'].max()
        sdl = self.ev['socd'].min()
        conds = [x < sdl,
                 (x >= sdl) & (x < sdu),
                 x >= sdu]

        k1 = 0.5 / (sdu - self.sdl)
        k2 = 0.5 / (1 - sdu)
        k3 = 0.5 / sdl
        k4 = 0.5 / sdu

        fun_cs = [lambda x: 1.,
                  lambda x: 1 - k1*(x-sdl),
                  lambda x: 0.5 - k2*(x-sdu)]
        fun_is = [lambda x: k3*x,
                  lambda x: 0.5 + k1*(x-sdl),
                  lambda x: 1]
        fun_ds = [lambda x: k4*x,
                  lambda x: k4*x,
                  lambda x: 0.5]

        fcs = np.piecewise(x, conds, fun_cs)
        fis = np.piecewise(x, conds, fun_is)
        fds = np.piecewise(x, conds, fun_ds)

        return np.array([fcs, fis, fds])

    def ict(self, scale=1):
        """
        Calculate increased charging time, using ``na``

        Parameters
        ----------
        scaler: int
            scaler for time, default unit is hour
        """
        start_time = self.tsd['ts'].iloc[0]
        end_time = self.tsd['ts'].iloc[-1]

        cts = self.ev['ts']
        masks = self.ev['ts'] < start_time
        cts[masks] = start_time

        cte = self.ev['tf']
        maske = self.ev['tf'] > end_time
        cte[maske] = end_time

        ct = cte - cts

        dsoc0 = self.ev['Pc'] * self.ev['nc'] * ct / self.ev['Q']
        dsoc = self.ev['soc'] - self.ev['soc0']
        td = (dsoc0 - dsoc) * self.ev['Q'] / self.ev['Pc'] / self.ev['nc']
        self.ev['ict'] = td * scale
        self.ev['ict'][self.ev['soc'] > self.ev['socd']] = 0
        self.ev['ict'][self.ev['ict'] < 0] = 0
        return True
