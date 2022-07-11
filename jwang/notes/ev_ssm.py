"""EV Aggregator based on State Space Modeling"""

import itertools
from tqdm import tqdm
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)


def r_agc_sev(evs, us, vs):
    """
    Single EV reaction `ev.c` to ctrl signal `us` `vs`.

    evs columns:
    ['u', 'c', 'socx']
      0,    1,    2
    """
    ctrl = evs[1]
    if us[-1]*vs[-1] == 0:  # no AGC
        ctrl = evs[1]
    elif evs[0] == 1:  # online
        if (evs[1] == 0) & (vs[-1] == 1):  # response with vs1 [I to D]
            ctrl = np.random.choice(a=[-1, 0], p=[vs[evs[2]], 1-vs[evs[2]]],
                                    size=1, replace=True)[0]
        if (evs[1] == 1) & (us[-1] == 1):  # response with us1, [C to I]
            ctrl = np.random.choice(a=[0, 1], p=[us[evs[2]], 1-us[evs[2]]],
                                    size=1, replace=True)[0]
        if (evs[1] == 0) & (us[-1] == -1):  # response with us-1 [I to C]
            ctrl = np.random.choice(a=[1, 0], p=[us[evs[2]], 1-us[evs[2]]],
                                    size=1, replace=True)[0]
        if (evs[1] == -1) & (vs[-1] == -1):  # response with vs-1 [D to I]
            ctrl = np.random.choice(a=[0, -1], p=[vs[evs[2]], 1-vs[evs[2]]],
                                    size=1, replace=True)[0]
    elif evs[0] == 0:  # offline
        ctrl = 0
    return ctrl


def find_x(x, soc_intv):
    """
    Find soc interval of a single EV.
    """
    out = -1
    for idx in soc_intv.keys():
        if x > soc_intv[idx][0] and x <= soc_intv[idx][1]:
            out = idx
    return out


def update_xl(inl_input):
    """
    Update x series.
    columns:
    ['sx', 'xl', 'u', 'u0', 'c2', 'ts', 'c0', 'bd']
      0,    1,    2,   3,    4,    5    6,    7
    """
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
    if y == 0:
        return 0
    else:
        return x/y


class ev_ssm():
    """
    EV State Space Model.

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
        current time (unit: 24H)

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
    """

    def find_sx(self):
        self.ev['sx'] = self.ev['soc'].apply(lambda x: find_x(x, self.soc_intv))

    def report(self, is_report=True):
        """
        Report EVA info.
        """
        # --- EV summary info ---
        self.Q = self.ev.Q.sum()/1e3
        self.ev['wQ'] = self.ev['Q'] * self.ev['u']
        self.wQ = self.ev.wQ.sum()/1e3
        self.ev['wsoc'] = self.ev['soc'] * self.ev['Q'] * self.ev['u']
        self.wsoc = self.ev['wsoc'].sum() / self.wQ / 1e3
        self.ev['Ps'] = self.ev[['u', 'c', 'Pc', 'Pd']].apply(
            lambda x: x[0]*x[1]*x[2] if x[1] >= 0 else x[0]*x[1]*x[3], axis=1)
        self.Pcc = -1 * self.ev.Ps[self.ev.Ps > 0].sum()/1e3
        self.Pdc = -1 * self.ev.Ps[self.ev.Ps < 0].sum()/1e3
        self.Ptc = -1 * self.ev.Ps.sum()/1e3
        self.ev.drop(columns=['wQ', 'wsoc', 'Ps'], inplace=True)

        cid = self.ev[self.ev.u == 1].c.value_counts().to_dict()
        msg_c = "Ctrl: "
        for key in cid.keys():
            msg_c += f"{key}={cid[key]}; "

        if is_report:
            # --- report info ---
            msg_time = f'{self.name}: ts={np.round(self.ts, 4)}[H], {self.N} EVs, Total Q={self.Q.round(2)} MWh\n'
            msg_soc = f"Online {self.ne}, Q={self.wQ.round(2)} MWh, SoC={self.wsoc.round(4)}\n"
            msg_p = f"Power(MW): Pt={self.Ptc.round(4)}, Pc={self.Pcc.round(4)}, Pd={self.Pdc.round(4)}\n"
            logger.warning(msg_time + msg_soc + msg_p + msg_c)

    def g_ts(self, ts):
        """
        Update time and time series.
        """
        if abs(ts - self.ts) < 0.1 * self.step/3600:
            logger.warning(f"{self.name}: {ts}[H] is too close to current time={self.ts} [H]")
        else:
            self.tss.append(ts)
            ts = ts if ts < 24 else ts-24
        return ts

    def __init__(self, ts=0, N=10000, step=1, tp=100,
                 lr=0.1, lp=100, seed=None, name="EVA",
                 is_report=True):
        """
        Note:
        
        For efficienct, the EVs that are not in the range of [ts, ts+1] will be droped
        after initialization.

        Parameters
        ----------
        ts: float
            Current time in hour, [0, 24].
        N: int
            Number of EVs
        tp: int
            SSM update period (second).
        lr: float
            learning rate of updating SSM A
        lp: int
            SSM update length, data length for update A.
        step: int
            Step size in seconds.
        seed: int
            Random seed. ``None`` for random.
        """
        # --- 1. init ---
        self.name = name
        self.N = N
        self.step = step
        self.tp = tp
        self.Np = int(tp / step)  # update cycle
        self.lr = lr
        self.lp = lp
        self.seed = seed
        np.random.seed(self.seed)
        # --- 1a. uniform distribution parameters range ---
        self.ev_ufparam = dict(Ns=20,
                               Pl=5.0, Pu=7.0,
                               nl=0.88, nu=0.95,
                               Ql=20.0, Qu=30.0,
                               socl=0, socu=1)
        #  --- 1b. normal distribution parameters range ---
        self.ev_pdf_name = ['soci', 'socd', 'ts1', 'ts2', 'tf1', 'tf2']
        self.ev_pdf_data = {'mean':     [0.3,    0.8,    -6.5,  17.5,   8.9,  32.9],
                            'var':      [0.05,   0.03,   3.4,   3.4,    3.4,  3.4],
                            'lb':       [0.2,    0.7,    0.0,   5.5,    0.0,  20.9],
                            'ub':       [0.4,    0.9,    5.5,   24.0,   20.9, 24.0],
                            'info':  ['initial SoC', 'demanded SoC',
                                      'start charging time 1', 'start charging time 2',
                                      'finish charging time 1', 'finish charging time 2']}
        self.ts = ts
        self.tss = [ts]
        self.build(ts=ts)
        self.nel = [self.ne]
        self.report(is_report=is_report)

        self.Ptl = [self.Ptc]
        self.Pcl = [self.Pcc]
        self.Pdl = [self.Pdc]

        self.ep()
        self.Pr = 0  # estimated response to AGC
        self.Prl = [self.Pr]
        self.Prc = 0  # calculated response to AGC
        self.Prcl = [self.Prc]
        self.Per = 0  # error of AGC power
        self.Perl = [self.Per]
        self.y = self.ep()[0]
        self.yl = [self.y] # estiamted output power

        # --- SSM ---
        self.Pdbd = 1e-4  # Pdbd: deadband of Power

        # --- output: estimated FRC ---
        self.prumax = 0
        self.prdmax = 0
        self.g_frc()

        # --- adaptive soc coeff ---
        self.sdu = self.ev.socd.max()
        self.sdl = self.ev.socd.min()
        Qave = self.ev.Q.mean() / 1000
        ncave = self.ev.nc.mean()
        self.Th = (self.sdu - self.sdl) * Qave / (2 * self.Pave * ncave)

    def build(self, ts):
        """
        Build the ev DataFrame.

        Returns
        -------
        ev: pandas.DataFrame
            EV dataset
        """
        self.socl = self.ev_ufparam['socl']
        self.socu = self.ev_ufparam['socu']

        # --- soc charging/dicharging boundary ---
        self.sl = 0.005
        self.su = 0.995

        #  --- 1a. uniform distribution parameters range ---
        unit = self.ev_ufparam['socu']/self.ev_ufparam['Ns']
        self.soc_intv = {}
        decimal = 4
        for i in range(self.ev_ufparam['Ns']):
            intv_single = [np.around(i*unit, decimal), np.around((i+1)*unit, decimal)]
            self.soc_intv[i] = intv_single
        self.Ns = self.ev_ufparam['Ns']
        ev_pdf = pd.DataFrame(data=self.ev_pdf_data, index=self.ev_pdf_name).transpose()
        self.ev_nfparam = ev_pdf.to_dict()

        # --- 1c. generate EV dataset ---
        self.ev = pd.DataFrame()

        # data from uniform distribution
        cols = ['Pc', 'Pd', 'nc', 'nd', 'Q']
        cols_bound = {'Pc':   ['Pl', 'Pu'],
                      'Pd':   ['Pl', 'Pu'],
                      'nc':   ['nl', 'nu'],
                      'nd':   ['nl', 'nu'],
                      'Q':    ['Ql', 'Qu']}
        for col in cols:
            idxl = cols_bound[col][0]
            idxh = cols_bound[col][1]
            self.ev[col] = np.random.uniform(
                low=self.ev_ufparam[idxl],
                high=self.ev_ufparam[idxh],
                size=self.N)

        # Assumption: (1) Pd == Pc; (2) nd == nc;
        self.ev.Pd = self.ev.Pc
        self.ev.nd = self.ev.nc

        # data from normal distribution
        # soci, socd
        for col in self.ev_pdf_name:
            self.ev[col] = stats.truncnorm(
                (ev_pdf[col]['lb'] - ev_pdf[col]['mean']) / ev_pdf[col]['var'],
                (ev_pdf[col]['ub'] - ev_pdf[col]['mean']) / ev_pdf[col]['var'],
                loc=ev_pdf[col]['mean'], scale=ev_pdf[col]['var']).rvs(self.N)

        # ts1, ts2, tf1, tf2
        et = self.ev.copy()
        r1 = 0.5  # ratio of t1
        tp1 = self.ev[['ts1', 'tf1']].sample(n=int(et.shape[0]*r1), random_state=2021)
        tp2 = self.ev[['ts2', 'tf2']].sample(n=int(et.shape[0]*(1-r1)), random_state=2021)
        tp = pd.concat([tp1, tp2], axis=0).reset_index(drop=True).fillna(0)
        tp['ts'] = tp['ts1'] + tp['ts2']
        tp['tf'] = tp['tf1'] + tp['tf2']
        check = tp.ts > tp.tf
        row_idx = tp[check].index
        mid = tp.tf.iloc[row_idx].values
        tp.tf.iloc[row_idx] = tp.ts.iloc[row_idx]
        tp.ts.iloc[row_idx] = mid
        check = tp.ts > tp.tf
        self.ev['ts'] = tp['ts']
        self.ev['tf'] = tp['tf']
        self.ev['u'] = 1

        # Initialize delta power
        self.ev['dP'] = 0

        self.states = list(itertools.product([1, 0, -1], self.soc_intv.keys()))
        self.states_str = [str(s[1])+'S'+str(s[0]) for s in self.states]

        # --- update soc interval and online status ---
        # soc is initialized considering random behavior
        self.ev['soc'] = self.ev[['soci', 'ts', 'Pc', 'nc', 'Q', 'tf']].apply(
            lambda x: x[0] + (min(ts-x[1], x[5]-x[1]))*x[2]*x[3]/x[4] if ts > x[1] else x[0], axis=1)
        self.ev['soc'] = self.ev['soc'].apply(lambda x: min(x, self.socu))
        self.ev['soc'] = self.ev['soc'].apply(lambda x: max(x, self.socl))
        self.ev['bd'] = self.ev.soc.apply(lambda x: 1 if (x <= self.sl) | (x >= self.su) else 0)

        # --- ev online status: u0 as u ---
        self.ev['u0'] = 0
        self.g_u()
        self.ev['u0'] = self.ev.u

        # Initialize control signal, randomly assign c/dc
        self.ev['c'] = 1
        self.ev['c2'] = 0
        self.g_c(is_test=True, Pi=0)
        # `IS` for demanded charging SoC EVs
        self.ev['c'] = self.ev[['soc', 'c', 'socd']].apply(
            lambda x: 0 if (x[0] >= x[2]) & (x[1] == 1) else x[1], axis=1)
        self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)

        self.find_sx()
        self.g_x()
        self.r_state()

        # initialize x series
        self.ev['xl'] = [[[], [], []]] * self.N
        self.ne = self.ev.u.sum()

        ev_cols = ['u', 'u0',  'soc', 'bd', 'c', 'c2', 'c0', 'sx', 'dP', 'xl',
                   'soci', 'socd', 'Pc', 'Pd', 'nc', 'nd', 'Q', 'ts', 'tf']
        self.ev = self.ev[ev_cols]
        self.ev['agc'] = 0  # `agc` is indicator of participation of AGC
        self.ev['mod'] = 0  # `mod` is indicator of modified of over or under charge
        self.g_BCD()
        self.n_step = 1
        # --- drop EVs that are not in time range ---
        set0 = set(self.ev.index)
        set1 = set(self.ev[(self.ev.tf <= ts)].index)
        set2 = set(self.ev[(self.ev.ts >= ts + 1)].index)
        set_in = set0 - set1 - set2
        self.ev = self.ev.iloc[list(set_in)].reset_index(drop=True)
        return True

    def g_u(self):
        """
        Update online status of EV at given time.
        """
        # --- online status ---
        self.ev['u0'] = self.ev.u.astype(int)
        self.ev['u'] = (self.ev.ts <= self.ts) & (self.ev.tf >= self.ts)
        self.ev['u'] = self.ev['u'].astype(int)
        self.ne = self.ev.u.sum()
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
        self.ev['sx'] = self.ev['soc'].apply(lambda x: find_x(x, self.soc_intv))

        # --- update c2 ---
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})

        # --- output tab ---
        # TODO: consider low charged path?
        states = self.ev[['c2', 'sx', 'u']].apply(lambda x: (x[0], x[1]) if x[2] else (-1, -1), axis=1)
        res = dict(states.value_counts())
        self.xtab = pd.DataFrame(columns=range(self.Ns), index=[0, 1, 2], data=0)
        for key in res.keys():
            if key[1] > -1:
                self.xtab.loc[key[0], key[1]] = res[key]
        self.xtab.fillna(0, inplace=True)

        # TODO: use estimated ne rather than actual ne?
        self.rtab = self.xtab.div(self.ne)
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
        logger.warning(f'{self.name}: Save A to %s.' % csv)
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
        logger.warning(f'{self.name}: Load A from %s.' % csv)

    def plot_agc(self, figsize=(6, 3), style='default'):
        """
        Plot the AGC results.

        Parameters
        ----------
        figsize: tuple
            Figure size.
        style: str
            Plt style.
        """
        plt.style.use(style)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.plot(self.tss, self.Prl, label="Control")
        ax.plot(self.tss, self.Prcl, label="Response")

        ax.set_xlabel("Time [H]")
        ax.set_ylabel("Power (MW)")
        ax.set_title(f"AGC response")
        ax.set_xlim(self.tss[0], self.tss[-1])
        ax.grid()
        ax.legend()
        return fig, ax

    def plot(self, figsize=(6, 3), style='default'):
        """
        Plot the results.

        Parameters
        ----------
        figsize: tuple
            Figure size.
        style: str
            Plt style.
        """
        plt.style.use(style)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        p1 = ax.plot(self.tss, self.Ptl, label="Total")
        p2 = ax.plot(self.tss, self.Pcl, label="Charging")
        p3 = ax.plot(self.tss, self.Pdl, label="Discharging")
        ax2 = ax.twinx()
        p4 = ax2.plot(self.tss, self.nel, label='Online EVs', color='orange')
        ax2.set_ylabel("Number")

        ax.set_xlabel("Time [H]")
        ax.set_ylabel("Power (MW)")
        ax.set_title(f"{self.name}")
        ax.set_xlim(self.tss[0], self.tss[-1])
        ax.grid()

        lns = p1 + p2 + p3 + p4
        labels = [l.get_label() for l in lns]
        plt.legend(lns, labels)
        return fig, ax

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
        t_step = self.step / 3600
        ev_cols = ['u', 'u0', 'soc', 'sx', 'bd',
                   'ts', 'tf', 'c', 'c0', 'dP', 'c2']
        ev_copy = self.ev[ev_cols].copy()

        self.ev['tf'] = 23.9
        self.ev['ts'] = 0.1
        self.ev['soc'] = np.random.uniform(low=0.001, high=0.999, size=self.N)
        self.g_u()

        # initialize control signal: randomly assign C/I/D
        self.ev['c'] = np.random.choice([-1, 0, 1], p=[0.33, 0.33, 0.34], size=self.N)
        # revise
        # offline -> I
        self.ev['c'] = self.ev[['u', 'c']].apply(lambda x: x[1] if x[0] != 0 else 0, axis=1)
        # sl, DC -> C
        self.ev['c'] = self.ev[['u', 'c', 'soc']].apply(
            lambda x: 1 if (x[2] < self.sl) & (x[1] == -1) else x[1], axis=1)
        # demanded SoC, C -> I
        # self.ev['c'] = self.ev[['u', 'c', 'soc', 'socd']].apply(
        #     lambda x: 1 if (x[2] >= x[3]) & (x[1] == 1) else x[1], axis=1)

        # record control signal
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        self.ev['c0'] = self.ev['c2']
        self.ev[['c', 'c2', 'c0']] = self.ev[['c', 'c2', 'c0']].astype(int)

        # --- update x ---
        self.find_sx()
        self.g_x()
        self.g_xl()

        for t in tqdm(np.arange(self.ts+t_step, tf, t_step), desc=f'{self.name} MCS'):
            self.ts = self.g_ts(t)
            self.g_u()  # update online status
            self.g_c(is_test=True)  # update control signal
            # --- update soc interval and online status ---
            # charging/discharging power, kW
            self.ev['dP'] = self.ev[['Pc', 'Pd', 'nc', 'nd', 'c', 'u']].apply(
                lambda x: x[0]*x[2]*x[5]*x[4] if x[4] >= 0 else x[1]*x[3]*x[5]*x[4], axis=1)
            # --- update and modify SoC ---
            self.ev['soc'] = self.ev.soc + t_step * self.ev['dP'] / self.ev['Q']
            self.ev['soc'] = self.ev['soc'].apply(lambda x: x if x < self.socu else self.socu)
            self.ev['soc'] = self.ev['soc'].apply(lambda x: x if x > self.socl else self.socl)
            self.ev['soc'] = self.ev.soc + t_step * self.ev['dP'] / self.ev['Q']
            # --- boundary ---
            self.ev['bd'] = self.ev[['soc', 'socd']].apply(
                lambda x: 1 if (x[0] <= self.sl) | (x[0] >= x[1]) else 0)

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
                lambda x: 1 if (x[2] < self.sl) & (x[1] == -1) else x[1], axis=1)
            # demanded SoC, C -> I
            self.ev['c'] = self.ev[['u', 'c', 'soc', 'socd']].apply(
                lambda x: 1 if (x[2] >= x[3]) & (x[1] == 1) else x[1], axis=1)

            # --- update x ---
            self.find_sx()
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
            end time [H]
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
        t_step = self.step / 3600
        if tf - self.ts < 1e-5:
            logger.warning(f"{self.name}: end time {tf}[H] is too close to start time {self.ts}[H],"
                           "simulation will not start.")
        else:
            for t in tqdm(np.arange(self.ts, tf, t_step), desc=f'{self.name} MCS', disable=disable):
                if abs(t - self.ts) < 1e-6:
                    continue
                # --- update SSM A ---
                if self.n_step % self.Np == 0:
                    if is_updateA:
                        self.g_A(is_update=True)
                    if is_rstate:
                        self.r_state()
                        self.Per = self.Prc - self.Pr  # error of AGC response
                        # self.Per = 0  # turn off the error correction
                    self.report(is_report=False)

                self.ts = self.g_ts(t)
                self.g_u()  # update online status
                # TODO: add warning when Pi is 0
                self.g_c(Pi=Pi - self.Per, is_test=is_test)  # update control signal
                # --- update soc interval and online status ---
                # charging/discharging power, kW
                self.ev['dP'] = self.ev[['Pc', 'Pd', 'nc', 'nd', 'c', 'u']].apply(
                    lambda x: x[0]*x[2]*x[4] if x[4] >= 0 else -1*x[1]*x[3]*x[4], axis=1)
                # --- update and modify SoC ---
                self.ev['soc'] = self.ev.soc + t_step * self.ev['dP'] / self.ev['Q']
                self.ev['soc'] = self.ev['soc'].apply(lambda x: x if x < self.socu else self.socu)
                self.ev['soc'] = self.ev['soc'].apply(lambda x: x if x > self.socl else self.socl)

                # --- update x ---
                self.find_sx()
                self.g_x()

                if is_updateA:
                    self.g_xl()

                # record power
                self.report(is_report=False)
                # Actual AGC response: AGC switched power if not modified
                self.Prc += np.sum(self.ev.agc * self.ev.Pc * (1 - self.ev['mod'])) * 1e-3  # to MW
                self.Prl.append(self.Pr + self.Per)
                self.Perl.append(self.Per)
                self.Prcl.append(self.Prc)
                self.Ptl.append(self.Ptc)
                self.Pcl.append(self.Pcc)
                self.Pdl.append(self.Pdc)
                self.yl.append(self.y)
                self.nel.append(self.ne)

                self.n_step += 1

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
                rx = c * self.Ns + s
                data0.append(rx)
            data.append(data0)

        A0 = np.zeros((3*self.Ns, 3*self.Ns))
        for d in data:
            for i in range(len(d)-1):
                A0[d[i+1], d[i]] += 1

        # TODO: Consider data length limit: self.lp
        if is_update:
            n = int(self.lr / (1-self.lr))
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
        self.ev['tnow'] = self.ts
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
        # record last ctrl signal
        self.ev['c0'] = self.ev['c2']

        if is_test:
            pass
        else:
            pass
            if abs(Pi) > 1e-4:  # deadband
                self.r_agc(Pi=Pi)
            else:
                self.r_agc(Pi=0)
            # --- revise control ---
            # `CS` for low charged EVs
            self.ev['c'] = self.ev[['soc', 'c']].apply(
                lambda x: 1 if x[0] <= self.sl else x[1], axis=1)
            self.ev['mod'] = self.ev[['soc', 'c']].apply(
                lambda x: 1 if x[0] <= self.sl else 0, axis=1)
            # `CS` for just arrived EVs
            self.ev['c'] = self.ev[['u0', 'u', 'c']].apply(
                lambda x: 1 if (x[0] == 0) & (x[1] == 1) else x[2], axis=1)
            # `IS` for demanded charging [96%:x[2]] SoC EVs
            self.ev['c'] = self.ev[['soc', 'c', 'socd']].apply(
                lambda x: 0 if (x[0] >= 0.96) & (x[1] == 1) else x[1], axis=1)
            self.ev['mod'] = self.ev[['soc', 'c', 'socd']].apply(
                lambda x: 1 if (x[0] >= 0.96) & (x[1] == 1) else 0, axis=1)
        # `IS` for offline EVs
        self.ev['c'] = self.ev[['c', 'u']].apply(
            lambda x: x[0]*x[1], axis=1)
        # reformat to [0, 1, 2]
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
            raise ValueError("Matrix A is not defined")
        # --- out ---
        An = np.linalg.matrix_power(self.A, n=n)
        self.x0 = np.matmul(An, x0)
        self.etab = pd.DataFrame(data=self.x0.reshape(3, 20),
                                 columns=range(20),
                                 index=range(3))

    def reset(self, tnow=10, clean_xl=False):
        """
        Reset the pd.DataFrame ``ev`` to the initial time.

        Parameters
        ----------
        tnow: int
            the current time [H].
        clean_xl: bool
            `True` to clean `ev['xl']`, `False` to keep it.
        """
        if tnow:
            self.tss = [tnow]
            self.ts = self.tss[0]
        else:
            self.tss = [self.tss[0]]
            self.ts = self.tss[0]
        col = ['u', 'u0',  'soc', 'bd', 'c', 'c2', 'c0', 'sx', 'dP', 'xl',
               'soci', 'socd', 'Pc', 'Pd', 'nc', 'nd', 'Q', 'ts', 'tf']
        self.ev = self.ev[col]
        # TODO: c0?
        self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        self.ev['c0'] = self.ev['c2']
        self.ev['dP'] = 0
        if clean_xl:
            self.ev['xl'] = [[[], [], []]] * self.N
        self.ev['u0'] = 0
        self.g_u()
        self.ev['u0'] = self.ev.u

        self.find_sx()
        self.g_x()
        self.r_state()
        self.ne = self.ev.u.sum()
        self.n_step = 1
        logger.warning(f"{self.name}: Reset to {self.ts}[H]")
        self.report()
        self.Ptl = [self.Ptc]
        self.Pcl = [self.Pcc]
        self.Pdl = [self.Pdc]
        self.nel = [self.ne]

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
            logger.warning(f"{self.name}: 'rtab' is not available!")
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
        self.ev['Ps'] = self.ev[['Pc', 'c', 'Pd']].apply(lambda x: x[0]*x[1] if x[1] >= 0 else x[2]*x[1], axis=1)
        Pt = - self.ev['Ps'].sum()
        Pc = - self.ev['Ps'][self.ev['Ps'] > 0].sum()
        Pd = self.ev['Ps'][self.ev['Ps'] < 0].sum()
        Pl = - self.ev[self.ev.u == 1]['Pc'].sum()
        Pu = self.ev[self.ev.u == 1]['Pd'].sum()

        self.ev.drop(['Ps'], axis=1, inplace=True)

        out = [Pt, Pu, Pl, Pc, Pd]
        out = [x/1000 for x in out]  # kW to MW
        return out

    def ep(self):
        """
        Estimate output power (MW).
        Total, upper, lower

        Returns
        -------
        out: list of float
            [Pt, Pu, Pl] (MW)
        """
        # TODO: `ne` should be replaced with an recorded value
        self.Pt = np.matmul(self.ne * self.D, self.x0)[0]
        self.Pu = np.matmul(self.ne * self.Db, self.x0)[0]
        self.Pl = np.matmul(self.ne * self.Dd, self.x0)[0]
        self.Pa = self.ne * np.matmul(self.Da - self.D, self.x0)[0]
        self.Pb = self.ne * np.matmul(self.Db - self.Da, self.x0)[0]
        self.Pc = self.ne * np.matmul(self.Dc - self.D, self.x0)[0]
        self.Pd = self.ne * np.matmul(self.Dd - self.Dc, self.x0)[0]
        out = [self.Pt, self.Pu, self.Pl, self.Pa, self.Pc]
        return out

    def g_frc(self, adjust_num=False, Ne=1):
        """
        Estimate frequency regulation capacity (FRC)[MW].

        Parameters
        ----------
        adjust_num: bool
            True to adjust FRC by given number of EVs,
            False to estimate FRC by the number of online EVs at
            current time.
        Ne: int
            number of EVs to adjust FRC

        Returns
        -------
        out: list of float
            [prumax, prdmax] (MW)
        """
        self.report(is_report=False)
        self.ep()
        RU = self.Pu - self.Pt
        RD = self.Pt - self.Pl
        # --- demanded-SOC limited FRC is not that reasonable ---
        # TODO: 0.8 is the estimated soc demanded, may need revision
        # self.prumax = max(min((self.wsoc - 0.8) * self.wQ / T, RU), 0)
        # self.prdmax = min((1 - self.wsoc) * self.wQ / T, RD)
        # --- SSM FRC is not that reasonable ---
        k = 1
        if adjust_num:
            k = Ne / self.ne
        self.prumax = k * RU
        self.prdmax = k * RD
        return [self.prumax, self.prdmax]

    def g_BCD(self):
        """
        Build SSM B, C, D matrix.
        Different from the paper, the D matrix does not contain `ne`.

        The output of EVA (MW) `Pt=ne*D*x0`

        The output of mode {m} ``P{m}=ne*D{m}*x0``, where `m` in [a, b, c, d]
        """
        # --- B ---
        B1 = -1 * np.eye(self.Ns)
        B2 = np.eye(self.Ns)
        B3 = np.zeros((self.Ns, self.Ns))
        self.B = np.vstack((B1, B2, B3))

        # --- C ---
        C1 = np.zeros((self.Ns, self.Ns))
        C2 = -1 * np.eye(self.Ns)
        C3 = np.eye(self.Ns)
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
        D1 = -1 * np.ones((1, self.Ns))
        D2 = np.zeros((1, self.Ns))
        D3 = np.ones((1, self.Ns))
        self.D = Pave * np.hstack((D1, D2, D3))

        # --- D a,b,c,d ---
        D1 = np.zeros((1, self.Ns))
        D2 = np.zeros((1, self.Ns))
        D3 = np.ones((1, self.Ns))
        self.Da = Pave * np.hstack((D1, D2, D3))

        D1 = np.ones((1, self.Ns))
        D2 = np.ones((1, self.Ns))
        D3 = np.ones((1, self.Ns))
        self.Db = Pave * np.hstack((D1, D2, D3))
        self.Db[0, self.Ns] = 0  # low charged EV don't DC

        D1 = -1 * np.ones((1, self.Ns))
        D2 = np.zeros((1, self.Ns))
        D3 = np.zeros((1, self.Ns))
        self.Dc = Pave * np.hstack((D1, D2, D3))

        D1 = -1 * np.ones((1, self.Ns))
        D2 = -1 * np.ones((1, self.Ns))
        D3 = -1 * np.ones((1, self.Ns))
        self.Dd = Pave * np.hstack((D1, D2, D3))
        self.Dd[0, 2*self.Ns-1] = 0  # high charged EV don't C

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
        u: numpy.NDArray
            vector `u`, (Ns,); mode (a), (d): CS - IS
        v: numpy.NDArray
            vector `v`, (Ns,); mode (b), (c): IS - DS
        us: numpy.NDArray
            vector `us`, (Ns+1,); probability mode (a), (d): CS - IS
        vs: numpy.NDArray
            vector `vs`, (Ns+1,); probability mode (b), (c): IS - DS
        """
        if Pi >= 0:
            Pi_cap = min(Pi, self.prumax)
        elif Pi < 0:
            Pi_cap = max(Pi, -1*self.prdmax)
        # initialize output
        u = np.zeros(self.Ns)
        v = np.zeros(self.Ns)
        us = np.zeros(self.Ns+1)
        vs = np.zeros(self.Ns+1)
        # corrected control
        error = Pi_cap - self.Pr
        iter = 0
        c0 = self.ev.c.copy()
        while (abs(error) >= 0.005) & (iter < 10):
            error0 = error
            u, v, us, vs = self.g_agc(Pi_cap - self.Pr)
            self.ev.c = self.ev[['u', 'c', 'sx']].apply(lambda x: r_agc_sev(x, us, vs), axis=1)
            self.g_x()
            # --- record output ---
            self.y0 = self.ep()[0]  # self.Pt
            # TODO: modification of random traveling behavior
            self.x0 = self.x0 + np.matmul(self.B, u) + np.matmul(self.C, v)
            self.y = self.ep()[0]   # self.Pt
            # dx0 = np.matmul(self.B, u) + np.matmul(self.C, v)
            # self.Pr = np.matmul(self.D, dx0)[0]
            self.Pr += np.matmul(self.ne * self.D,
                                 np.matmul(self.B, u) + np.matmul(self.C, v))[0]
            error = Pi_cap - self.Pr
            iter += 1
            if abs(error0 - error) < 0.005:
                break
        self.ev['agc'] = c0 - self.ev.c
        self.uv = [u, v, us, vs]
        return u, v, us, vs

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

        u = np.zeros(self.Ns)
        v = np.zeros(self.Ns)
        us = np.zeros(self.Ns+1)
        vs = np.zeros(self.Ns+1)

        deadband = 1e-6

        if Pi >= deadband:  # deadband0:  # RegUp
            # --- step I ---
            ru = min(Pi, self.Pa) / (self.Pave * self.ne)  # total RegUp power
            u = np.zeros(self.Ns)   # C->I
            v = np.zeros(self.Ns)   # I->D
            for j in range(self.Ns):
                a = ru - np.sum(x[j+1: self.Ns]) - np.sum(x[j+1+self.Ns: 2*self.Ns])
                u[j] = min(max(a, 0), x[j])
                v[j] = min(max(a-x[j], 0), x[j+self.Ns])
            # --- step II ---
            for j in range(self.Ns):
                us[j] = min(safe_div(u[j], x[j]), 1)
                vs[j] = min(safe_div(v[j], x[j+self.Ns]+u[j]), 1)   # TODO: switch two times?
            # --- step III ---
            us[-1] = 1
            vs[-1] = 1

        elif Pi <= -deadband:  # RegDn
            # --- step I ---
            rv = max(Pi, self.Pc) / (self.Pave * self.ne)
            v = np.zeros(self.Ns)
            for j in range(self.Ns):
                a = rv + np.sum(x[2*self.Ns: 2*self.Ns+j])
                v[j] = max(min(a, 0), -1 * x[j+2*self.Ns])

            ru = min(Pi-self.Pc, 0) / (self.Pave * self.ne)
            u = np.zeros(self.Ns)
            for j in range(self.Ns):
                a = ru - np.sum(v[0:j]) - np.sum(u[0:j-1])
                u[j] = max(min(a, 0), -1 * x[j+self.Ns])
            # --- step II ---
            for j in range(self.Ns):
                vs[j] = min(safe_div(-1*v[j], x[j+2*self.Ns]), 1)
                us[j] = min(safe_div(-1*u[j], x[j+self.Ns]-v[j]), 1)
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
        conds = [x < self.sdl,
                 (x >= self.sdl) & (x < self.sdu),
                 x >= self.sdu]

        k1 = 0.5 / (self.sdu - self.sdl)
        k2 = 0.5 / (1 - self.sdu)
        k3 = 0.5 / self.sdl
        k4 = 0.5 / self.sdu

        fun_cs = [lambda x: 1.,
                  lambda x: 1 - k1*(x-self.sdl),
                  lambda x: 0.5 - k2*(x-self.sdu)]
        fun_is = [lambda x: k3*x,
                  lambda x: 0.5 + k1*(x-self.sdl),
                  lambda x: 1]
        fun_ds = [lambda x: k4*x,
                  lambda x: k4*x,
                  lambda x: 0.5]

        fcs = np.piecewise(x, conds, fun_cs)
        fis = np.piecewise(x, conds, fun_is)
        fds = np.piecewise(x, conds, fun_ds)

        return np.array([fcs, fis, fds])
