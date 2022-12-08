import itertools
import sys
import time
from collections import OrderedDict
from tqdm import tqdm
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import pprint
import io
from contextlib import redirect_stdout

import logging
logger = logging.getLogger(__name__)


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


class DictAttr():
    """Class for attributes stored in OrderedDict"""

    def __init__(self, attr):
        """
        Base class for attribtues stored in OrderedDict

        Parameters
        ----------
        attr: OrderedDict
            Data attribute dictionary
        """
        for key, val in attr.items():
            setattr(self, key, val)
        self._dict = self.as_dict()

    def as_dict(self) -> OrderedDict:
        """
        Return the config fields and values in an ``OrderedDict``.
        """
        out = []
        for key, val in self.__dict__.items():
            if not key.startswith('_'):
                out.append((key, val))
        return OrderedDict(out)

    def __repr__(self):
        self._dict = self.as_dict()
        return pprint.pformat(self._dict)


class EVData():
    """EV data class"""

    def __init__(self, s, d, ts) -> None:
        """
        Two pandas.DataFrame are used to store EV data, one for static data,
        the other for dynamic data.

        Static data is initialized once, and dynamic data is updated at each
        time step by ``MCS.run()``.

        After initialization, Static data is not allowed to be changed.

        Parameters
        ----------
        s: pd.DataFrame
            Static data.
        d: pd.DataFrame
            Dynamic data.
        """
        self.s = s
        self.d = d
        self.ts = ts

    def __repr__(self) -> str:
        info = f'EVData: {self.s.shape[0]} EVs'
        return pprint.pformat(info)


class EVStation():
    """
    EV Station class to hold EV data, control EV status, and collecte EV info.
    """

    def __init__(self,
                 config, mcs_config,
                 ud_param, nd_param,
                 name='EVS') -> None:
        """
        Parameters
        ----------
        config: dict
            EV station configuration.
        ud_param: Dict of Dict
            Uniform distribution parameters.
        nd_param: Dict of Dict
            Normal distribution parameters.
        t: float
            Simulation start time in 24 hour.
        name: str
            EV station name.

        config
        ------
        N: int
            Number of EVs
        Ns: int
            Number of SOC intervals
        Tagc: float
            AGC time period in second
        socf: float
            SOC level that will be switched to force charging
        seed: int
            random seed
        r: float
            Ratio of time param type1 to type2, [0, 1].

        nd_param
        --------
        soci: float
            Initial SOC
        socd: float
            Demanded SOC
        ts1: float
            Start charging time 1
        ts2: float
            Start charging time 2
        tf1: float
            Finish charging time 1
        tf2: float
            Finish charging time 2
        tt: float
            Tolerance of increased charging time

        ud_param
        --------
        Pc: float
            Rated charging power
        Pd: float
            Rated discharging power
        nc: float
            Rated charging efficiency
        nd: floar
            Rated discharging efficiency
        Q: float
            Rated battery capacity
        """
        self.name = name
        # --- config ---
        self.config = DictAttr(config)
        # --- declear data variable ---
        dcols = {'s': ['ts', 'tf', 'tt', 'soc0', 'na0',
                       'soci', 'socd', 'Pc', 'Pd',
                       'nc', 'nd', 'Q'],
                 'd': ['u', 'u0', 'soc', 'c', 'lc', 'sx',
                       'na', 'ama', 'agc', 'mod']}
        sdata = pd.DataFrame()  # static data
        sdata['idx'] = [i for i in range(self.config.N)]
        sdata[[dcols['s']]] = np.nan
        ddata = pd.DataFrame()  # dynamic data
        ddata['idx'] = sdata['idx']
        ddata[[dcols['d']]] = np.nan

        # --- initialize MCS ---
        # add current timestamp `t` to config
        mcs_config = {**{'t': 0.0, 'tf': 0.0, 'socf': self.config.socf,
                         'agc': self.config.agc,
                         'ict': self.config.ict},
                      **mcs_config}
        # put data into MCS
        self.MCS = MCS(config=mcs_config, sdata=sdata, ddata=ddata)

        datas = self.MCS.data.s  # pointer to static data
        datad = self.MCS.data.d  # pointer to dynamic data

        # --- initialize data ---
        # --- 1. uniform distribution parameters ---
        ud_cols = ['Pc', 'Pd', 'nc', 'nd', 'Q']
        np.random.seed(self.config.seed)
        for col in ud_cols:
            datas[col] = np.random.uniform(size=self.config.N,
                                           low=ud_param[col]['lb'],
                                           high=ud_param[col]['ub'])
        datas['nd'] = datas['nc']  # NOTE: assumtpion: nc = nd

        # --- 2. normal distribution parameters ---
        # --- 2.1 ---
        nd_cols = ['soci', 'socd', 'tt']
        for col in nd_cols:
            a = (nd_param[col]['lb'] - nd_param[col]['mu']) / nd_param[col]['var']
            b = (nd_param[col]['ub'] - nd_param[col]['mu']) / nd_param[col]['var']
            distribution = stats.truncnorm(a, b,
                                           loc=nd_param[col]['mu'],
                                           scale=nd_param[col]['var'])
            datas[col] = distribution.rvs(self.config.N,
                                          random_state=self.config.seed)

        # --- 2.2 time parameters ---
        nd_cols = ['ts1', 'ts2', 'tf1', 'tf2']
        tparam = pd.DataFrame()
        for col in nd_cols:
            a = (nd_param[col]['lb'] - nd_param[col]['mu']) / nd_param[col]['var']
            b = (nd_param[col]['ub'] - nd_param[col]['mu']) / nd_param[col]['var']
            distribution = stats.truncnorm(a, b,
                                           loc=nd_param[col]['mu'],
                                           scale=nd_param[col]['var'])
            tparam[col] = distribution.rvs(self.config.N,
                                           random_state=self.config.seed)

        r = self.config.r  # ratio of ts1 to ts2
        tp1 = tparam[['ts1', 'tf1']].sample(n=int(self.config.N * r),
                                            random_state=self.config.seed)
        tp2 = tparam[['ts2', 'tf2']].sample(n=int(self.config.N * (1 - r)),
                                            random_state=self.config.seed)
        tp = pd.concat([tp1, tp2], axis=0).reset_index(drop=True).fillna(0)
        tp['ts'] = tp['ts1'] + tp['ts2']
        tp['tf'] = tp['tf1'] + tp['tf2']

        check = tp['ts'] > tp.tf
        row_idx = tp[check].index
        mid = tp['tf'].iloc[row_idx].values
        tp['tf'].iloc[row_idx] = tp['ts'].iloc[row_idx]
        tp['ts'].iloc[row_idx] = mid
        datas['ts'] = tp['ts']
        datas['tf'] = tp['tf']

        def _check_mem(df):
            """Check memory usage of a dataframe"""
            buffer = io.StringIO()
            with redirect_stdout(buffer):
                df.info(verbose=False, memory_usage=True, buf=buffer)
            s = buffer.getvalue()
            mem_use = s.split(' ')[-2] + ' ' + s.split(' ')[-1].strip('\n')
            return mem_use

        # --- memory save settings ---
        # --- check memory usage ---
        sd0 = _check_mem(datad)
        ss0 = _check_mem(datas)
        info_mem = f'{self.name} Memory usage:\n'\
            f'Static data: {ss0}, Dynamic data: {sd0}\n'
        info_mem_save = ''
        if self.config.memory_save:
            # TODO: mask is not correct
            mask_u = datas[(datas['ts'] > self.MCS.config.ts + self.MCS.config.th)
                           | (datas['tf'] < self.MCS.config.ts)].index
            datad.drop(mask_u, axis=0, inplace=True)
            datas.drop(mask_u, axis=0, inplace=True)
            # --- reset index ---
            datad.reset_index(drop=True, inplace=True)
            datas.reset_index(drop=True, inplace=True)
            datad['idx'] = range(len(datad))
            datas['idx'] = datad['idx']
            # --- check memory usage ---
            sd1 = _check_mem(datad)
            ss1 = _check_mem(datas)
            # --- info ---
            info_mem_save = f'Memory save is turned on, EVs out of time range '\
                f'[{self.MCS.config.ts}, {self.MCS.config.ts + self.MCS.config.th}] are dropped. \n'\
                f'Static data: {ss1}, Dynamic data: {sd1}'
        logger.warning(info_mem + info_mem_save)

        # --- 3. online status ---
        self.MCS.g_u()
        datad['u0'] = datad['u']

        # --- 4. initialize SOC ---
        # TODO: do we need to consider the AGC participation?
        # time required to charge to demanded SOC
        tr = (datas['socd'] - datas['soci']) * datas['Q'] / datas['Pc'] / datas['nc']
        # stay time
        tc = self.MCS.config.ts - datas['ts']
        tc[tc < 0] = 0  # reset negative time to 0
        # charge
        datad['soc'] = datas['soci'] + tc * datas['Pc'] * datas['nc'] / datas['Q']
        # ratio of stay/required time
        kt = tc/tr
        kt[kt < 1] = 1
        mask = kt[kt > 1].index
        # higher than required charging time, log scale higher than socd
        socp = datas['socd'] + np.log(kt) * (1 - datas['socd'])
        datad.loc[mask, 'soc'] = socp.iloc[mask]
        # reset out of range EV soc
        datad.loc[datad['soc'] >= 1, 'soc'] = 1.0
        datad.loc[datad['soc'] <= 0, 'soc'] = 0.0
        datas['soc0'] = datad['soc']
        datad['sx'] = np.ceil(datad['soc'] / (1 / self.config.Ns)) - 1

        # --- 5. initialize control signal ---
        datad['c'] = 0.0
        mask_c = datad[(datad['soc'] < datas['socd']) & datad['u'] > 0].index
        datad.loc[mask_c, 'c'] = 1.0

        # --- 6. initialize other parameters ---
        datad['agc'] = 0.0
        datad['mod'] = 0.0

        # --- 7. initialize na [number of action] ---
        # load history data of ina
        if self.config.ict:
            ina = np.genfromtxt('ev_ina_ict.csv', delimiter=',')
        else:
            ina = np.genfromtxt('ev_ina.csv', delimiter=',')

        # initialization of number of actions; truncated normal distribution;
        sx0 = np.ceil(datad['soc'] / (1 / self.config.Ns)) - 1
        # size of each sx
        sx0d = sx0.value_counts().sort_index()
        for i in sx0d.index:
            i = int(i)
            a, b = ina[i, 2], ina[i, 3]
            if a == b:
                b = a + 0.01
            pdf = stats.norm(loc=0, scale=ina[i, 1])
            res = pdf.rvs(sx0d[float(i)], random_state=self.config.seed).round(0)
            mask = datad[sx0 == i].index
            datad.loc[mask, 'na'] = ina[i, 0] * (datas['tf'].iloc[mask] - self.MCS.config.ts) + res
        mask = datad[(datad['soc'] < datas['socd']) & (datad['na'] < 0)].index
        na0 = 1000 * (datas['socd'] - datad['soc'])
        datad.loc[mask, 'na'] = na0.iloc[mask]
        # DEBUG: scale up soc [0.6, 0.7] na0
        mask = datad[(datad['soc'] < 0.7) & (datad['soc'] > 0.5)].index
        datad.loc[mask, 'na'] = na0.iloc[mask] * 10
        # for fully charged EVs, reset their na to 0
        mask = datad[(datad['soc'] >= datas['socd'])].index
        datad.loc[mask, 'na'] = 0.0
        datad['na'] = datad['na'].round(0).astype(float)
        datas['na0'] = datad['na']
        # TODO: calc number of action mileage
        datad['ama'] = 0.0  # `nama` is the number of action mileage

        # --- 8. initialize nam [max number of action] ---
        pcn = datas['Pc'] * datas['nc']
        datas['nam'] = ((datas['tf'].mean() - datas['ts'].mean() + datas['tt']) * pcn
                        - datas['socd'] * datas['Q']) / (pcn * self.config.Tagc / 3600)
        datas['nam'] = datas['nam'].round(0).astype(float)

        # --- 9. initialize lc ---
        datad['lc'] = 0.0
        if self.config.ict:
            mask = datad[(datad['na'] >= datas['nam'])].index
            datad.loc[mask, 'na'] = datas['nam'].iloc[mask]
            datad.loc[mask, 'lc'] = 1.0
        cond2 = datad['soc'] <= self.config.socf  # force charging SOC level
        datad.loc[cond2, 'lc'] = 1.0

        # --- 10. initialize MCS data ---
        self.MCS.g_ts()

        # --- 10. data dict ---
        # TODO: how to organize?

        # --- report info ---
        init_info = f'{self.name}: Initialized successfully with:\n'\
            f'Capacity: {self.config.N}, r: {self.config.r}\n'\
            + self.__repr__()
        logger.warning(init_info)

    def __repr__(self) -> str:
        # TODO: how to organize?
        datad = self.MCS.data.d
        total = datad.shape[0]
        online = int(datad['u'].sum())
        info = f'{self.name}: clock time: {self.MCS.config.ts + self.MCS.config.t / 3600}, Online: {online}, Total: {total}'
        return info

    def rctrl(self):
        """Response to control signal"""
        pass


class MCS():
    """Monte-Carlo simulation class"""

    def __init__(self, config, sdata, ddata) -> None:
        """
        Parameters
        ----------
        config: dict
            Monte-Carlo simulation configuration.
        data: EVData
            Data for Monte-Carlo simulation.

        config
        ------
        t: float
            Current timestamp in seconds.
        tf: float
            Simulation end time in seconds.
        ts: float
            Simulation start time in 24 hour.
        th: float
            Simulation tiem horizon time in hour.
            If memory_save is True, all the EVs out of the
            time horizon ``[ts, ts + th]`` will be dropped.
        h: float
            Simulation time step in second.
        no_tqdm: bool
            Disable tqdm progress bar.
        """
        self.config = DictAttr(config)

        # declear DataFrame for time series data
        ts = pd.DataFrame()
        # time in seconds
        ts['t'] = np.arange(0,
                            self.config.th * 3600 + 0.1,
                            self.config.h)
        cols = ['Pi', 'Prc', 'Ptc']
        ts[cols] = np.nan
        self.data = EVData(s=sdata, d=ddata, ts=ts)

        # --- declear info dict ---
        info = {'t': self.config.t, 'Pi': 0,
                'Prc': 0.0, 'Ptc': 0}
        self.info = DictAttr(info)

    def g_ts(self) -> True:
        """Update info into time series data"""
        datas = self.data.s
        datad = self.data.d
        # NOTE: `Ptc`, `Prc`, are converted from kW to MW, seen from the grid
        Prc = datad['agc'] * datad['u'] * datas['Pc'] * datas['nc']
        Ptc = datad['c'] * datad['u'] * datas['Pc'] * datas['nc']
        info = {'t': self.config.t, 'Pi': 0,
                'Prc': -1 * Prc.sum() * 1e-3,
                'Ptc': -1 * Ptc.sum() * 1e-3}
        self.info = DictAttr(info)
        datats = self.data.ts
        # TODO: this might be slow and wrong
        rid = datats[datats['t'] >= info['t']].index[0]
        datats.loc[rid, info.keys()] = info.values()

    def __repr__(self) -> str:
        # TODO; any other info?
        t0 = self.data.ts['t'].iloc[0]
        t1 = self.config.tf
        info = f'MCS: start from {t0}s, end at {t1}s, beginning from {self.config.ts}[H], '
        return info

    def run(self) -> bool:
        """
        Run Monte-Carlo simulation

        Parameters
        ----------
        datas: pd.DataFrame
            EV static data.
        datad: pd.DataFrame
            EV dynamic data.
        """
        # TODO: extend variable if self.config.tf > self.config.th * self.config.h
        # --- variable ---
        datas = self.data.s
        datad = self.data.d
        t0 = self.config.t  # start time of this run
        resume = t0 > 0  # resume flag, true means not start from zero
        perc_t = 0  # total percentage
        perc_add = 0  # incremental percentage
        pbar = tqdm(total=100, unit='%', file=sys.stdout,
                    disable=self.config.no_tqdm)
        # --- loop ---
        while self.config.t < self.config.tf:
            # --- computation ---
            # --- 1. update timestamp ---
            self.config.t += self.config.h

            # --- 2. update EV online status ---
            self.g_u()

            # --- 3. update control ---
            self.g_c()
            # --- 4. update EV dynamic data ---
            # --- 4.1 update soc interval and online status ---
            # charging/discharging power, kW
            datad['soc'] += datad['c'] * datas['nc'] * datas['Pc'] \
                / datas['Q'] * self.config.h / 3600
            # --- 4.2 modify outranged SoC ---
            masku = datad[datad['soc'] >= 1.0].index
            maskl = datad[datad['soc'] <= 0.0].index
            datad.loc[masku, 'soc'] = 1.0
            datad.loc[maskl, 'soc'] = 0.0

            # --- log info ---
            self.g_ts()
           # --- 2. update progress bar ---
            if resume:
                perc = 100 * self.config.h / (self.config.tf - t0)
                perc_add = 100 * t0 / self.config.tf
                resume = False  # reset resume flag
            else:
                perc = 100 * self.config.h / self.config.tf
            perc_update = perc + perc_add
            perc_update = round(perc_update, 2)
            # --- limit pbar not exceed 100 ---
            perc_update = min(100 - perc_t, perc_update)
            pbar.update(perc_update)
            perc_t += perc_update

        pbar.close()
        return True

    def g_c(self, Pi=0, is_test=False) -> bool:
        """
        Generate EV control signal.
        EV start to charge with rated power as soon as it plugs in.
        The process will not be interrupted until receive control signal
        or achieved demanmded SoC level.


        Test mode is used to build SSM A.

        Parameters
        ----------
        is_test: bool
            `True` to turn on test mode.
            test mode: TBD
        """
        # --- variable pointer ---
        datas = self.data.s
        datad = self.data.d
        if is_test:
            # --- test mode ---
            # TODO: add test mode
            return True
        if False:  # deadband
            # TODO: if control not zero, response to control signal
            # NOTE: from EVC
            pass
        else:
            # --- revise control if no signal ---
            # `CS` for low charged EVs, and set 'lc' to 1
            mask_lc = datad[(datad['soc'] <= self.config.socf) & (datad['u']) >= 1.0].index
            datad.loc[mask_lc, ['lc', 'c']] = 1.0
            datad.loc[mask_lc, 'mod'] = 1.0
            # `IS` for full EVs
            mask_full = datad[(datad['soc'] >= datas['socd'])].index
            datad.loc[mask_full, 'c'] = 0.0
            # `IS` for offline EVs
            datad['c'] = datad['c'] * datad['u']
        # TODO: is this necessary? reformatted control signal c2
        # self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        return True

    def g_u(self) -> bool:
        """
        Update EV online status.
        """
        # --- variable pointer ---
        datas = self.data.s
        datad = self.data.d
        datad['u0'] = datad['u']  # log previous online status
        # --- check time range ---
        # TODO: seems wrong, double check
        u_check1 = datas['ts'] <= self.config.ts
        u_check2 = datas['tf'] >= self.config.ts
        u_check = u_check1 & u_check2
        # --- update value ---
        datad['u'] = u_check
        datad['u'] = datad['u'].astype(float)
        return True
