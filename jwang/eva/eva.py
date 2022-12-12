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


class DataUnit():
    """Base class for data unit"""

    def __init__(self, data, col, is_idx=True, idx=None) -> None:
        self.v = data
        self.col = col
        self.is_idx = is_idx
        self.idx = idx

    def as_df(self) -> pd.DataFrame:
        """Return data as a pandas.DataFrame"""
        if self.is_idx:
            pass
            df = pd.DataFrame()
            df['idx'] = self.idx
            df[self.col] = self.v
        else:
            df = pd.DataFrame(data=self.v, columns=self.col)
        return df

    def set(self, src, data, cond=None) -> bool:
        """
        Set value to a column.

        Parameters
        ----------
        src: str
            Column name.
        data: np.array
            Data to be set.
        cond: np.array, optional
            An array of bool, condition to set data.
            ``data`` must have the length of the length of ``True``.

            Only positions where ``cond`` is ``True`` will be set.

        Example
        -------
        ``cond = (self.v[:, col_idx] == 0)``

        ``set(src='u', data=u_data, cond=cond)``
        """
        col_idx = self.col.index(src)
        if cond is not None:
            if np.sum(cond) != len(data):
                raise ValueError(
                    f'Data shape incompatible: data length: {len(data)}, True cond: {np.sum(cond)}')
            self.v[:, col_idx][cond] = data
        else:
            self.v[:, col_idx] = data
        return True

    def get(self, src, cond=None) -> np.array:
        """
        Get the copy of value of a column

        Parameters
        ----------
        src: str
            Column name.

        cond: np.array, optional
            An array of bool, condition to access data.
            The length of ``cond`` must be the same as the length of the accessed data.

            Only positions where ``cond`` is ``True`` will be retrieved.
        """
        col_idx = self.col.index(src)
        if cond is not None:
            if len(cond) != self.v.shape[0]:
                raise ValueError(
                    f'Data shape incompatible: data length: {self.v.shape[0]}, True cond: {np.sum(cond)}')
        return self.v[:, col_idx].copy()


class EVData():
    """EV data class"""

    def __init__(self, idx, s, d, ts, dcs, dcd, tsca) -> None:
        """
        Two pandas.DataFrame are used to store EV data, one for static data,
        the other for dynamic data.

        Static data is initialized once, and dynamic data is updated at each
        time step by ``MCS.run()``.

        After initialization, Static data is not allowed to be changed.

        Parameters
        ----------
        s: np.array
            Static data.
        d: np.array
            Dynamic data.
        ts: np.array
            Time series data.
        col: Dict
            Data columns.
        """
        self.idx = idx
        self.s = DataUnit(data=s, col=[c for c in dcs.as_dict().keys()],
                          is_idx=True, idx=idx)
        self.d = DataUnit(data=d, col=[c for c in dcd.as_dict().keys()],
                          is_idx=True, idx=idx)
        self.ts = DataUnit(data=ts, col=[c for c in tsca.as_dict().keys()],
                           is_idx=False, idx=None)
        self.col = DictAttr({'s': dcs, 'd': dcd, 'ts': tsca})

    # def __repr__(self) -> str:
    #     ss = _check_mem(self.s)
    #     sd = _check_mem(self.d)
    #     st = _check_mem(self.ts)
    #     info = f'EVData: {self.s.shape[0]} EVs, Sdata: '+ss + ', Ddata: '+sd+', Tdata: '+st
    #     return info


def _check_mem(df):
    """Check memory usage of a dataframe"""
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        df.info(verbose=False, memory_usage=True, buf=buffer)
    s = buffer.getvalue()
    mem_use = s.split(' ')[-2] + ' ' + s.split(' ')[-1].strip('\n')
    return mem_use


class EVStation():
    """
    EV Station class to hold EV data, control EV status, and collecte EV info.
    """

    # def ctype(self):
    #     """Change data type"""
    #     ddata = self.MCS.data.d
    #     sdata = self.MCS.data.s
    #     ts = self.MCS.data.ts
    #     dcols = {'s': ['ts', 'tf', 'tt', 'soc0', 'na0',
    #                    'soci', 'socd', 'Pc', 'Pd',
    #                    'nc', 'nd', 'Q'],
    #              'd': ['u', 'u0', 'soc', 'c', 'lc', 'sx',
    #                    'na', 'ama', 'agc', 'mod']}
    #     sdata['idx'] = sdata['idx'].astype('int32')
    #     sdata[dcols['s']] = sdata[dcols['s']].astype('float32')
    #     ddata['idx'] = ddata['idx'].astype('int32')
    #     ddata[dcols['d']] = ddata[dcols['d']].astype('float32')
    #     ts = ts.astype('float32')
    #     return True

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
        # data columns
        dc = {'s': ['nam', 'ts', 'tf', 'tt', 'soc0', 'na0',
                    'soci', 'socd', 'Pc', 'Pd',
                    'nc', 'nd', 'Q'],
              'd': ['u', 'u0', 'soc', 'c', 'lc', 'sx',
                    'na', 'ama', 'agc', 'mod']}
        dcs = DictAttr({c: dc['s'].index(c) for c in dc['s']})
        dcd = DictAttr({c: dc['d'].index(c) for c in dc['d']})
        # data columns as nested DictAttr
        self.dca = DictAttr({'s': dcs, 'd': dcd})
        dca = self.dca  # pointer to data columns
        self.idx = np.arange(self.config.N)  # EV index
        # static data
        datas = np.full((self.config.N, len(dc['s'])),
                        fill_value=np.nan,
                        dtype=np.float64)
        # dynamic data
        datad = np.full((self.config.N, len(dc['d'])),
                        fill_value=np.nan,
                        dtype=np.float64)

        # --- initialize MCS ---
        # add current timestamp `t` to config
        mcs_config = {**{'t': 0.0, 'tf': 0.0,
                         'socf': self.config.socf,
                         'agc': self.config.agc,
                         'ict': self.config.ict},
                      **mcs_config}
        # put data into MCS
        self.MCS = MCS(config=mcs_config, idx=self.idx,
                       datas=datas, datad=datad,
                       dcs=dcs, dcd=dcd)

        dsp = self.MCS.data.s  # pointer to static data
        ddp = self.MCS.data.d  # pointer to dynamic data

        # --- initialize data ---
        # --- 1. uniform distribution parameters ---
        ud_cols = ['Pc', 'Pd', 'nc', 'nd', 'Q']
        np.random.seed(self.config.seed)
        for col in ud_cols:
            data = np.random.uniform(size=self.config.N,
                                     low=ud_param[col]['lb'],
                                     high=ud_param[col]['ub'])
            dsp.set(src=col, data=data)
        # NOTE: assumtpion: nc = nd
        dsp.set(src='nd', data=dsp.get(src='nc'))

        # --- 2. normal distribution parameters ---
        # --- 2.1 ---
        nd_cols = ['soci', 'socd', 'tt']
        for col in nd_cols:
            a = (nd_param[col]['lb'] - nd_param[col]['mu']) / nd_param[col]['var']
            b = (nd_param[col]['ub'] - nd_param[col]['mu']) / nd_param[col]['var']
            distribution = stats.truncnorm(a, b,
                                           loc=nd_param[col]['mu'],
                                           scale=nd_param[col]['var'])
            data = distribution.rvs(self.config.N,
                                    random_state=self.config.seed)
            dsp.set(src=col, data=data)

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
        dsp.set(src='ts', data=tp['ts'])
        dsp.set(src='tf', data=tp['tf'])

        # --- memory save settings ---
        if self.config.memory_save:
            mask_u1 = dsp.get(src='ts') > (self.MCS.config.ts + self.MCS.config.th)
            mask_u2 = dsp.get(src='tf') < self.MCS.config.ts
            mask = mask_u1 | mask_u2
            drop_id = np.where(mask)[0]
            for obj in [self, self.MCS.data,
                        self.MCS.data.s, self.MCS.data.d]:
                idx = getattr(obj, 'idx')
                idx = np.delete(idx, drop_id, axis=0)
                setattr(obj, 'idx', idx)
            for obj in [self.MCS.data.s, self.MCS.data.d]:
                value = getattr(obj, 'v')
                value = np.delete(value, drop_id, axis=0)
                setattr(obj, 'v', value)
            # --- info ---
            info_mem_save = f'Memory save is turned on, EVs out of time range '\
                f'[{self.MCS.config.ts}, {self.MCS.config.ts + self.MCS.config.th}] are dropped.'
            logger.warning(info_mem_save)

        # --- 3. online status ---
        self.MCS.g_u()
        ddp.set(src='u0', data=ddp.get(src='u0'))

        # --- 4. initialize SOC ---
        # TODO: do we need to consider the AGC participation?
        # time required to charge to demanded SOC
        factor1 = dsp.get(src='socd') - dsp.get(src='soci')
        factor2 = dsp.get(src='Q') / dsp.get(src='Pc') / dsp.get(src='nc')
        tr = factor1 * factor2
        # time that has charged
        tc = self.MCS.config.ts - dsp.get(src='ts')
        tc[tc < 0] = 0  # reset negative time to 0
        # charge
        dsp.set(src='soc0',
                data=dsp.get(src='soci') + tc / factor2)
        # ratio of stay/required time
        # kt = tc / tr
        # kt[kt < 1] = 1
        # mask = kt > 1
        # # higher than required charging time, log scale higher than socd
        # socp = dsp.get(src='socd') + np.log(kt) * (1 - dsp.get(src='socd'))
        # ddp.set(src='soc', data=socp[mask], cond=mask)
        # # reset out of range EV soc
        # soc = ddp.get(src='soc')
        # mask_u = soc > 1.0
        # mask_l = soc < 0.0
        # soc[mask_u] = 1.0
        # soc[mask_l] = 0.0
        # ddp.set(src='soc', data=soc)
        # sx_data = np.ceil(ddp.get(src='soc') / (1 / self.config.Ns)) - 1
        # ddp.set(src='sx', data=sx_data)

        # # --- 5. initialize control signal ---
        # ddp.set(src='c', data=0.0 * np.ones(len(self.idx)))
        # mask1 = (ddp.get(src='soc') < dsp.get(src='socd'))
        # mask2 = (ddp.get(src='u') > 0)
        # ddp.set(src='c', data=np.float64(mask1 & mask2))

        # # --- 6. initialize other parameters ---
        # ddp.set(src='agc', data=0.0 * np.ones(len(self.idx)))
        # ddp.set(src='mod', data=0.0 * np.ones(len(self.idx)))

        # # --- 7. initialize na [number of action] ---
        # # TODO: this part seems wrong
        # # load history data of ina
        # if self.config.ict:
        #     ina = np.genfromtxt('ev_ina_ict.csv', delimiter=',')
        # else:
        #     ina = np.genfromtxt('ev_ina.csv', delimiter=',')

        # # initialization of number of actions; truncated normal distribution;
        # soc = ddp.get(src='soc')
        # sx0 = np.ceil(soc / (1 / self.config.Ns)) - 1
        # # size of each sx
        # sx0 = pd.Series(sx0)
        # sx0d = sx0.value_counts().sort_index()
        # for i in sx0d.index:
        #     i = int(i)
        #     a, b = ina[i, 2], ina[i, 3]
        #     if a == b:
        #         b = a + 0.01
        #     pdf = stats.norm(loc=0, scale=ina[i, 1])
        #     res = pdf.rvs(sx0d[float(i)], random_state=self.config.seed).round(0)
        #     mask = (sx0 == i)
        #     tf = dsp.get(src='tf')
        #     na_data = ina[i, 0] * (tf[mask] - self.MCS.config.ts) + res
        #     na_data0 = ddp.get(src='na')
        #     na_data0[mask] = na_data
        #     ddp.set(src='na', data=na_data0)
        # soc = ddp.get(src='soc')
        # socd = dsp.get(src='socd')
        # na = ddp.get(src='na')
        # mask = (soc < socd) & (na < 0)
        # na0_data = 1000 * (socd - soc)
        # na0_data[mask] = na[mask]
        # ddp.set(src='na', data=na0_data)
        # # DEBUG: scale up soc [0.6, 0.7] na0
        # mask = (soc < 0.7) & (soc > 0.5)
        # na_data = ddp.get(src='na')
        # na_data[mask] = na0_data[mask] * 10
        # ddp.set(src='na', data=na_data)
        # # for fully charged EVs, reset their na to 0
        # mask = (soc >= socd)
        # na_data = ddp.get(src='na')
        # na_data[mask] = 0.0
        # ddp.set(src='na', data=na_data)
        # dsp.set(src='na0', data=na_data)
        # # TODO: calc number of action mileage
        # ddp.set(src='ama', data=0.0)  # `nama` is the number of action mileage
        # na_data = ddp.get(src='na')
        # na_data[na_data < 0] = 0.0
        # na_data = na_data.round(0).astype(float)
        # ddp.set(src='na', data=na_data)

        # # --- 8. initialize nam [max number of action] ---
        # Pc = dsp.get(src='Pc')
        # nc = dsp.get(src='nc')
        # tf = dsp.get(src='tf')
        # ts = dsp.get(src='ts')
        # tt = dsp.get(src='tt')
        # socd = dsp.get(src='socd')
        # Q = dsp.get(src='Q')
        # pcn = Pc * nc
        # nam_data = ((tf - ts + tt) * pcn - socd * Q) / (pcn * self.config.Tagc / 3600)
        # nam_data = nam_data.round(0).astype(float)
        # dsp.set(src='nam', data=nam_data)

        # # --- 9. initialize lc ---
        # ddp.set(src='lc', data=0.0 * np.ones(self.config.N))
        # if self.config.ict:
        #     na = ddp.get(src='na')
        #     nam = dsp.get(src='nam')
        #     lc = ddp.get(src='lc')
        #     mask = na >= nam
        #     na_data[mask] = nam[mask]
        #     lc_data = lc
        #     lc_data[mask] = 1.0
        #     ddp.set(src='na', data=na_data)
        #     ddp.set(src='lc', data=lc_data)
        # soc = ddp.get(src='soc')
        # mask = soc <= self.config.socf  # force charging SOC level
        # lc = ddp.get(src='lc')
        # lc_data = lc
        # lc_data[mask] = 1.0
        # ddp.set(src='lc', data=lc_data)

        # --- 10. initialize MCS data ---
        # self.MCS.g_ts()
    #     self.ctype()

    #     # --- 11. data dict ---
    #     # TODO: how to organize?
    #     # NOTE: include MCS info and online EV info

    #     # --- report info ---
    #     init_info = f'{self.name}: Initialized successfully with:\n'\
    #         f'Capacity: {self.config.N}, r: {self.config.r}\n'\
    #         + self.__repr__()
    #     logger.warning(init_info)

    # def __repr__(self) -> str:
    #     # TODO: how to organize?
    #     datad = self.MCS.data.d
    #     total = datad.shape[0]
    #     online = int(datad['u'].sum())
    #     info = f'{self.name}: clock time: {self.MCS.config.ts + self.MCS.config.t / 3600}, Online: {online}, Total: {total}'
    #     return info

    # def rctrl(self):
    #     """Response to control signal"""
    #     pass


class MCS():
    """
    Class for Monte-Carlo simulation.
    Store EV data and timeseries data.
    Control simulation.
    """

    def __init__(self, config, idx, datas, datad, dcs, dcd) -> None:
        """
        Parameters
        ----------
        config: dict
            Monte-Carlo simulation configuration.
        sdata: pd.DataFrame
            Static data
        ddata: pd.DataFrame
            Dynamic data

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

        # --- declear EV data ---
        # time series columns
        tsc = ['t', 'Pi', 'Prc', 'Ptc']
        # time series columns as DictAttr
        self.tsca = DictAttr({c: tsc.index(c) for c in tsc})
        # timestamp in seconds
        # NOTE: this might need to be extended if sim time is longer than th
        t = np.arange(0, self.config.th * 3600 + 0.1,
                      self.config.h)
        ts = np.full((len(t), len(tsc)),
                     fill_value=np.nan,
                     dtype=np.float64)
        ts[:, 0] = t
        self.data = EVData(idx=idx, s=datas, d=datad, ts=ts,
                           dcs=dcs, dcd=dcd, tsca=self.tsca)

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
        # TODO: exit_code
        return True

    def g_c(self, cvec=None, is_test=False) -> bool:
        """
        Generate EV control signal.
        EV start to charge with rated power as soon as it plugs in.
        The process will not be interrupted until receive control signal
        or achieved demanmded SoC level.


        Test mode is used to build SSM A.

        Parameters
        ----------
        cvec: np.array
            EV control vector from EVCenter
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
        if not cvec:
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
        else:
            # --- response with control vector ---
            # TODO: if control not zero, response to control signal
            # NOTE: from EVC
            pass
        # TODO: is this necessary? reformatted control signal c2
        # self.ev['c2'] = self.ev['c'].replace({1: 0, 0: 1, -1: 2})
        return True

    def g_u(self) -> bool:
        """
        Update EV online status.
        """
        # --- variable pointer ---

        dsp = self.data.s  # pointer to static data
        ddp = self.data.d  # pointer to dynamic data
        ddp.set(src='u0', data=ddp.get(src='u'))  # log previous online status
        # --- check time range ---
        # TODO: seems wrong, double check
        u_check1 = dsp.get(src='ts') <= self.config.ts
        u_check2 = dsp.get(src='tf') >= self.config.ts
        u_check = u_check1 & u_check2
        # --- update value ---
        ddp.set(src='u', data=u_check)
        return True
