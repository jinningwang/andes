import numpy as np

from andes.core import (Algeb, ConstService, ExtAlgeb, ExtParam, ExtService,
                        HardLimiter, IdxParam, Lag, LeadLag, NumParam,)
from andes.core.block import IntegratorAntiWindup, Piecewise
from andes.core.service import (FlagValue, InitChecker, NumSelect, ParamCalc,
                                PostInitService,)
from andes.models.governor.tgbase import TGBase, TGBaseData


class IEEEG1Data(TGBaseData):

    def __init__(self):
        TGBaseData.__init__(self)

        self.syn2 = IdxParam(model='SynGen',
                             info='Optional SynGen idx',
                             )
        self.K = NumParam(default=20, tex_name='K',
                          info='Gain (1/R) in mach. base',
                          unit='p.u. (power)',
                          power=True,
                          vrange=(5, 30),
                          )
        self.T1 = NumParam(default=1, tex_name='T_1',
                           info='Gov. lag time const.',
                           vrange=(0, 5),
                           )
        self.T2 = NumParam(default=1, tex_name='T_2',
                           info='Gov. lead time const.',
                           vrange=(0, 10),
                           )
        self.T3 = NumParam(default=0.1, tex_name='T_3',
                           info='Valve controller time const.',
                           vrange=(0.04, 1),
                           )
        # "UO" is "U" and capitalized "o" character
        self.UO = NumParam(default=0.1, tex_name='U_o',
                           info='Max. valve opening rate',
                           unit='p.u./sec', vrange=(0.01, 0.3),
                           )
        self.UC = NumParam(default=-0.1, tex_name='U_c',
                           info='Max. valve closing rate',
                           unit='p.u./sec', vrange=(-0.3, 0),
                           )
        self.PMAX = NumParam(default=5, tex_name='P_{MAX}',
                             info='Max. turbine power',
                             vrange=(0.5, 2), power=True,
                             )
        self.PMIN = NumParam(default=0., tex_name='P_{MIN}',
                             info='Min. turbine power',
                             vrange=(0.0, 0.5), power=True,
                             )

        self.T4 = NumParam(default=0.4, tex_name='T_4',
                           info='Inlet piping/steam bowl time constant',
                           vrange=(0, 1.0),
                           )
        self.K1 = NumParam(default=0.5, tex_name='K_1',
                           info='Fraction of power from HP',
                           vrange=(0, 1.0),
                           )
        self.K2 = NumParam(default=0, tex_name='K_2',
                           info='Fraction of power from LP',
                           vrange=(0,),
                           )
        self.T5 = NumParam(default=8, tex_name='T_5',
                           info='Time constant of 2nd boiler pass',
                           vrange=(0, 10),
                           )
        self.K3 = NumParam(default=0.5, tex_name='K_3',
                           info='Fraction of HP shaft power after 2nd boiler pass',
                           vrange=(0, 0.5),
                           )
        self.K4 = NumParam(default=0.0, tex_name='K_4',
                           info='Fraction of LP shaft power after 2nd boiler pass',
                           vrange=(0,),
                           )

        self.T6 = NumParam(default=0.5, tex_name='T_6',
                           info='Time constant of 3rd boiler pass',
                           vrange=(0, 10),
                           )
        self.K5 = NumParam(default=0.0, tex_name='K_5',
                           info='Fraction of HP shaft power after 3rd boiler pass',
                           vrange=(0, 0.35),
                           )
        self.K6 = NumParam(default=0, tex_name='K_6',
                           info='Fraction of LP shaft power after 3rd boiler pass',
                           vrange=(0, 0.55),
                           )

        self.T7 = NumParam(default=0.05, tex_name='T_7',
                           info='Time constant of 4th boiler pass',
                           vrange=(0, 10),
                           )
        self.K7 = NumParam(default=0, tex_name='K_7',
                           info='Fraction of HP shaft power after 4th boiler pass',
                           vrange=(0, 0.3),
                           )
        self.K8 = NumParam(default=0, tex_name='K_8',
                           info='Fraction of LP shaft power after 4th boiler pass',
                           vrange=(0, 0.3),
                           )


class IEEEG1SpeedControl:
    def __init__(self):
        # check if K1-K8 sums up to 1
        self._sumK18 = ConstService(v_str='K1+K2+K3+K4+K5+K6+K7+K8',
                                    info='summation of K1-K8',
                                    tex_name=r"\sum_{i=1}^8 K_i"
                                    )

        self._Kcoeff = ConstService(v_str='1/_sumK18',
                                    info='normalization factor to be multiplied to K1-K8',
                                    tex_name='K_{coeff}',
                                    )
        self.K1n = ConstService(v_str='K1 * _Kcoeff',
                                info='normalized K1',
                                tex_name='K_{1n}',
                                )
        self.K2n = ConstService(v_str='K2 * _Kcoeff',
                                info='normalized K2',
                                tex_name='K_{2n}',
                                )
        self.K3n = ConstService(v_str='K3 * _Kcoeff',
                                info='normalized K3',
                                tex_name='K_{3n}',
                                )
        self.K4n = ConstService(v_str='K4 * _Kcoeff',
                                info='normalized K4',
                                tex_name='K_{4n}',
                                )
        self.K5n = ConstService(v_str='K5 * _Kcoeff',
                                info='normalized K5',
                                tex_name='K_{5n}',
                                )
        self.K6n = ConstService(v_str='K6 * _Kcoeff',
                                info='normalized K6',
                                tex_name='K_{6n}',
                                )
        self.K7n = ConstService(v_str='K7 * _Kcoeff',
                                info='normalized K7',
                                tex_name='K_{7n}',
                                )
        self.K8n = ConstService(v_str='K8 * _Kcoeff',
                                info='normalized K8',
                                tex_name='K_{8n}',
                                )

        # check if  `tm0 * (K2 + k4 + K6 + K8) = tm02 *(K1 + K3 + K5 + K7)
        self._tm0K2 = PostInitService(info='mul of tm0 and (K2n+K4n+K6n+K8n)',
                                      v_str='zsyn2*tm0*(K2n + K4n + K6n + K8n)',
                                      )
        self._tm02K1 = PostInitService(info='mul of tm02 and (K1n+K3n+K5n+K7n)',
                                       v_str='tm02*(K1n + K3n + K5n + K7n)',
                                       )
        self._Pc = InitChecker(u=self._tm0K2,
                               info='proportionality of tm0 and tm02',
                               equal=self._tm02K1,
                               )

        self.Sg2 = ExtParam(src='Sn',
                            model='SynGen',
                            indexer=self.syn2,
                            allow_none=True,
                            default=0.0,
                            tex_name='S_{n2}',
                            info='Rated power of Syn2',
                            unit='MVA',
                            export=False,
                            )
        self.Sg12 = ParamCalc(self.Sg, self.Sg2, func=np.add,
                              tex_name="S_{g12}",
                              info='Sum of generator power ratings',
                              )
        self.Sn = NumSelect(self.Tn,
                            fallback=self.Sg12,
                            tex_name='S_n',
                            info='Turbine or Gen rating',
                            )

        self.zsyn2 = FlagValue(self.syn2,
                               value=None,
                               tex_name='z_{syn2}',
                               info='Exist flags for syn2',
                               )

        self.tm02 = ExtService(src='tm',
                               model='SynGen',
                               indexer=self.syn2,
                               tex_name=r'\tau_{m02}',
                               info='Initial mechanical input of syn2',
                               allow_none=True,
                               default=0.0,
                               )
        self.tm012 = ConstService(info='total turbine power',
                                  v_str='tm0 + tm02',
                                  )

        # Note: the following applies `zsyn2` to disable the syn2
        self.tm2 = ExtAlgeb(src='tm',
                            model='SynGen',
                            indexer=self.syn2,
                            allow_none=True,
                            tex_name=r'\tau_{m2}',
                            e_str='zsyn2 * ue * (PLP - tm02)',
                            info='Mechanical power to syn2',
                            ename='tm2',
                            tex_ename=r'\tau_{m2}',
                            )

        self.wd = Algeb(info='Generator under speed',
                        unit='p.u.',
                        tex_name=r'\omega_{dev}',
                        v_str='0',
                        e_str='ue * (wref - omega) - wd',
                        )

        self.LL = LeadLag(u=self.wd,
                          T1=self.T2,
                          T2=self.T1,
                          K=self.K,
                          info='Signal conditioning for wd',
                          )

        # `P0` == `tm0`
        self.vs = Algeb(info='Valve speed',
                        tex_name='V_s',
                        v_str='0',
                        e_str='ue * (LL_y + v0 + paux - IAW_y) / T3 - vs',
                        )

        self.HL = HardLimiter(u=self.vs,
                              lower=self.UC,
                              upper=self.UO,
                              info='Limiter on valve speed',
                              )

        self.vsl = Algeb(info='Valve move speed after limiter',
                         tex_name='V_{sl}',
                         v_str='vs * HL_zi + UC * HL_zl + UO * HL_zu',
                         e_str='vs * HL_zi + UC * HL_zl + UO * HL_zu - vsl',
                         )

        self.v0 = ConstService(info='Initial valve position')

        self.IAW = IntegratorAntiWindup(u=self.vsl,
                                        T=1,
                                        K=1,
                                        y0=self.v0,
                                        lower=self.PMIN,
                                        upper=self.PMAX,
                                        info='Valve position integrator')


class IEEEG1ValvePosition:
    def __init__(self):

        self.v0.v_str = 'tm012'

        self.GV = Algeb(info='steam flow',
                        tex_name='G_{V}',
                        v_str='tm012',
                        e_str='IAW_y - GV')

        self.L4 = Lag(u=self.GV, T=self.T4, K=1,
                      info='first process',
                      )


class IEEEG1Turbine:
    def __init__(self):

        self.L5 = Lag(u=self.L4_y, T=self.T5, K=1,
                      info='second (reheat) process',
                      )

        self.L6 = Lag(u=self.L5_y, T=self.T6, K=1,
                      info='third process',
                      )

        self.L7 = Lag(u=self.L6_y, T=self.T7, K=1,
                      info='fourth (second reheat) process',
                      )

        self.PHP = Algeb(info='HP output',
                         tex_name='P_{HP}',
                         v_str='ue * (K1n*L4_y + K3n*L5_y + K5n*L6_y + K7n*L7_y)',
                         e_str='ue * (K1n*L4_y + K3n*L5_y + K5n*L6_y + K7n*L7_y) - PHP',
                         )

        self.PLP = Algeb(info='LP output',
                         tex_name='P_{LP}',
                         v_str='ue * (K2n*L4_y + K4n*L5_y + K6n*L6_y + K8n*L7_y)',
                         e_str='ue * (K2n*L4_y + K4n*L5_y + K6n*L6_y + K8n*L7_y) - PLP',
                         )

        self.pout.e_str = 'ue * PHP - pout'


class IEEEG1Model(TGBase):
    def __init__(self, system, config):
        TGBase.__init__(self, system, config, add_sn=False)
        IEEEG1SpeedControl.__init__(self)
        IEEEG1ValvePosition.__init__(self)
        IEEEG1Turbine.__init__(self)


class IEEEG1(IEEEG1Data, IEEEG1Model):
    """
    IEEE Type 1 Speed-Governing Model.

    If only one generator is connected, its `idx` must be given to `syn`, and
    `syn2` must be left blank. Each generator must provide data in its `Sn`
    base.

    `syn` is connected to the high-pressure output (PHP) and the optional `syn2`
    is connected to the low- pressure output (PLP).

    The speed deviation of generator 1 (syn) is measured. If the turbine rating
    `Tn` is not specified, the sum of `Sn` of all connected generators will be
    used.

    Normally, K1 + K2 + ... + K8 = 1.0. If the second generator is not
    connected, K1 + K3 + K5 + K7 = 1, and K2 + K4 + K6 + K8 = 0. If K1 to K8 do
    not sum up to 1.0, they will be normalized. The normalized parameters are
    called ``K1n`` through ``K8n``.

    Developer Notes
    ----------------
    After a refactoring, the `IEEEG1Model` was separated into `IEEEG1SpeedControl`,
    `IEEEG1ValvePosition`, and `IEEEG1Turbine`.

    In `IEEEG1ValvePosition`, an Algeb `GV` is developed to represent the dynamic
    of valve position to steam flow.
    Remember to define the initial valve position `v0`, which is left intentionally
    black in `IEEEG1SpeedControl`.

    In `IEEEG1`, valve position to steam flow is considered linear, and the flow
    and the valve position maintain a 1:1 linear ratio.

    To model the nonlinear process, the equation of `GV` can be replaced with other
    nonlinear equations.
    """

    def __init__(self, system, config):
        IEEEG1Data.__init__(self)
        IEEEG1Model.__init__(self, system, config)


class IEEEG1ValvePositionPL:
    def __init__(self):

        self.GV = Piecewise(u=self.IAW_y,
                            points=('PMIN', 'PMAX'),
                            funs=('PMIN',
                                  'IAW_y',
                                  'PMAX'),
                            info='steam flow',
                            tex_name='G_{V}',
                            )
        self.GV.y.v_iter = self.GV.y.e_str

        self.v0.v_str = 'tm012'

        self.L4 = Lag(u=self.GV_y, T=self.T4, K=1,
                      info='first process',
                      )


class IEEEG1PLModel(TGBase):
    def __init__(self, system, config):
        TGBase.__init__(self, system, config, add_sn=False)
        IEEEG1SpeedControl.__init__(self)
        IEEEG1ValvePositionPL.__init__(self)
        IEEEG1Turbine.__init__(self)


class IEEEG1PL(IEEEG1Data, IEEEG1PLModel):
    def __init__(self, system, config):
        IEEEG1Data.__init__(self)
        IEEEG1PLModel.__init__(self, system, config)


class IEEEG1ValvePositionNL:
    def __init__(self):

        k = 0.8

        self.GV = Algeb(info='steam flow',
                        tex_name='G_{V}',
                        v_str='tm012',
                        e_str=f'IAW_y * {k} - GV')

        self.v0.v_str = f'tm012 / {k}'

        self.L4 = Lag(u=self.GV, T=self.T4, K=1,
                      info='first process',
                      )


class IEEEG1NLModel(TGBase):
    """
    In this implementation, initialization issue still exists, just like in IEEG1PW.
    After TDS.init(), IAW.y != IAW.y0.

    It is suspected that it will not take effect if we define IAW.y0 outside of the
    definition of IAW.
    """

    def __init__(self, system, config):
        TGBase.__init__(self, system, config, add_sn=False)
        IEEEG1SpeedControl.__init__(self)
        IEEEG1ValvePositionNL.__init__(self)
        IEEEG1Turbine.__init__(self)


class IEEEG1NL(IEEEG1Data, IEEEG1NLModel):
    def __init__(self, system, config):
        IEEEG1Data.__init__(self)
        IEEEG1NLModel.__init__(self, system, config)


class IEEEG1PWData(IEEEG1Data):
    def __init__(self):
        IEEEG1Data.__init__(self)

        # Define Gv1-Gv6 and Pgv1-Pgv6
        for i in range(1, 7):
            setattr(self, f'Gv{i}', NumParam(
                info=f'Gate value-steam flow pair (point {i}), nominal gate value',
                tex_name=f'G_{{v{i}}}',
                default=None,
                power=True,
            ))
            setattr(self, f'Pgv{i}', NumParam(
                info=f'Gate value-steam flow pair (point {i}), nominal steam flow',
                tex_name=f'P_{{gv{i}}}',
                unit='p.u.',
                default=None,
                power=True,
            ))


class IEEEG1PWModel(IEEEG1Model):
    def __init__(self, system, config):
        IEEEG1Model.__init__(self, system, config)

        # Define Kgp1-Kgp6 and Pgv1p-Pgv6p. Is this a bad practice?
        self.Kgp1 = ConstService(
            v_str='(Pgv1 - PMIN) / (Gv1 - 0)',
            info='Gain of gate value-steam flow pair (point 1)',
            tex_name='K_{gp1}',
        )
        for i in range(2, 7):
            setattr(self, f'Kgp{i}', ConstService(
                v_str=f'(Pgv{i} - Pgv{i-1}) / (Gv{i} - Gv{i-1})',
                info=f'Gain of gate value-steam flow pair (point {i})',
                tex_name=f'K_{{gp{i}}}',
            ))

        self.GV = Piecewise(u=self.IAW_y,
                            points=('PMIN', 'Gv1', 'Gv2', 'Gv3', 'Gv4', 'Gv5', 'Gv6'),
                            funs=('PMIN',
                                  '(IAW_y - 0) * Kgp1 + 0',
                                  '(IAW_y - Gv1) * Kgp2 + Pgv1',
                                  '(IAW_y - Gv2) * Kgp3 + Pgv2',
                                  '(IAW_y - Gv3) * Kgp4 + Pgv3',
                                  '(IAW_y - Gv4) * Kgp5 + Pgv4',
                                  '(IAW_y - Gv5) * Kgp6 + Pgv5',
                                  '(IAW_y - Gv6) * Kgp6 + Pgv6',
                                  'PMAX'),
                            tex_name='G_{V}',
                            info='steam flow',
                            )
        self.GV.y.v_iter = self.GV.y.e_str

        self.L4.u = self.GV_y


class IEEEG1PW(IEEEG1):
    """
    UNDER DEVELOPMENT, DO NOT USE!

    IEEE Type 1 Speed-Governing Model in PowerWorld.

    The IEEEG1 implementation in PowerWorld and GE PSLF models the gate-steam as a nonlinear
    process.

    In ANDES implementation, the nonlinear gate-steam `GV` is represented as
    a piecewise linear function, defined by six points: (Gv1, Pgv1), ..., (Gv6, Pgv6).

    If left unspecified, these points will be linearly interpolated between `PMIN` and `PMAX`,
    resulting in a straight line.

    References:

    [1] PowerWorld, Governor IEEEG1, IEEEG1D, and IEEEG1_GE.
    Available at:
    https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20IEEEG1,%20IEEEG1D%20and%20IEEEG1_GE.htm

    Developer Notes
    ----------------
    In this implementation, the initialization process has not been fully tuned, and several notes
    are left here for future reference:
    1. When Gv1-Gv6, Pgv1-Pgv6 are linearly given, IEEEG1PW works as expected.
    2. When they are set as a nonlinear function (e.g., sqrt), initialization may encounter issues:
       - L4.u != L4.y, with no error message in TDS.init()
       - Cannot find a value for IAW_y different from tm012 that satisfies GV.y.v_iter
    """

    def __init__(self, system, config):
        IEEEG1PWData.__init__(self)
        IEEEG1PWModel.__init__(self, system, config)
