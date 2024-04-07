"""
CLPI model for test.
"""
from andes.core.var import ExtAlgeb
from andes.core.service import ConstService
from andes.core.block import PIDController
from andes.core.model import Model

from andes.models.experimental.olpi import OLPIData


class CLPIModel(Model):
    """
    Implementation for close-loop PI controller.
    """

    def __init__(self, system, config):
        Model.__init__(self, system, config)
        self.group = 'Experimental'
        self.flags.tds = True

        self.wd = ExtAlgeb(model='TurbineGov', src='pout', indexer=self.gov,
                           info='Generator speed deviation',
                           unit='p.u.',
                           tex_name=r'\omega_{dev}',
                           )
        self.pout = ExtAlgeb(model='TurbineGov', src='pout', indexer=self.gov,
                             tex_name='P_{out}',
                             info='Turbine governor output',
                             )
        self.pout0 = ConstService(v_str='pout',
                                  tex_name='P_{out0}',
                                  info='initial turbine governor output',
                                  )
        self.PID = PIDController(u=self.wd, kp=self.kP, ki=self.kI,
                                 kd=self.kD, Td=self.tD,
                                 tex_name='PID', info='PID', name='PID',
                                 ref=self.pout0,
                                 )
        self.pref = ExtAlgeb(indexer=self.gov,
                             tex_name='P_{ref}',
                             info='Turbine governor output',
                             model='TurbineGov',
                             src='pref',
                             e_str='u * PID_y',
                             v_str='u * PID_y',
                             )


class CLPI(OLPIData, CLPIModel):
    r"""
    Close-loop PI controller that takes Generator speed deviation as input.

    ```
        ┌────────────────────┐
        │      ki     skd    │
    u ->│kp + ─── + ───────  │ -> y
        │      s    1 + sTd  │
        └────────────────────┘
    ```
    """

    def __init__(self, system, config):
        OLPIData.__init__(self)
        CLPIModel.__init__(self, system, config)
