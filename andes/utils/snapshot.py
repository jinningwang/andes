"""
Utility functions for saving and loading snapshots.
"""

import dill
import andes


def save_ss(path, system):
    """
    Save a system with all internal states as a snapshot.
    """

    system.remove_pycapsule()

    with open(path, 'wb') as file:
        dill.dump(system, file, recurse=True)


def load_ss(path):
    """
    Load an ANDES snapshot and return a System object.
    """

    system = andes.System()  # process `pycode`
    with open(path, 'rb') as file:
        system = dill.load(file)

    system.fix_address()

    return system
