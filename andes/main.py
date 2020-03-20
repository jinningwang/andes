#!/bin/bash python
import glob
import logging
import os
import io
import sys
import platform
import pprint
import cProfile
import pstats
from subprocess import call
from time import sleep
from typing import Optional, Union

import andes
from andes.system import System
from andes.utils.misc import elapsed, is_interactive
from andes.utils.paths import get_config_path, tests_root
from andes.shared import coloredlogs, Process, unittest

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def config_logger(logger=None,
                  stream=True,
                  file=False,
                  stream_level=logging.INFO,
                  log_file='andes.log',
                  log_path='',
                  file_level=logging.DEBUG,
                  ):
    """
    Configure a logger for the andes package with options for a `FileHandler`
    and a `StreamHandler`. This function is called at the beginning of
    ``andes.main.main()``.

    Parameters
    ----------
    logger
        Existing logger.
    stream : bool, optional
        Create a `StreamHandler` for `stdout` if ``True``.
        If ``False``, the handler will not be created.
    file : bool, optionsl
        True if logging to ``log_file``.
    log_file : str, optional
        Logg file name for `FileHandler`, ``'andes.log'`` by default.
        If ``None``, the `FileHandler` will not be created.
    log_path : str, optional
        Path to store the log file. By default, the path is generated by
        get_log_dir() in utils.misc.
    stream_level : {10, 20, 30, 40, 50}, optional
        `StreamHandler` verbosity level.

    Returns
    -------
    None

    """
    if not logger:
        logger = logging.getLogger('andes')
        logger.setLevel(logging.DEBUG)

    if not len(logger.handlers):
        if stream is True:
            sh_formatter = logging.Formatter('%(message)s')
            sh = logging.StreamHandler()

            sh.setFormatter(sh_formatter)
            sh.setLevel(stream_level)
            logger.addHandler(sh)

        # file handler for level DEBUG and up
        if file is True and (log_file is not None):
            log_full_path = os.path.join(log_path, log_file)
            fh_formatter = logging.Formatter(
                '%(process)d: %(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh = logging.FileHandler(log_full_path)
            fh.setLevel(file_level)
            fh.setFormatter(fh_formatter)
            logger.addHandler(fh)
            logger.debug(f'Logging to file {log_full_path}')

        globals()['logger'] = logger

    if not is_interactive():
        coloredlogs.install(logger=logger, level=stream_level, fmt='%(message)s')


def edit_conf(edit_config: Optional[Union[str, bool]] = ''):
    """
    Edit the Andes config file which occurs first in the search path.

    Parameters
    ----------
    edit_config : bool
        If ``True``, try to open up an editor and edit the config file. Otherwise returns.

    Returns
    -------
    bool
        ``True`` is a config file is found and an editor is opened. ``False`` if ``edit_config`` is False.
    """
    ret = False

    # no `edit-config` supplied
    if edit_config == '':
        return ret

    conf_path = get_config_path()

    if conf_path is not None:
        logger.info('Editing config file {}'.format(conf_path))

        if edit_config is None:
            # use the following default editors
            if platform.system() == 'Linux':
                editor = os.environ.get('EDITOR', 'gedit')
            elif platform.system() == 'Darwin':
                editor = os.environ.get('EDITOR', 'vim')
            elif platform.system() == 'Windows':
                editor = 'notepad.exe'
        else:
            # use `edit_config` as default editor
            editor = edit_config

        call([editor, conf_path])
        ret = True

    else:
        logger.info('Config file does not exist. Save config with \'andes '
                    '--save-config\'')
        ret = True

    return ret


def save_conf(config_path=None):
    """
    Save the Andes config to a file at the path specified by ``save_config``.
    The save action will not run if ``save_config = ''``.

    Parameters
    ----------
    config_path : None or str, optional, ('' by default)

        Path to the file to save the config file. If the path is an emtpy
        string, the save action will not run. Save to
        `~/.andes/andes.conf` if ``None``.

    Returns
    -------
    bool
        ``True`` is the save action is run. ``False`` otherwise.
    """
    ret = False

    # no ``--save-config ``
    if config_path == '':
        return ret

    if config_path is not None and os.path.isdir(config_path):
        config_path = os.path.join(config_path, 'andes.rc')

    ps = System()
    ps.save_config(config_path)
    ret = True

    return ret


def remove_output():
    """
    Remove the outputs generated by Andes, including power flow reports
    ``_out.txt``, time-domain list ``_out.lst`` and data ``_out.dat``,
    eigenvalue analysis report ``_eig.txt``.

    Returns
    -------
    bool
        ``True`` is the function body executes with success. ``False``
        otherwise.
    """
    found = False
    cwd = os.getcwd()

    for file in os.listdir(cwd):
        if file.endswith('_eig.txt') or \
                file.endswith('_out.txt') or \
                file.endswith('_out.lst') or \
                file.endswith('_out.npy') or \
                file.endswith('_out.csv') or \
                file.endswith('_prof.prof') or \
                file.endswith('_prof.txt'):
            found = True
            try:
                os.remove(file)
                logger.info('<{:s}> removed.'.format(file))
            except IOError:
                logger.error('Error removing file <{:s}>.'.format(file))
    if not found:
        logger.info('No output file found in the working directory.')

    return True


def print_license():
    with open(os.path.join(os.path.dirname(__file__), '..', 'LICENSE'), 'r') as f:
        print(f.read())
    return True


def load(case, **kwargs):
    """
    Load a case and set up without running. Return a system
    """
    system = System(case=case, options=kwargs)
    system.undill()

    if not andes.io.parse(system):
        return

    system.setup()
    return system


def run_case(case, routine=None, profile=False, convert='', convert_all='', add_book=None, **kwargs):
    """
    Run a single simulation case.
    """
    pr = cProfile.Profile()
    # enable profiler if requested
    if profile is True:
        pr.enable()

    system = load(case, **kwargs)

    # convert to the requested format
    if convert != '':
        andes.io.dump(system, convert)
        return system
    # convert to xlsx with all model templates
    elif convert_all != '':
        andes.io.xlsx.write(system, system.files.dump, skip_empty=False)
        return system
    # add template workbook
    elif add_book is not None:
        andes.io.xlsx.write(system, system.files.dump, skip_empty=True,
                            overwrite=True, add_book=add_book)
        return system

    system.PFlow.run()

    if routine is not None:
        routine = routine.lower()
        if routine == 'tds':
            system.TDS.run()
        elif routine == 'eig':
            system.EIG.run()

    # Disable profiler and output results
    if profile:
        pr.disable()

        if system.files.no_output:
            nlines = 40
            s = io.StringIO()
            ps = pstats.Stats(pr, stream=sys.stdout).sort_stats('cumtime')
            ps.print_stats(nlines)
            logger.info(s.getvalue())
            s.close()
        else:
            nlines = 999
            with open(system.files.prof, 'w') as s:
                ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
                ps.print_stats(nlines)
                ps.dump_stats(system.files.prof_raw)
            logger.info(f'cProfile text data written to <{system.files.prof}>.')
            logger.info(f'cProfile raw data written to <{system.files.prof_raw}. View it with \'snakeviz\'.')

    return system


def _find_cases(filename, path):
    """
    Find valid cases using the provided names and path

    Parameters
    ----------
    filename

    Returns
    -------

    """
    if len(filename) == 0:
        logger.info('info: no input file. Use \'andes run -h\' for help.')
    if isinstance(filename, str):
        filename = [filename]

    cases = []
    for file in filename:
        # use absolute path for cases which will be respected by FileMan
        full_paths = os.path.abspath(os.path.join(path, file))
        found = glob.glob(full_paths)
        if len(found) == 0:
            logger.error('error: file {} does not exist.'.format(full_paths))
        else:
            cases += found

    # remove folders and make cases unique
    unique_cases = list(set(cases))
    valid_cases = []
    for case in unique_cases:
        if os.path.isfile(case):
            valid_cases.append(case)
    if len(valid_cases):
        logger.debug('Found files: ' + pprint.pformat(valid_cases))

    return valid_cases


def _run_multiprocess(cases, ncpu, **kwargs):
    """
    Run multiprocessing jobs

    Returns
    -------

    """
    logger.info('Processing {} jobs on {} CPUs'.format(len(cases), ncpu))
    logger.handlers[0].setLevel(logging.WARNING)

    # start processes
    jobs = []
    for idx, file in enumerate(cases):
        job = Process(name=f'Process {idx:d}', target=run_case, args=(file,), kwargs=kwargs)
        jobs.append(job)
        job.start()

        start_msg = f'Process {idx:d} <{file:s}> started.'
        print(start_msg)
        logger.debug(start_msg)

        if (idx % ncpu == ncpu - 1) or (idx == len(cases) - 1):
            sleep(0.1)
            for job in jobs:
                job.join()
            jobs = []

    # restore command line output when all jobs are done
    logger.handlers[0].setLevel(logging.INFO)


def run(filename, input_path='', ncpu=1, verbose=20, **kwargs):
    if is_interactive():
        config_logger(file=False, stream_level=verbose)

    cases = _find_cases(filename, input_path)
    system = None

    t0, _ = elapsed()
    if len(cases) == 0:
        pass
    elif len(cases) == 1:
        system = run_case(cases[0], **kwargs)
    else:
        _run_multiprocess(cases, ncpu, **kwargs)

    t0, s0 = elapsed(t0)

    if len(cases) == 1:
        logger.info(f'-> Single process finished in {s0}.')
    elif len(cases) >= 2:
        logger.info(f'-> Multiple processes finished in {s0}.')

    return system


def plot(**kwargs):
    from andes.plot import tdsplot
    tdsplot(**kwargs)


def misc(edit_config='', save_config='', show_license=False, clean=True, **kwargs):
    if edit_conf(edit_config):
        return True
    if show_license:
        print_license()
        return True
    if save_config != '':
        save_conf(save_config)
        return True
    if clean is True:
        remove_output()
        return True

    logger.info('info: no option specified. Use \'andes misc -h\' for help.')


def prepare(quick=False, **kwargs):
    t0, _ = elapsed()
    logger.info('Numeric code preparation started...')
    system = System()
    system.prepare(quick=quick)
    _, s = elapsed(t0)
    logger.info(f'Successfully generated numerical code in {s}.')
    return system


def selftest(**kwargs):
    """
    Run unit tests
    """
    logger.handlers[0].setLevel(logging.WARNING)
    sys.stdout = open(os.devnull, 'w')  # suppress print statements

    test_directory = tests_root()
    suite = unittest.TestLoader().discover(test_directory)
    unittest.TextTestRunner(verbosity=2).run(suite)
    sys.stdout = sys.__stdout__


def doc(attribute=None, list_supported=False, **kwargs):
    system = System()
    if attribute is not None:
        if attribute in system.__dict__ and hasattr(system.__dict__[attribute], 'doc'):
            logger.info(system.__dict__[attribute].doc())
        else:
            logger.error(f'Model <{attribute}> does not exist.')

    elif list_supported is True:
        logger.info(system.supported_models())

    else:
        logger.info('info: no option specified. Use \'andes doc -h\' for help.')
