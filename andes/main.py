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
from time import sleep
from subprocess import call
from typing import Optional, Union

import andes
from andes.system import System
from andes.routines import routine_cli
from andes.utils.misc import elapsed, is_interactive
from andes.utils.paths import get_config_path, tests_root, get_log_dir
from andes.shared import coloredlogs, unittest
from andes.shared import Pool, Process

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def config_logger(stream=True,
                  file=True,
                  stream_level=logging.INFO,
                  log_file='andes.log',
                  log_path=None,
                  file_level=logging.DEBUG,
                  ):
    """
    Configure a logger for the andes package with options for a `FileHandler`
    and a `StreamHandler`. This function is called at the beginning of
    ``andes.main.main()``.

    Parameters
    ----------
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
    file_level : {10, 20, 30, 40, 50}, optional
        `FileHandler` verbosity level.
    Returns
    -------
    None

    """
    logger = logging.getLogger('andes')
    logger.setLevel(logging.DEBUG)

    if log_path is None:
        log_path = get_log_dir()

    sh_formatter_str = '%(message)s'
    if stream_level == 1:
        sh_formatter_str = '%(name)s:%(lineno)d - %(levelname)s - %(message)s'
        stream_level = 10

    sh_formatter = logging.Formatter(sh_formatter_str)
    if not len(logger.handlers):
        if stream is True:
            sh = logging.StreamHandler()
            sh.setFormatter(sh_formatter)
            sh.setLevel(stream_level)
            logger.addHandler(sh)

        # file handler for level DEBUG and up
        if file is True and (log_file is not None):
            log_full_path = os.path.join(log_path, log_file)
            fh_formatter = logging.Formatter('%(process)d: %(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh = logging.FileHandler(log_full_path)
            fh.setLevel(file_level)
            fh.setFormatter(fh_formatter)
            logger.addHandler(fh)

        globals()['logger'] = logger

    if not is_interactive():
        coloredlogs.install(logger=logger, level=stream_level, fmt=sh_formatter_str)


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

    if conf_path is None:
        logger.info('Config file does not exist. Automatically saving.')
        system = System()
        conf_path = system.save_config()

    logger.info('Editing config file {}'.format(conf_path))

    editor = ''
    if edit_config is not None:
        # use `edit_config` as default editor
        editor = edit_config
    else:
        # use the following default editors
        if platform.system() == 'Linux':
            editor = os.environ.get('EDITOR', 'vim')
        elif platform.system() == 'Darwin':
            editor = os.environ.get('EDITOR', 'vim')
        elif platform.system() == 'Windows':
            editor = 'notepad.exe'

    editor_cmd = editor.split()
    editor_cmd.append(conf_path)
    call(editor_cmd)
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
    print(f"""
    ANDES version {andes.__version__}

    Copyright (c) 2015-2020 Hantao Cui

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is included below.
    For further information, see <http://www.gnu.org/licenses/>.
    """)
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


def run_case(case, routine='pflow', profile=False, convert='', convert_all='', add_book=None,
             remove_pycapsule=False, **kwargs):
    """
    Run a single simulation case.
    """
    pr = cProfile.Profile()
    # enable profiler if requested
    if profile is True:
        pr.enable()

    system = load(case, **kwargs)
    if system is None:
        return

    # convert to the requested format
    if convert != '':
        andes.io.dump(system, convert, overwrite=None)
        return system
    # convert to xlsx with all model templates
    elif convert_all != '':
        andes.io.xlsx.write(system, system.files.dump, skip_empty=False, overwrite=None)
        return system
    # add template workbook
    elif add_book is not None:
        andes.io.xlsx.write(system, system.files.dump, skip_empty=True,
                            overwrite=True, add_book=add_book)
        return system

    if routine is not None:
        if isinstance(routine, str):
            routine = [routine]
        if 'pflow' in routine:
            routine = list(routine)
            routine.remove('pflow')

        system.PFlow.run(**kwargs)
        for r in routine:
            system.__dict__[routine_cli[r]].run(**kwargs)

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
            logger.info(f'cProfile text data written to "{system.files.prof}".')
            logger.info(f'cProfile raw data written to "{system.files.prof_raw}". View with tool `snakeviz`.')

    if remove_pycapsule is True:
        system.remove_pycapsule()

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
    logger.info(f'Working directory: "{os.getcwd()}"')
    for file in filename:
        full_paths = os.path.join(path, file)
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
        valid_cases = sorted(valid_cases)
        logger.debug('Found files: ' + pprint.pformat(valid_cases))

    return valid_cases


def set_logger_level(lg, type_to_set, level):
    """Set logging level for the given type of handler."""
    for ii, h in enumerate(lg.handlers):
        if isinstance(h, type_to_set):
            h.setLevel(level)


def find_log_path(lg):
    """Find the file paths of the FileHandlers."""
    out = []
    for ii, h in enumerate(lg.handlers):
        if isinstance(h, logging.FileHandler):
            out.append(h.baseFilename)
    return out


def _run_multiprocess_proc(cases, ncpu=os.cpu_count(), **kwargs):
    """
    Run multiprocessing with `Process`.

    Return values from `run_case` are not preserved. Always return `True` when done.
    """
    # start processes
    jobs = []
    for idx, file in enumerate(cases):
        job = Process(name=f'Process {idx:d}', target=run_case, args=(file,), kwargs=kwargs)
        jobs.append(job)
        job.start()
        start_msg = f'Process {idx:d} for "{file:s}" started.'
        print(start_msg)
        logger.debug(start_msg)
        if (idx % ncpu == ncpu - 1) or (idx == len(cases) - 1):
            sleep(0.1)
            for job in jobs:
                job.join()
            jobs = []

    return True


def _run_multiprocess_pool(cases, ncpu=os.cpu_count(), verbose=logging.INFO, **kwargs):
    """
    Run multiprocessing jobs using Pool.

    This function returns all System instances in a list, but requires longer computation time.

    Parameters
    ----------
    ncpu : int, optional = os.cpu_cout()
        Number of cpu cores to use in parallel
    mp_verbose : 10 - 50
        Verbosity level during multiprocessing
    verbose : 10, 20, 30, 40, 50
        Verbosity level outside multiprocessing
    """
    from functools import partial
    pool = Pool(ncpu)
    print("Cases are processed in the following order:")
    print('\n'.join([f'"{name}"' for name in cases]))
    ret = pool.map(partial(run_case, verbose=verbose, remove_pycapsule=True, **kwargs), cases)

    return ret


def run(filename, input_path='', verbose=20, mp_verbose=30, ncpu=os.cpu_count(), pool=False,
        cli=False, **kwargs):
    """
    Entry point to run ANDES routines.

    Parameters
    ----------
    filename : str
        file name (or pattern)
    input_path : str, optional
        input search path
    verbose : int, 10 (DEBUG), 20 (INFO), 30 (WARNING), 40 (ERROR), 50 (CRITICAL)
        Verbosity level
    mp_verbose : int
        Verbosity level for multiprocessing tasks
    ncpu : int, optional
        Number of cpu cores to use in parallel
    pool: bool, optional
        Use Pool for multiprocessing to return a list of created Systems.
    kwargs
        Other supported keyword arguments

    Other Parameters
    ----------------
    return_code : bool, optional
        Return exit code instead of system instances

    Returns
    -------
    System
        An instance

    """
    if is_interactive():
        config_logger(file=False, stream_level=verbose)

    cases = _find_cases(filename, input_path)
    system = None

    t0, _ = elapsed()
    if len(cases) == 1:
        system = run_case(cases[0], **kwargs)
    elif len(cases) > 1:

        # suppress logging output during multiprocessing
        logger.info('-> Processing {} jobs on {} CPUs.'.format(len(cases), ncpu))
        set_logger_level(logger, logging.StreamHandler, mp_verbose)
        set_logger_level(logger, logging.FileHandler, logging.DEBUG)

        kwargs['disable_pbar'] = True
        if pool is True:
            system = _run_multiprocess_pool(cases, ncpu=ncpu, verbose=verbose, mp_verbose=mp_verbose, **kwargs)
        else:
            system = _run_multiprocess_proc(cases, ncpu=ncpu, verbose=verbose, mp_verbose=mp_verbose, **kwargs)

        # restore command line output when all jobs are done
        set_logger_level(logger, logging.StreamHandler, verbose)

        log_files = find_log_path(logger)
        if len(log_files) > 0:
            log_paths = '\n'.join(log_files)
            print(f'Log saved to "{log_paths}".')

    t0, s0 = elapsed(t0)

    ex_code = 0
    if len(cases) == 1:
        ex_code = system.exit_code
    elif len(cases) >= 2:
        if isinstance(system, list):
            for s in system:
                ex_code += s.exit_code

    if len(cases) == 1:
        if ex_code == 0:
            logger.info(f'-> Single process finished in {s0}.')
        else:
            logger.error(f'-> Single process finished with error in {s0}.')
    else:
        if ex_code == 0:
            print(f'-> Multiprocessing finished in {s0}.')
        else:
            print(f'-> Multiprocessing finished with error in {s0}.')

    if cli is True:
        return ex_code
    else:
        return system


def plot(**kwargs):
    """Wrapper for the plot tool."""
    from andes.plot import tdsplot
    tdsplot(**kwargs)


def misc(edit_config='', save_config='', show_license=False, clean=True, **kwargs):
    """
    Misc functions.
    """
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
    """
    Run code generation.

    Returns
    -------
    System object
    """
    t0, _ = elapsed()
    logger.info('Numeric code preparation started...')
    system = System()
    system.prepare(quick=quick)
    _, s = elapsed(t0)
    logger.info(f'Successfully generated numerical code in {s}.')
    return system


def selftest(**kwargs):
    """
    Run unit tests.
    """
    logger.handlers[0].setLevel(logging.WARNING)
    sys.stdout = open(os.devnull, 'w')  # suppress print statements

    test_directory = tests_root()
    suite = unittest.TestLoader().discover(test_directory)
    unittest.TextTestRunner(verbosity=2).run(suite)
    sys.stdout = sys.__stdout__


def doc(attribute=None, list_supported=False, config=False, **kwargs):
    """
    Quick documentation from command-line.
    """
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
