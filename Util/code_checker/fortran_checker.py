"""
A code checker to make sure that there are no simple errors.
"""

import re
import pytest
import os
from pathlib import Path

def find_fortran_files():
    # find Castro Fortran source files 
    try:
        castro_home = os.environ['CASTRO_HOME']
    except KeyError:
        # this assumes this file is run from the Castro/Util/code_checker directory
        castro_home = '../..'

    p = Path(castro_home + '/Source')
    return p.glob(r'**/*.[fF]90')

def idfunc(argvalue):
    return '/'.join(str(argvalue).split('/')[-2:])

def pytest_generate_tests(metafunc):
    files = find_fortran_files()
    metafunc.parametrize('filename', files, ids=idfunc, scope="module")


def test_check_rt(filename):
    """
    make sure that all of the numerical constants use _rt and are defined as real(rt)
    """
    with open(filename, 'r') as file_dat:

        r = re.compile(r'\W(\d*\.\d*[de]?-?\d+(?:_rt)?)')
        rt = re.compile(r'(\d*\.\d*e?-?\d+_rt)')

        double_prec = re.compile(r'(double precision)')

        for l in file_dat.readlines():
            for m in re.finditer(r, l.split('!')[0]):
                assert re.fullmatch(rt, m.group(1)) is not None

            assert re.search(double_prec, l.split('!')[0]) is None

def test_check_eos_inputs(filename):
    """
    when we call the EOS with something other than eos_input_rt, we need to make sure that eos_state % rho and eos_state % T are both defined (possibly with reasonable guesses) to ensure Newton is happy.
    """

    with open(filename, 'r') as file_dat:
        r = re.compile(r'(?:subroutine|function)(.*?)end (?:subroutine|function)', re.M | re.S)

        for m in re.finditer(r, file_dat.read()):
            # does it contain an EOS call with something other than `eos_input_rt`?
            s = re.compile(r' eos\(eos_input_([a-z]{2}),\s?(\w+)\s?\)')
            function_body = m.group(1)

            for n in re.finditer(s, function_body):
                if n.group(1) != 'rt':
                    # now we need to check that rho and T are both defined
                    # NOTE: this bit is not clever enough to be able to deal with 
                    # branches. It's just checking that somewhere in the routine 
                    # prior to the EOS call, rho and T were set
                    try:
                        assert re.search(re.compile(fr'{n.group(2)}\s?%\s?rho\s*='), function_body[:n.start()]) is not None
                        assert re.search(re.compile(fr'{n.group(2)}\s?%\s?T\s*='), function_body[:n.start()]) is not None
                    except AssertionError:
                        # it might be that the state is burned, in which case 
                        # the state being initialized will have a different 
                        # name to the one having the eos called on it
                        assert re.search(re.compile(fr'%\s?rho\s*='), function_body[:n.start()]) is not None
                        assert re.search(re.compile(fr'%\s?T\s*='), function_body[:n.start()]) is not None

def test_eos_aux_initialized(filename):
    """
    make sure all the EOS calls have aux initialized
    """

    with open(filename, 'r') as file_dat:
        r = re.compile(r'(?:subroutine|function)(.*?)end (?:subroutine|function)', re.M | re.S)

        for m in re.finditer(r, file_dat.read()):
            s = re.compile(r' eos\(eos_input_([a-z]{2}),\s?(\w+)\s?\)')
            function_body = m.group(1)

            for n in re.finditer(s, function_body):
                # check that aux has been initialized
                # NOTE: this bit is not clever enough to be able to deal with 
                # branches. It's just checking that somewhere in the routine 
                # prior to the EOS call aux was set
                try:
                    assert re.search(re.compile(fr'{n.group(2)}\s?%\s?aux[ ():\w]*='), function_body[:n.start()]) is not None
                except AssertionError:
                    assert re.search(re.compile(fr'%\s?aux[ ():\w]*='), function_body[:n.start()]) is not None

def test_state_vector_species(filename):
    """
    we should not assume that species are at the end of the state vector, e.g., q(QFS:) to get only species is unsafe.
    """
    with open(filename, 'r') as file_dat:
        r = re.compile(r'(?:Q|U)FS:\s*[,)]')

        assert re.search(r, file_dat.read()) is None

def test_fortran_saved_variables(filename):
    """
    check that there are no variables intialized as 
        type :: var = value
    instead of 
        type, save :: var = value
    """

    fortran_datatypes = ['integer', 'real(rt)', 'logical', 'double precision', 'char']

    with open(filename, 'r') as file_dat:
        
        r = re.compile(r'(\w[\w \t,]+)::\s*\w+\s*=')

        for m in re.finditer(r, file_dat.read()):
            assert 'save' in m.group(1)

        # sometimes variable definitions don't include colons. We'll check that here
        r = re.compile(f'(?:{"|".join(fortran_datatypes)})([\w \t,]+)\s\w+\s*=')
        for m in re.finditer(r, file_dat.read()):
            assert 'save' in m.group(1)
