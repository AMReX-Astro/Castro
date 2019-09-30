"""
We should write a code checker to make sure that there are no simple errors. In particular:

when we call the EOS with something other than eos_input_rt, we need to make sure that eos_state % rho and eos_state % T are both defined (possibly with reasonable guesses) to ensure Newton is happy.

make sure all the EOS calls have aux initialized

make sure that all of the numerical constants use _rt and are defined as real(rt)

we should not assume that species are at the end of the state vector, e.g., q(QFS:) to get only species is unsafe.
"""

import re
import pytest

@pytest.fixture(scope="module")
def file_dat(filename):
    # open file
    return open(filename, 'r').read()

def test_check_rt(file_dat):
    r = re.compile(r'(\d*\.\d*[de]?\d+(?:_rt)?)')
    # d = re.compile(r'(\d*\.\d*d\d+)')
    rt = re.compile(r'(\d*\.\d*e?\d+_rt)')
    # no_d_or_rt = re.compile(r'(\d*\.\d+)')

    double_prec = re.compile(r'(double precision)')

    for l in file_dat.split('\n'):
        for m in re.finditer(r, l.split('!')[0]):
            assert re.fullmatch(rt, m.group(1)) is not None

        assert re.search(double_prec, l.split('!')[0]) is None

