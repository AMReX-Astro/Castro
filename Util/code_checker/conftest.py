import os
from pathlib import Path

def find_fortran_files():
    # find Castro Fortran source files 
    castro_home = os.environ['CASTRO_HOME']
    p = Path(castro_home + '/Source')
    return p.glob(r'**/*.[fF]90')

def pytest_addoption(parser):
    parser.addoption(
        "--filename",
        action="append",
        default=find_fortran_files(),
        help="list of filenames to pass to test functions",
    )

def idfunc(argvalue):
    return '/'.join(str(argvalue).split('/')[-2:])

def pytest_generate_tests(metafunc):
    if "filename" in metafunc.fixturenames:
        metafunc.parametrize("filename", metafunc.config.getoption("filename"), ids=idfunc, scope="module")
