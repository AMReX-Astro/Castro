import argparse
import sys

import numpy as np

def sci_not(num):
    exp = int(np.log10(num))
    mant = num/10.0**exp

    if mant < 1:
        mant *= 10.0
        exp -= 1

    return fr"${round(mant, 3):5.3f} \times 10^{{{exp}}}$"

class Variable():
    def __init__(self, name, lo, o1, med, o2, hi):
        self.name = name
        self.lo = float(lo)
        self.o1 = float(o1)
        self.med = float(med)
        self.o2 = float(o2)
        self.hi = float(hi)

    def get_table_line(self, pretty_name=None, simple=False):
        if pretty_name is not None:
            name = pretty_name
        else:
            name = self.name

        if simple:
            _str = r" {:27}   {:15.10g}   {:5.3f}    {:15.10g}   {:5.3f}    {:15.10g}"
            return _str.format(name, self.lo, round(self.o1, 3), self.med, round(self.o2, 3), self.hi)

        else:
            _str = r" {:27} & {:23} & {:5.3f}  & {:23} & {:5.3f}  & {:23} \\"
            return _str.format(name, sci_not(self.lo), round(self.o1, 3), sci_not(self.med), round(self.o2, 3), sci_not(self.hi))

class ConvergenceData():
    def __init__(self):
        self.data = []

    def add_variable(self, name, lo, order1, med, order2, hi):
        self.data.append(Variable(name, lo, order1, med, order2, hi))

def read_convergence(file_lo, file_hi):

    # we'll wait until we find the L1 data

    lines_lo = []
    found_l1 = False
    with open(file_lo) as flo:
        for line in flo:
            if "L1 norm" in line:
                found_l1 = True
                continue
            if not found_l1:
                continue
            # value data lines have 4 columns
            if len(line.split()) == 4:
                lines_lo.append(line.strip())

    lines_hi = []
    found_l1 = False
    with open(file_hi) as fhi:
        for line in fhi:
            if "L1 norm" in line:
                found_l1 = True
                continue
            if not found_l1:
                continue
            # value data lines have 4 columns
            if len(line.split()) == 4:
                lines_hi.append(line.strip())

    cd = ConvergenceData()

    for llo, lhi in zip(lines_lo, lines_hi):

        vlo, elo, o1, emed1 = llo.split()
        vhi, emed2, o2, ehi = lhi.split()

        if "---" in o1 or "---" in o2:
            print(f"skipping {vlo}")
            continue

        if vlo != vhi:
            sys.exit("error: variable mismatch")

        if emed1.strip() != emed2.strip():
            print(emed1, emed2)
            sys.exit("error: error mismatch")

        cd.add_variable(vlo, elo, o1, emed1, o2, ehi)

    return cd

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--simple', action="store_true", help='no latex output')
    parser.add_argument("lofile", type=str, nargs=1,
                        help="name of the low resolution convergence output file")
    parser.add_argument("hifile", type=str, nargs=1,
                        help="name of the high resolution convergence output file")

    args = parser.parse_args()

    good_vars = {"density": r"$\rho$",
                 "xmom": r"$\rho u$",
                 "ymom": r"$\rho v$",
                 "rho_E": r"$\rho E$",
                 "rho_e": r"$\rho e$",
                 "Temp": r"$T$",
                 "rho_He4": r"$\rho X(\isotm{He}{4})$",
                 "rho_C12": r"$\rho X(\isotm{C}{12})$",
                 "rho_O16": r"$\rho X(\isotm{O}{16})$",
                 "rho_Fe56": r"$\rho X(\isotm{Fe}{56})$",
                 "rho_Ye": r"$\rho Y_e$",
                 "rho_abar": r"$\rho \bar{A}$",
                 "rho_bea": r"$\rho (B/A)$"
                 }

    # sdc4
    file_lo = args.lofile[0]
    file_hi = args.hifile[0]

    sdc4 = read_convergence(file_lo, file_hi)

    for v in sdc4.data:
        if v.name in good_vars.keys():
            print(v.get_table_line(pretty_name=good_vars[v.name], simple=args.simple))

