import sys

import numpy as np

def sci_not(num):
    exp = int(np.log10(num))
    mant = num/10.0**exp

    if mant < 1:
        mant *= 10.0
        exp -= 1

    return r"${:5.3f} \times 10^{{{}}}$".format(round(mant, 3), exp)

class Variable():
    def __init__(self, name, lo, o1, med, o2, hi):
        self.name = name
        self.lo = float(lo)
        self.o1 = float(o1)
        self.med = float(med)
        self.o2 = float(o2)
        self.hi = float(hi)

    def get_table_line(self, pretty_name=None):
        if pretty_name is not None:
            name = pretty_name
        else:
            name = self.name

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
    with open(file_lo, "r") as flo:
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
    with open(file_hi, "r") as fhi:
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
            print("skipping {}".format(vlo))
            continue

        if vlo != vhi:
            sys.exit("error: variable mismatch")

        if emed1.strip() != emed2.strip():
            print(emed1, emed2)
            sys.exit("error: error mismatch")

        cd.add_variable(vlo, elo, o1, emed1, o2, ehi)

    return cd

if __name__ == "__main__":

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
    file_lo = sys.argv[-2]
    file_hi = sys.argv[-1]

    sdc4 = read_convergence(file_lo, file_hi)

    for v in sdc4.data:
        if v.name in good_vars.keys():
            print(v.get_table_line(pretty_name=good_vars[v.name]))

