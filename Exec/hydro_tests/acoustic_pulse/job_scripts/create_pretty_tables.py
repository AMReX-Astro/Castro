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

    lines_lo = []
    found_l1 = False
    with open(file_lo) as flo:
        for line in flo:
            # skip everything until we see the L1 norm
            if line.find("L1 norm") > 0:
                found_l1 = True
                continue
            elif not found_l1:
                continue

            if line.startswith("#") or len(line.strip()) == 0:
                continue

            if line.startswith("Variable"):
                continue

            if len(line.replace(r"\\", "").split("&")) != 4:
                continue

            lines_lo.append(line.replace(r"\\", "").strip())


    lines_hi = []
    found_l1 = False
    with open(file_hi) as fhi:
        for line in fhi:
            # skip everything until we see the L1 norm
            if line.find("L1 norm") > 0:
                found_l1 = True
                continue
            elif not found_l1:
                continue

            if line.startswith("#") or len(line.strip()) == 0:
                continue

            if line.startswith("Variable"):
                continue

            if len(line.replace(r"\\", "").split("&")) != 4:
                continue

            lines_hi.append(line.replace(r"\\", "").strip())

    cd = ConvergenceData()

    for llo, lhi in zip(lines_lo, lines_hi):

        vlo, elo, o1, emed1 = llo.split("&")
        vhi, emed2, o2, ehi = lhi.split("&")

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

    good_vars = {"density": r"$\rho$",
                 "xmom": r"$\rho u$",
                 "ymom": r"$\rho v$",
                 "rho_E": r"$\rho E$",
                 "rho_e": r"$\rho e$",
                 "Temp": r"$T$",
                 "rho_He4": r"$\rho X(\isotm{He}{4})$",
                 "rho_C12": r"$\rho X(\isotm{C}{12})$",
                 "rho_O16": r"$\rho X(\isotm{O}{16})$",
                 "rho_Fe56": r"$\rho X(\isotm{Fe}{56})$"}

    # sdc4
    file_lo = "convergence.2d.lo.sdc4.out"
    file_hi = "convergence.2d.hi.sdc4.out"

    sdc4 = read_convergence(file_lo, file_hi)

    print("\n SDC 4 \n\n")

    for v in sdc4.data:
        if v.name in good_vars.keys():
            print(v.get_table_line(pretty_name=good_vars[v.name]))

