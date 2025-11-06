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
    def __init__(self, name, lo, o1, med, o2, hi, o3=None, vhi=None):
        self.name = name
        self.lo = float(lo)
        self.o1 = float(o1)
        self.med = float(med)
        self.o2 = float(o2)
        self.hi = float(hi)
        if o3 is not None:
            self.o3 = float(o3)
        else:
            self.o3 = None
        if vhi is not None:
            self.vhi = float(vhi)
        else:
            self.vhi = None

    def get_table_line(self, pretty_name=None, simple=False):
        if pretty_name is not None:
            name = pretty_name
        else:
            name = self.name

        if simple:
            if self.o3 is None:
                return rf" {name:27}   {self.lo:14.10g}   {round(self.o1, 3):5.3f}    {self.med:14.10g}   {round(self.o2, 3):5.3f}    {self.hi:14.10g}"
            else:
                return rf" {name:27}   {self.lo:14.10g}   {round(self.o1, 3):5.3f}    {self.med:14.10g}   {round(self.o2, 3):5.3f}    {self.hi:14.10g}    {round(self.o3, 3):5.3f}    {self.vhi:14.10g}"

        else:
            if self.o3 is None:
                return rf" {name:27} & {sci_not(self.lo):23} & {round(self.o1, 3):5.3f}  & {sci_not(self.med):23} & {round(self.o2, 3):5.3f}  & {sci_not(self.hi):23} \\"
            else:
                return rf" {name:27} & {sci_not(self.lo):23} & {round(self.o1, 3):5.3f}  & {sci_not(self.med):23} & {round(self.o2, 3):5.3f}  & {sci_not(self.hi):23}  & {round(self.o3, 3):5.3f}  & {sci_not(self.vhi):23} \\"

class ConvergenceData2():
    def __init__(self):
        self.data = []

    def add_variable(self, name, lo, order1, med, order2, hi):
        self.data.append(Variable(name, lo, order1, med, order2, hi))

class ConvergenceData3():
    def __init__(self):
        self.data = []

    def add_variable(self, name, lo, order1, med, order2, hi, order3, vhi):
        self.data.append(Variable(name, lo, order1, med, order2, hi, order3, vhi))

def read_convergence(file_lo, file_hi, file_vhi):

    # we'll wait until we find the L1 data

    lines_lo = []
    lines_hi = []
    lines_vhi = []

    fdata = [(lines_lo, file_lo), (lines_hi, file_hi)]
    if file_vhi is not None:
        fdata.append((lines_vhi, file_vhi))

    for lines, filec in fdata:
        found_l1 = False
        with open(filec) as fc:
            for line in fc:
                if "L1 norm" in line:
                    found_l1 = True
                    continue
                if not found_l1:
                    continue
                # value data lines have 4 columns
                if len(line.split()) == 4:
                    lines.append(line.strip())

    if file_vhi is None:

        cd = ConvergenceData2()

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

    else:

        cd = ConvergenceData3()

        for llo, lhi, lvhi in zip(lines_lo, lines_hi, lines_vhi):

            vlo, elo, o1, emed1 = llo.split()
            vhi, emed2, o2, ehi1 = lhi.split()
            vvhi, ehi2, o3, evhi = lvhi.split()

            if "---" in o1 or "---" in o2 or "---" in o3:
                print(f"skipping {vlo}")
                continue

            if vlo != vhi or vlo != vvhi:
                sys.exit("error: variable mismatch")

            if emed1.strip() != emed2.strip() or ehi1.strip() != ehi2.strip():
                print(emed1, emed2, ehi1, ehi2)
                print(llo)
                print(lhi)
                print(lvhi)
                sys.exit("error: error mismatch")

            cd.add_variable(vlo, elo, o1, emed1, o2, ehi1, o3, evhi)

    return cd

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--simple', action="store_true", help='no latex output')
    parser.add_argument("lofile", type=str, nargs=1,
                        help="name of the low resolution convergence output file")
    parser.add_argument("hifile", type=str, nargs=1,
                        help="name of the high resolution convergence output file")
    parser.add_argument("veryhifile", type=str, nargs="?", default=None,
                        help="(optional) name of the very high resolution convergence output file")

    args = parser.parse_args()

    good_vars = {"density": r"$\rho$",
                 "xmom": r"$\rho u$",
                 "ymom": r"$\rho v$",
                 "rho_E": r"$\rho E$",
                 "rho_e": r"$\rho e$",
                 "Temp": r"$T$",
                 "rho_n": r"$\rho X(\mathrm{n})$",
                 "rho_p": r"$\rho X(\mathrm{p})$",
                 "rho_He4": r"$\rho X(\isotm{He}{4})$",
                 "rho_Cr48": r"$\rho X(\isotm{Cr}{48})$",
                 "rho_Fe52": r"$\rho X(\isotm{Fe}{52})$",
                 "rho_Fe54": r"$\rho X(\isotm{Fe}{54})$",
                 "rho_Ni56": r"$\rho X(\isotm{Ni}{56})$",
                 "rho_Fe56": r"$\rho X(\isotm{Fe}{56})$",
                 "rho_Ye": r"$\rho Y_e$",
                 "rho_abar": r"$\rho \bar{A}$",
                 "rho_bea": r"$\rho (B/A)$",
                 "rho_enuc": r"$\rho \dot{S}$"
                 }

    # sdc4
    file_lo = args.lofile[0]
    file_hi = args.hifile[0]
    file_vhi = args.veryhifile
    print(file_vhi)

    sdc4 = read_convergence(file_lo, file_hi, file_vhi)

    for v in sdc4.data:
        if v.name in good_vars.keys():
            if args.simple:
                name = v.name
            else:
                name = good_vars[v.name]
            print(v.get_table_line(pretty_name=name, simple=args.simple))

