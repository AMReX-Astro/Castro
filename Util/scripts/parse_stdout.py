#!/usr/bin/env python

import re
import sys

RE_STEP = r"^STEP = (\d+)"
RE_LEVEL = r"^\[Level (\d+) step (\d+)\]"
RE_SDC_ITER = r"^Beginning SDC iteration (\d+) of (\d+)"
RE_VODE = r"DVODE:"

ofile = sys.argv[-1]

with open(ofile, "r") as of:

    step = 0
    level = 0
    sdc_iter = 0

    while (line := of.readline()):

        if g := re.search(RE_STEP, line):
            step = g.groups()[0]
            continue

        if g := re.search(RE_LEVEL, line):
            level = g.groups()[0]
            continue

        if g := re.search(RE_SDC_ITER, line):
            sdc_iter = g.groups()[0]
            continue

        if g := re.search(RE_VODE, line):
            print(f"step: {step}, level: {level}, sdc iter: {sdc_iter} ", line)

