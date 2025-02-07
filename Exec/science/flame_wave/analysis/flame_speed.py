#!/usr/bin/env python3

import argparse
import itertools
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.stats import linregress

plt.rcParams.update(
    {
        "font.family": "serif",
    }
)


def measure_and_plot(dat, args):  # pylint: disable=too-many-locals, too-many-statements
    """
    Plots front_location vs. time on the current pyplot figure, as well as any
    reasonable linear fits. *tmin* and *tmax* should give a range of times for
    which the flame front has stabilized, and the slope obtained from linear
    regression will yield an accurate measurement of the front speed.
    """

    slopes = {}

    tarr = dat["time"]
    # sort by time
    indx = tarr.argsort()
    tarr = tarr[indx]

    # pylint: disable-next=use-dict-literal
    common_props = dict(linewidth=2, alpha=0.8)

    plt.figure(figsize=(8, 6))

    # give the same color to <column> and <column>_gmax, regardless of where
    # they occur in columns
    color_iter = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    color_lookup = {}
    for column in args.columns:
        m = re.match(r"(?P<field>.+?)\[(?P<percent>.+?)%\](?P<gmax>_gmax)?", column)
        assert m is not None
        field = m["field"]
        percent = float(m["percent"])
        global_max = m["gmax"] is not None

        key = (field, percent)
        if key not in color_lookup:
            color_lookup[key] = next(color_iter)
        kwargs = {"color": color_lookup[key]}
        if not global_max and len(args.columns) > 1:
            kwargs["linestyle"] = "dashed"

        radarr = dat[column][indx]

        plt.plot(tarr * 1000, radarr, **common_props, **kwargs)

        reg_mask = ~np.isnan(radarr)
        if args.tmin is not None:
            reg_mask &= tarr >= args.tmin / 1000
        if args.tmax is not None:
            reg_mask &= tarr <= args.tmax / 1000

        m, b, r, _, sd = linregress(tarr[reg_mask], radarr[reg_mask])
        print(f"{column:<30}\t{m=:.2e}\t{r**2=:.3f}\t{sd=:.2e}")

        if r**2 > 0.8:
            if len(args.columns) > 1:
                color = kwargs["color"]
            else:
                color = "black"
            plt.plot(
                tarr[reg_mask] * 1000,
                m * tarr[reg_mask] + b,
                **common_props,
                color=color,
                linestyle="dotted",
            )
            slopes[column] = m

    ax = plt.gca()
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.major.formatter.set_useMathText(True)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # build a custom legend
    legend_elements = []
    if len(args.columns) > 1:
        legend_elements.extend(
            [
                Line2D([0], [0], **common_props, c="k", ls="dashed", label="local max"),
                Line2D([0], [0], **common_props, c="k", ls="solid", label="global max"),
            ]
        )
    legend_elements.extend(
        Patch(facecolor=color, label=f"{field} {percent:g}%")
        for (field, percent), color in color_lookup.items()
    )
    plt.legend(handles=legend_elements)

    plt.xlabel("t [ms]")
    plt.ylabel("r [cm]")
    plt.title("Front Position vs. Time", fontsize=24)

    for item in (ax.xaxis.label, ax.yaxis.label):
        item.set_fontsize(20)
        item.set_color("dimgrey")

    for item in (ax.get_xticklabels() + ax.get_yticklabels()) + [ax.yaxis.offsetText]:
        item.set_fontsize(16)
        item.set_color("dimgrey")

    plt.tight_layout()

    return slopes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", type=argparse.FileType("r"))
    parser.add_argument(
        "--tmin",
        type=float,
        help="ignore everything before this time for the linear regression",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        help="ignore everything after this time for the linear regression",
    )
    parser.add_argument("columns", nargs="+")
    args = parser.parse_args()

    dat = pd.read_csv(args.data_file, sep=r"\s+")

    plt.style.use("ggplot")
    measure_and_plot(dat, args)
    plt.show()


if __name__ == "__main__":
    main()
