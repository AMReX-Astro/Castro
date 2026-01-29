#!/usr/bin/env python
# pylint: disable=too-many-branches, too-many-locals, too-many-statements
# pylint: disable=protected-access

import argparse
import concurrent.futures
import csv
import math
import re
import sys
from enum import Enum
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yt
from yt.units import cm
from yt.visualization._handlers import ColorbarHandler


class Transform(Enum):
    SLC = 0
    AVG = 1


def lookup_datasets():
    plotfile_paths = {}
    with open("ftime.out") as f:
        for line in f:
            pltfile, time_str = line.rstrip("\n").split()
            if Path(pltfile).exists():
                ms = round(float(time_str) * 1000)
                plotfile_paths[ms] = pltfile

    def get_dataset(ms: int) -> Any:
        assert ms in plotfile_paths, f"don't have plotfile for {ms}ms"
        print(f"loading {plotfile_paths[ms]} ({ms}ms)...")
        return yt.load(plotfile_paths[ms])

    return get_dataset


try:
    _get_dataset = lookup_datasets()
except FileNotFoundError:
    _get_dataset = None


def load_dataset(fname):
    path = Path(fname.rstrip("/"))
    if _get_dataset is not None and not path.exists():
        # try looking up by time
        if fname.endswith("ms"):
            return _get_dataset(int(fname.removesuffix("ms")))
    print(f"loading {path}...", flush=True)
    return yt.load(path)


def parse_args():
    ################################
    # set up parser and parse args #
    ################################

    description = "Script for tracking the position of a flame or shock front as a function of time."

    datasets_help = "Datasets to track the front position over."
    out_help = "Output filename for the tracking information."
    metrics_help = """A list of metrics to use for tracking. Should be fields followed by floating point
        numbers in the range (0, 1], indicating percent of maximum. For example, enuc 1 1e-3 Temp 1
        will track the position of max enuc and 1 / 1000th max enuc and the position of max
        temperature."""
    xlim_help = "Limits on the first dimension."
    ylim_help = "Limits on the second dimension."
    zlim_help = "Limits on the third dimension."
    res_help = "FRB resolution to use."
    transform_help = """Operation to apply to each extra dimension. Can be of format
        ax_ind:transform or just a sequence of transforms, with ax_ind assumed to be in descending
        order. Transforms: 0 => slice, 1 => average."""
    branch_help = """Whether to use the upper branch or lower branch when computing location. The upper
        branch is everything past the first instance of the local maximum, while the lower branch
        is everything before that. 0 => lower, 1 => upper. Default is upper."""
    tmax_help = "Maximum time to consider."
    jobs_help = "Use this many workers to process the plotfiles in parallel."
    show_plots_help = """Display plots of each field along with the corresponding metrics.
        If passed the output file from a previous run, the plots will use the
        global maxes from that file."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("datasets", nargs="+", help=datasets_help)
    parser.add_argument("-o", "--out", default="front_tracking.dat", help=out_help)
    parser.add_argument(
        "-m", "--metrics", nargs="*", default=["enuc", "1e-3"], help=metrics_help
    )
    parser.add_argument(
        "-x", "--xlim", nargs=2, type=float, metavar=("LOWER", "UPPER"), help=xlim_help
    )
    parser.add_argument(
        "-y", "--ylim", nargs=2, type=float, metavar=("LOWER", "UPPER"), help=ylim_help
    )
    parser.add_argument(
        "-z", "--zlim", nargs=2, type=float, metavar=("LOWER", "UPPER"), help=zlim_help
    )
    parser.add_argument("-r", "--res", nargs=2, type=int, help=res_help)
    parser.add_argument("-t", "--transform", nargs="+", help=transform_help)
    parser.add_argument("-b", "--branch", default=1, type=int, help=branch_help)
    parser.add_argument("--tmax", type=float, default=None, help=tmax_help)
    parser.add_argument("-j", "--jobs", type=int, default=1, help=jobs_help)
    parser.add_argument(
        "--show-plots", nargs="?", default=None, const=True, help=show_plots_help
    )

    args = parser.parse_args(sys.argv[1:])

    if args.show_plots and args.jobs > 1:
        parser.error("--show-plots is incompatible with --jobs > 1")

    # sort datasets by the numbers in the filename (natural sort)
    digit_pat = re.compile(r"\d+")
    args.datasets.sort(key=lambda x: tuple(map(int, digit_pat.findall(x))))

    # Assume same domain for all datasets
    ds = load_dataset(args.datasets[0])

    # used below for adjusting the default FRB resolution
    custom_limits = [False] * 3
    if args.xlim is None:
        args.xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else:
        args.xlim = tuple(val * cm for val in args.xlim)
        custom_limits[0] = True
    if args.ylim is None:
        args.ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else:
        args.ylim = tuple(val * cm for val in args.ylim)
        custom_limits[1] = True
    if args.zlim is None:
        args.zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
    else:
        args.zlim = tuple(val * cm for val in args.zlim)
        custom_limits[2] = True

    metrics = {}

    ls = []
    for item in args.metrics:
        try:
            ls.append(float(item))
        except ValueError:
            ls = metrics.setdefault(item, [])

    args.metrics = metrics

    if args.transform is None:
        args.transform = {Transform.SLC: [2], Transform.AVG: [1]}
    else:
        transform = {tr: [] for tr in Transform}

        for i, item in enumerate(args.transform):
            item = item.split(":")
            if len(item) == 1:
                (item,) = item
                transform[Transform(int(item))].append(args.spacedim - i - 1)
            elif len(item) == 2:
                ind, t = map(int, item)
                transform[Transform(t)].append(ind)
            else:
                parser.error("Invalid transform format.")

        args.transform = transform

    max_level = ds.index.max_level
    ref_ratio = int(np.prod(ds.ref_factors[0:max_level]))
    default_res = ds.domain_dimensions * ref_ratio
    # adjust the default resolution if using custom limits
    if any(custom_limits):
        limits = [args.xlim, args.ylim, args.zlim]
        for i, (is_custom, lim) in enumerate(zip(custom_limits, limits)):
            if not is_custom:
                continue
            new_width = lim[1] - lim[0]
            if new_width != ds.domain_width[i]:
                new_res = (default_res[i] / ds.domain_width[i] * new_width).value
                print(
                    "adjusting {}-axis resolution to match custom limits, from {} to {} ({})".format(
                        ds.coordinates.axis_order[i],
                        default_res[i],
                        math.ceil(new_res),
                        new_res,
                    )
                )
                default_res[i] = math.ceil(new_res)
    if args.res is None:
        args.res = default_res
    if len(args.res) < 3:
        args.res += [1] * (3 - len(args.res))
    for i, res in enumerate(args.res):
        if res <= 0:
            args.res[i] = default_res[i]

    # Eventually may want to generalize this to allow multiple axes
    # Then we would just return a point in 2D or 3D space

    transformed = sum(args.transform.values(), [])
    (args.axis,) = (ax for ax in range(3) if ax not in transformed)

    args.axis_name = ds.coordinates.axis_order[args.axis]

    if len(args.res) != 3:
        parser.error("wrong number of arguments to --res")

    del ds
    return args


#################
# Metrics class #
#################


def get_window_parameters(ds, axis, width=None, center="c"):
    """Some parameters controlling the frb window."""
    width = ds.coordinates.sanitize_width(axis, width, None)
    center, display_center = ds.coordinates.sanitize_center(center, axis)
    xax = ds.coordinates.x_axis[axis]
    yax = ds.coordinates.y_axis[axis]
    bounds = (
        display_center[xax] - width[0] / 2,
        display_center[xax] + width[0] / 2,
        display_center[yax] - width[1] / 2,
        display_center[yax] + width[1] / 2,
    )
    return bounds, center, display_center


def get_width(ds, xlim=None, ylim=None, zlim=None):
    """Get the width of the frb."""
    if xlim is None:
        xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else:
        xlim = xlim[0], xlim[1]

    if ylim is None:
        ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else:
        ylim = ylim[0], ylim[1]

    if zlim is None:
        zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
    else:
        zlim = zlim[0], zlim[1]

    xwidth = (xlim[1] - xlim[0]).in_cgs()
    ywidth = (ylim[1] - ylim[0]).in_cgs()
    zwidth = (zlim[1] - zlim[0]).in_cgs()

    return xwidth, ywidth, zwidth


def get_center(ds, xlim=None, ylim=None, zlim=None):
    """Get the coordinates of the center of the frb."""
    if xlim is None:
        xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else:
        xlim = xlim[0], xlim[1]

    if ylim is None:
        ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else:
        ylim = ylim[0], ylim[1]

    if zlim is None:
        zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
    else:
        zlim = zlim[0], zlim[1]

    xctr = 0.5 * (xlim[0] + xlim[1]).in_cgs()
    yctr = 0.5 * (ylim[0] + ylim[1]).in_cgs()
    zctr = 0.5 * (zlim[0] + zlim[1]).in_cgs()

    return xctr, yctr, zctr


# default plot settings for the debug slice plots, per field
default_limits = {
    "Temp": (1e6, 2e9),
    "enuc": (1e14, 1e19),
    "density": (1e-3, 1e8),
    "abar": (1, 6),
    "z_velocity": (-3e8, 3e8),
    "X(H1)": (0, 0.1),
    "X(He4)": (0, 0.9),
    "ash_density": (1e-2, 2e5),
}
default_cmaps = {
    "Temp": "magma_r",
    "abar": "plasma_r",
    "z_velocity": "bwr",
    "ash_density": "plasma_r",
}
default_logs = {
    "Temp": False,
    "abar": False,
    "z_velocity": False,
    "X(H1)": False,
    "X(He4)": False,
}


def get_avg_label(
    fi: yt.DerivedField, avg_axis_names: list[str], unit_registry: Any
) -> str:
    # slightly modified from yt.DerivedField.get_label()
    data_label = rf"$\left\langle\rm{{{fi.display_name}}}\right\rangle_{{\rm{{"
    data_label += r",".join(avg_axis_names)
    data_label += r"}}"
    units = yt.units.Unit(fi.units, registry=unit_registry)
    if not units.is_dimensionless:
        data_label += yt.visualization._commons._get_units_label(
            units.latex_representation()
        ).strip("$")
    data_label += r"$"
    return data_label


class Metrics:
    """Class for defining different measurements of the position of the flame front."""

    def __init__(self, ds, args):
        self.frbs = self.makefrbs(ds, args, globmax={})
        self.time = ds.current_time
        self._metrics = args.metrics
        self._upper_branch = args.branch > 0

    @classmethod
    def makefrbs(cls, ds: yt.data_objects.static_output.Dataset, args, globmax):
        # plotting code is integrated into this function so nothing goes out of sync
        width = get_width(ds, args.xlim, args.ylim, args.zlim)
        center = get_center(ds, args.xlim, args.ylim, args.zlim)

        region = [
            slice(*args.xlim, complex(0, args.res[0])),
            slice(*args.ylim, complex(0, args.res[1])),
            slice(*args.zlim, complex(0, args.res[2])),
        ]

        for axis in args.transform[Transform.SLC]:
            _, center, _ = get_window_parameters(ds, axis, width, center)
            region[axis] = center[axis]

        # The resolution in yt FixedResolutionBuffers is backwards
        # If we're doing a slice, we need to swap the steps
        is_slice = [isinstance(bounds, slice) for bounds in region]
        dim = sum(is_slice)

        if dim == 2:
            normal = is_slice.index(False)
            xax = ds.coordinates.x_axis[normal]
            yax = ds.coordinates.y_axis[normal]
            xslc, yslc = region[xax], region[yax]

            # Note: this is a workaround for
            # https://github.com/yt-project/yt/pull/5112 and
            # https://github.com/yt-project/yt/issues/3429.
            # Once those get fixed/merged, replace with this commented code:
            # region[xax] = slice(xslc.start, xslc.stop, yslc.step)
            # region[yax] = slice(yslc.start, yslc.stop, xslc.step)
            # frr = ds.r[tuple(region)]
            region[xax] = slice(xslc.start, xslc.stop)
            region[yax] = slice(yslc.start, yslc.stop)

            sl = ds.r._create_slice(region)
            frr = sl.to_frb(
                width=width[xax],
                height=width[yax],
                # The resolution in yt FixedResolutionBuffers is backwards
                resolution=(args.res[yax], args.res[xax]),
                center=center,
            )
            # end workaround
            del sl

        else:
            frr = ds.r[tuple(region)]

        # Create an FRB for each field and average over the remaining extra dimensions
        frbs = {"pos": frr[args.axis_name]}
        fields = args.metrics.keys()
        if args.show_plots:
            fig: mpl.figure.Figure
            ax: mpl.axes.Axes
            fig, plt_axes = plt.subplots(
                2,
                len(fields),
                layout="constrained",
                sharex="all",
                figsize=(7 * len(fields), 4.8),
            )
            if len(fields) == 1:
                plt_axes = [[ax] for ax in plt_axes]
            finfo = {field: ds._get_field_info(field) for field in fields}
            # share y axis between all the slice plots (this, along with
            # sharex="all" above, keeps them in sync when zooming and panning
            # with the QtAgg backend)
            for i in range(1, len(fields)):
                plt_axes[0, i].sharey(plt_axes[0, 0])

        for i, field in enumerate(fields):
            frbs[field] = frr[field]
            # check for https://github.com/yt-project/yt/pull/5112
            assert not np.any(np.isnan(frbs[field]))

            if args.show_plots:
                ax = plt_axes[0][i]
                kwargs = {}
                kwargs["origin"] = "lower"
                kwargs["aspect"] = "auto"
                kwargs["extent"] = (
                    xslc.start.value,
                    xslc.stop.value,
                    yslc.start.value,
                    yslc.stop.value,
                )
                # use this helper object to deal with colormaps
                cbh = ColorbarHandler()
                if field in default_limits:
                    kwargs["vmin"], kwargs["vmax"] = default_limits[field]
                if field in default_cmaps:
                    cbh.cmap = default_cmaps[field]
                if default_logs.get(field, finfo[field].take_log):
                    kwargs["norm"] = "log"
                # use the bottom value of the colormap for the background color
                cbh.background_color = None
                kwargs["cmap"] = cbh.cmap
                im = ax.imshow(frbs[field].value, **kwargs)
                ax.set_xlim(*kwargs["extent"][:2])
                ax.set_ylim(*kwargs["extent"][2:])
                ax.set_facecolor(cbh.background_color)
                ax.set_ylabel(
                    ds._get_field_info(ds.coordinates.axis_name[yax]).get_label()
                )
                fig.colorbar(im, ax=ax, label=finfo[field].get_label(), location="top")

        del frr

        axes = sorted([args.axis] + args.transform[Transform.AVG])

        avg_axis_names = []
        for i, ax in enumerate(reversed(axes)):
            if ax == args.axis:
                continue
            for field in frbs:
                frbs[field] = frbs[field].mean(axis=i)
            avg_axis_names.append(ds.coordinates.axis_name[ax])

        if args.show_plots:
            pos_frb = frbs["pos"]
            for i, field in enumerate(fields):
                ax = plt_axes[1][i]
                frb = frbs[field]
                ax.plot(pos_frb.value, frb.value)
                if default_logs.get(field, finfo[field].take_log):
                    ax.set_yscale("log")
                ax.set_xlim(xslc.start.value, xslc.stop.value)
                ax.set_xlabel(
                    ds._get_field_info(ds.coordinates.axis_name[args.axis]).get_label()
                )
                ax.set_ylabel(
                    get_avg_label(finfo[field], avg_axis_names, ds.unit_registry)
                )

                # plot the specific points for each threshold
                maxind = frb.argmax()
                # vertical line at the maximum
                ax.axvline(pos_frb[maxind].value, color="k", linestyle="dashed")
                # local max as circles, global max as squares
                for maxval, marker in [(frb.max(), "o"), (globmax.get(field), "s")]:
                    if maxval is None:
                        # running without global maxes available
                        continue
                    x_values = []
                    y_values = []
                    colors = []
                    for i, fac in enumerate(args.metrics[field]):
                        indices = cls._locate_helper(
                            frb, maxind, maxval, args.branch > 0, fac
                        )
                        if len(indices) > 0:
                            x_val = pos_frb[indices[0]].value
                            y_val = frb[indices[0]].value
                        else:
                            # no matching values found; place at the threshold
                            # value on top of the vertical max line
                            x_val = pos_frb[maxind].value
                            y_val = maxval * fac
                        x_values.append(x_val)
                        y_values.append(y_val)
                        colors.append(f"C{i}")
                    # draw these above the line plot
                    ax.scatter(
                        x_values,
                        y_values,
                        c=colors,
                        marker=marker,
                        alpha=0.8,
                        zorder=2.5,
                    )

            fig.suptitle(f"t={ds.current_time.to('s')}")
            plt.show()

        return frbs

    @staticmethod
    def _locate_helper(frb, maxind, maxval, upper_branch, fac):
        """Returns indices for all positions where `frb` drops to or reaches `fac * maxval`."""
        thresh = maxval * fac
        if upper_branch:
            # all indices to the right of the maximum where the field is less
            # than or equal to the threshold
            indices = np.nonzero(frb[maxind:] <= thresh)[0] + maxind
        else:
            # all indices to the left of the maximum where the field is greater
            # than or equal to the threshold
            indices = np.nonzero(frb[: maxind + 1] >= thresh)[0]
        return indices

    def locate(self, field, fac, globmax):
        """Returns position where `field` drops to or reaches `fac * max(field)`."""

        frb = self.frbs[field]
        maxind = frb.argmax()
        if globmax:
            maxval = globmax[field]
        else:
            maxval = frb[maxind]

        indices = self._locate_helper(frb, maxind, maxval, self._upper_branch, fac)
        if indices.size == 0:
            # no value found
            return np.nan * self.frbs["pos"].units
        # get the first matching index
        return self.frbs["pos"][indices[0]]

    @staticmethod
    def tostring(field, fac):
        # use :g here to remove the trailing .0
        return f"{field}[{fac*100:g}%]"

    @property
    def fields(self):
        return list(self._metrics.keys())

    def get_maxes(self):
        return {field: self.frbs[field].max() for field in self.fields}

    def getall(self, globmax):
        """Return positions for all metrics, in a dictionary keyed by metric string."""
        locs = {}

        for field, facs in self._metrics.items():
            for fac in facs:
                locs[self.tostring(field, fac)] = self.locate(field, fac, {})
                locs[self.tostring(field, fac) + "_gmax"] = self.locate(
                    field, fac, globmax
                )

        return locs


#########################################
# Compute positions and write data file #
#########################################


def process_dataset(fname: str, args: Any) -> tuple[float, Metrics] | tuple[None, None]:
    ds = load_dataset(fname)
    if args.tmax is not None and ds.current_time > args.tmax + 1e-6:
        return None, None
    return ds.current_time, Metrics(ds, args)


def show_plots(args):
    # make plots of the intermediate steps for the user to inspect
    globmax = {}
    if isinstance(args.show_plots, str):
        # parse out global maxes from a previous run
        with open(args.show_plots) as f:
            reader = csv.DictReader(f, delimiter=" ", quotechar='"')
            for k in reader.fieldnames:
                field = k.removeprefix("max_")
                if field == k:
                    # doesn't start with max_
                    continue
                globmax[field] = float("-inf")
            for row in reader:
                if args.tmax is not None and float(row["time"]) > args.tmax:
                    break
                for field in globmax:
                    globmax[field] = max(globmax[field], float(row[f"max_{field}"]))
    for fname in args.datasets:
        ds = load_dataset(fname)
        Metrics.makefrbs(ds, args, globmax)
        del ds


def main():
    args = parse_args()

    times = [None] * len(args.datasets)
    metrics = [None] * len(args.datasets)

    if args.show_plots is not None:
        show_plots(args)
        return

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_index = {
            executor.submit(process_dataset, fname, args): i
            for i, fname in enumerate(args.datasets)
        }
        try:
            for future in concurrent.futures.as_completed(future_to_index):
                i = future_to_index.pop(future)
                try:
                    time, metric = future.result()
                except Exception as exc:  # pylint: disable=broad-exception-caught
                    # note the exception and keep going
                    print(
                        f"{args.datasets[i]} generated an exception: {exc}",
                        file=sys.stderr,
                        flush=True,
                    )
                else:
                    times[i] = time
                    metrics[i] = metric
        except KeyboardInterrupt:
            print(
                "\n*** got ctrl-c, cancelling remaining tasks and waiting for existing ones to finish...\n",
                flush=True,
            )
            executor.shutdown(wait=True, cancel_futures=True)
            sys.exit(1)

    # remove entries for skipped plotfiles
    times = [t for t in times if t is not None]
    metrics = [m for m in metrics if m is not None]

    # compute the global max for each field
    all_maxes = [m.get_maxes() for m in metrics]
    fields = metrics[0].fields
    globmax = {}
    for field in fields:
        globmax[field] = max(*(maxes[field] for maxes in all_maxes))

    # We need the global max to have been computed already to properly track the contour
    loclist = [m.getall(globmax) for m in metrics]
    cols = loclist[0].keys()

    with open(args.out, "w", newline="") as file:
        csvwriter = csv.writer(
            file, delimiter=" ", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        header = ["time"]
        header.extend(f"max_{field}" for field in fields)
        header.extend(cols)
        csvwriter.writerow(header)

        for time, maxes, locs in zip(times, all_maxes, loclist):
            row = [time.to_value("s")]
            row.extend(maxes[field].value for field in fields)
            row.extend(locs[col].to_value("cm") for col in cols)
            csvwriter.writerow(row)

    print("Task completed.")


if __name__ == "__main__":
    main()
