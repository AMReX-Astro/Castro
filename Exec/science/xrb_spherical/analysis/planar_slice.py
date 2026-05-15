#!/usr/bin/env python

import sys
import os
import argparse
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.pyplot as plt
from matplotlib.colors import (LogNorm, Normalize, SymLogNorm)
import numpy as np
import yt
from yt.frontends.boxlib.api import CastroDataset
from matplotlib.patches import Rectangle
import pynucastro as pyna
from yt.units import (cm, s, g)
from front_tracker import track_front

## Diverging colormap to use: "RdBu_r" or "coolwarm"
## Perceptually Uniform Sequential colormap to use: 'viridis', 'plasma', 'inferno', 'magma'
## Cyclic (rainbow) colormap (NOT IDEAL TO USE): 'jet', 'turbo' (better than jet for smoothing)

# Set some fontsize settings
SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)


FIELD_PRESETS = {
    "Temp": dict(cmap="magma_r" , norm=LogNorm(vmin=5e7, vmax=2.5e9), label="Temperature [K]"),
    "density": dict(cmap="viridis", norm=LogNorm(vmin=1e-1, vmax=1e8), label=r"Density [g cm$^{-3}$]"),
    "pressure": dict(cmap="inferno_r", norm=LogNorm(vmin=1e16, vmax=1e25), label=r"Pressure [dyn cm$^{-2}$]"),
    "enuc": dict(cmap="inferno_r", norm=LogNorm(vmin=1e16, vmax=5e19, clip=True),
                 label=r"$\dot{e}_{\mathrm{nuc}}$ [erg g${}^{-1}$ s${}^{-1}$]"),
    "eint_e": dict(cmap="inferno_r", norm=Normalize(vmin=5e16, vmax=1e18, clip=True),
                 label=r"$e$ [erg g${}^{-1}$]"),
    "diff_coeff": dict(cmap="plasma_r", norm=LogNorm(vmin=5e5, vmax=3e6, clip=True),
                 label=r"$D$ [cm${}^{2}$ s${}^{-1}$]"),
    "abar": dict(cmap="plasma_r", norm=Normalize(vmin=4, vmax=8), label=r"$\bar{A}$"),
    "ash": dict(cmap="plasma_r", norm=LogNorm(vmin=1e-2, vmax=1e6), label=r"$\rho X \left(ash\right) [g cm$^{-3}$]$"),
    "grad_rho": dict(cmap="plasma_r", norm=LogNorm(vmin=1e-2, vmax=1e6), label=r"$|\nabla \rho|$"),
    "MachNumber": dict(cmap="inferno_r", norm=LogNorm(vmin=1e-4, vmax=1.1), label="MachNumber"),
    "divu": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-1e5, vmax=1e5, linthresh=1e2), label=r"$\nabla \cdot U$"),
    "x_velocity": dict(cmap="RdBu_r", norm=Normalize(vmin=-1e6, vmax=1e6), label=r"v$_r$ [cm s$^{-1}$]"),
    "y_velocity": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-1e7, vmax=1e7, linthresh=1e4), label=r"v$_\theta$ [cm s$^{-1}$]"),
    "z_velocity": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-1e9, vmax=1e9, linthresh=1e7), label=r"v$_\phi$ [cm s$^{-1}$]"),
    "Dvr_Dt": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-5e13, vmax=5e13, linthresh=5e11), label=r"a$_r$ [cm s$^{-2}$]"),
    "Dvt_Dt": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-1e13, vmax=1e13, linthresh=1e11), label=r"a$_\theta$ [cm s$^{-2}$]"),
    "Dvp_Dt": dict(cmap="RdBu_r", norm=SymLogNorm(vmin=-1e11, vmax=1e11, linthresh=1e9), label=r"a$_\phi$ [cm s$^{-2}$]")
}

# Default Preset for plotting different massfractions
DEFAULT_PRESET = dict(cmap="turbo", norm=LogNorm(vmin=1e-6, vmax=1e-2), label="")

CONTOUR_LEVEL_PRESETS = {
    "pressure": np.logspace(np.log10(1e17), np.log10(1e24), 10),
    "density": np.logspace(np.log10(1), np.log10(1e7), 10),
    "enuc": np.logspace(np.log10(1e16), np.log10(5e19), 6),
    "abar": np.linspace(4, 8, 10),
}

def _ash(field, data):
    '''
    Computes the rho X_ash.
    Here ash is anything heavier than Oxygen, but exclude Fe and Ni
    '''

    ds = data.ds
    field_list = ds.field_list

    # If X(ash) is directly stored, then just use that
    if ("boxlib", "X(ash)") in field_list:
        return data["boxlib", "X(ash)"] * data["boxlib", "density"]

    # If we cannot find X(ash) as a available field then compute manually.
    # First check if the plotfile has massfractions -- i.e. using plt not smallplt
    has_species = any(f[1].startswith("X(") for f in field_list)
    if not has_species:
        raise RuntimeError(
            "Derived field 'ash' requires species mass fractions X(nuc), "
            "but no X(...) fields were found in this plotfile."
        )

    rho = data["boxlib", "density"]
    rhoAsh = rho * 0.0
    for f in field_list:
        # If the first two letters are "X(", then we're dealing with species massfractions
        if f[1][:2] == "X(":
            # Then extract out what species we have, assuming the format is "X(...)"
            speciesName = f[1][2:-1]
            nuc = pyna.Nucleus(speciesName)

            # Include elements beyond oxygen but don't include Ni56
            if nuc.Z > 8.0 and nuc.Z != 28:
                rhoAsh += rho * data[f]

    return rhoAsh

def _grad_rho(field, data):
    """
    absolution value Gradient of density -- Contour line to show shock structure.
    """
    rho   = data["boxlib", "density"]
    r     = data["index", "r"]
    theta = data["index", "theta"]

    gradrho_r = np.gradient(rho.to_ndarray(), r[:,0,0].to_ndarray(), axis=0)
    gradrho_theta = np.gradient(rho.to_ndarray(), theta[0,:,0].to_ndarray(), axis=1) / r.to_ndarray()
    gradrho = np.sqrt(gradrho_r*gradrho_r + gradrho_theta*gradrho_theta) * g/cm**3

    return gradrho

def _Dvr_Dt(field, data):
    """
    Material derivative of velocity along radial direction
    Dvr/Dt = -1/rho dp/dr + g + 1/r (vt^2 + vp^2) + 2 Omega vp sin(theta)
    """
    omega = 2.0 * np.pi / data.ds.parameters.get("castro.rotational_period") / s
    grav  = data.ds.parameters.get("gravity.const_grav") * cm / s**2
    P     = data["boxlib", "pressure"]
    vr    = data["boxlib", "x_velocity"]
    vt    = data["boxlib", "y_velocity"]
    vp    = data["boxlib", "z_velocity"]
    rho   = data["boxlib", "density"]
    r     = data["index", "r"]
    theta = data["index", "theta"]

    F_P = -np.gradient(P.to_ndarray(), r[:,0,0].to_ndarray(), axis=0) / rho.to_ndarray() *cm/s**2
    F_geom = (vt**2 + vp**2) / r
    F_coriolis = 2 * omega * vp * np.sin(theta)

    return F_P + grav + F_geom + F_coriolis

def _Dvt_Dt(field, data):
    """
    Material derivative of velocity along theta direction
    Dvt/Dt = -1/rho dp/rdtheta + 1/r (-vr vt + vp^2 cot(theta)) + 2 Omega vp cos(theta)
    """
    omega = 2.0 * np.pi / data.ds.parameters.get("castro.rotational_period") / s
    P     = data["boxlib", "pressure"]
    vr    = data["boxlib", "x_velocity"]
    vt    = data["boxlib", "y_velocity"]
    vp    = data["boxlib", "z_velocity"]
    rho   = data["boxlib", "density"]
    r     = data["index", "r"]
    theta = data["index", "theta"]

    gradP_theta = np.gradient(P, theta[0,:,0], axis=1)
    F_P = -gradP_theta / (r*rho)
    F_geom = (vp**2 / np.tan(theta) - vr * vt) / r
    F_coriolis = 2 * omega * vp * np.cos(theta)

    return F_P + F_geom + F_coriolis

def _Dvp_Dt(field, data):
    """
    Material derivative of velocity along phi direction
    Dvt/Dt = -1/rho dp/rsin(theta)dphi - 1/r (vr vp + vt vp cot(theta)) - 2 Omega (vt cos(theta) + vr sin(theta))
    """
    omega = 2.0 * np.pi / data.ds.parameters.get("castro.rotational_period") / s
    vr    = data["boxlib", "x_velocity"]
    vt    = data["boxlib", "y_velocity"]
    vp    = data["boxlib", "z_velocity"]
    rho   = data["boxlib", "density"]
    r     = data["index", "r"]
    theta = data["index", "theta"]

    F_P = 0  # Axisymmetric -- no pressure grad in phi
    F_geom = -(vt * vp / np.tan(theta) + vr * vp) / r
    F_coriolis = -2 * omega * (vt * np.cos(theta) + vr * np.sin(theta))

    return F_P + F_geom + F_coriolis

def get_level_bounds(ds, level):
    """Get the bounding box of all grids at a given level."""
    grids = [g for g in ds.index.grids if g.Level == level]
    if not grids:
        return None
    r_left   = min(g.LeftEdge[0].value  for g in grids)
    r_right  = max(g.RightEdge[0].value for g in grids)
    th_left  = min(g.LeftEdge[1].value  for g in grids)
    th_right = max(g.RightEdge[1].value for g in grids)
    return r_left, r_right, th_left, th_right


def plot_finest_available(ax, ds, field):
    """Plot finer level data by overplotting each level."""

    # Choose different settings based on what field parameter is plotted
    preset = FIELD_PRESETS.get(field, DEFAULT_PRESET)

    max_level = ds.max_level
    pcm = None

    for level in range(1, max_level + 1):
        # check if this level exists
        bounds = get_level_bounds(ds, level)
        if bounds is None:
            continue

        # restrict covering grid to just the refined region
        r_left, r_right, th_left, th_right = bounds
        left_edge  = ds.arr([r_left,  th_left,  ds.domain_left_edge[2].value],  "code_length")
        #right_edge = ds.arr([r_right, th_right, ds.domain_right_edge[2].value], "code_length")

        dims = ds.domain_dimensions * ds.refine_by**level
        dims[2] = 1
        cg = ds.smoothed_covering_grid(level=level, left_edge=left_edge, dims=dims)

        var   = cg["boxlib", field][:, :, 0].to_ndarray()
        r     = cg["index", "r"][:, :, 0].to("km").to_ndarray()
        theta = cg["index", "theta"][:, :, 0].to_ndarray()

        pcm = ax.pcolormesh(theta, r, var, cmap=preset["cmap"], norm=preset["norm"], shading="auto")

    return pcm

def extract_info(ds, data_level=0):
    '''
    Extract relevant info from dataset.
    Mainly returns vertical height info and data grid

    Parameters
    ----------
    ds: CastroDataset
    level: what AMR level dataset to get. Going higher level
        doesn't seem to give sharper image somehow.
    '''
    dims = ds.domain_dimensions * ds.refine_by**data_level

    # change phi-direction to have single level
    dims[2] = 1

    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    rl = left_edge[0] * 1e-5
    rr = right_edge[0] * 1e-5
    r = [rl, rr] # in km
    theta = [left_edge[1], right_edge[1]]

    # Get grid data
    cg = ds.smoothed_covering_grid(level=data_level, left_edge=left_edge, dims=dims)
    return (r, theta, cg)

def add_derived_fields(ds):
    """Given a dataset add the derived fields"""

    has_species = any(f[1].startswith("X(") for f in ds.field_list)
    if has_species:
        ds.add_field(("boxlib", "ash"), function=_ash,
                     display_name=r"\rho X\left(ash\right)",
                     units="auto", sampling_type="cell")

    if ("boxlib", "pressure") in ds.field_list:
        ds.add_field(("boxlib", "Dvr_Dt"),
                     function=_Dvr_Dt,
                     units="cm/s**2", sampling_type="cell",
                     display_name=r"$D v_r / Dt$")

        ds.add_field(("boxlib", "Dvt_Dt"),
                     function=_Dvt_Dt,
                     units="cm/s**2", sampling_type="cell",
                     display_name=r"$D v_\theta / Dt$")

        ds.add_field(("boxlib", "Dvp_Dt"),
                     function=_Dvp_Dt,
                     units="cm/s**2", sampling_type="cell",
                     display_name=r"$D v_\phi / Dt$")

    if ("boxlib", "density") in ds.field_list:
        ds.add_field(("boxlib", "grad_rho"),
                     function=_grad_rho,
                     units="g/cm**3", sampling_type="cell",
                     display_name=r"$|\nabla \rho|$")

def planar_slice(fnames:list[str], fields:list[str],
                 figsize=(16, 9),
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 contour_field=None,
                 overplot_fine_levels=False,
                 annotate_front=False,
                 annotate_velocity_streamlines=False,
                 annotate_acceleration_streamlines=False,
                 annotate_grids=False,
                 outName=None):

    # Get actual dataset from fnames
    ts = [CastroDataset(fname) for fname in fnames]

    # Allow only multiple fnames or multiple fields.
    if len(fields) > 1 and len(fnames) > 1:
        raise ValueError("Only multiple fields or multiple fnames are allowed")

    # If multiple data files passed in, assume time-series plot
    # assert that we're plotting a single field
    nrow = len(fnames)
    if len(fields) > 1:
        # If multiple fields, then plot multi-field plot
        nrow = len(fields)
        assert len(fnames) == 1

    fig = plt.figure(figsize=figsize)
    grid = ImageGrid(fig, 111, nrows_ncols=(nrow, 1),
                     axes_pad=0.5, label_mode="L", cbar_location="right",
                     cbar_mode="each", cbar_size="2.5%", cbar_pad="0%",
                     aspect=False) # do asepect=False to allow elongation in vertical dir

    # Loop over all possible dataset
    for j, ds in enumerate(ts):

        # Add derived fields
        add_derived_fields(ds)

        # Needed for smoothed_covering_grid to use fine level data
        ds.force_periodicity()

        # Find current time
        time = ds.current_time.in_units("ms")
        if time.value < 1e-1:
            time = ds.current_time.in_units("us")

        # Get level 0 data
        r, theta, cg0 = extract_info(ds, data_level=0)

        if ymin is None:
            ymin = r[0].value
        if ymax is None:
            ymax = 0.5*(r[1].value - r[0].value)
        if xmin is None:
            xmin = theta[0].value
        if xmax is None:
            xmax = theta[1].value

        # Compute necessary data for plotting front beforehand
        if annotate_front:
            front_tracking_data = track_front(ds)

        # Loop over all possible fields
        for i, field in enumerate(fields):

            # Get ax
            ax  = grid[i+j*len(fields)]
            cax = grid.cbar_axes[i+j*len(fields)]

            # Get plotting variables
            var    = cg0["boxlib", field][:, :, 0].to_ndarray()
            r_arr = cg0["index", "r"][:, :, 0].to("km").to_ndarray() # longitudinal distance in km
            theta_arr  = cg0["index", "theta"][:, :, 0].to_ndarray()

            # Choose different settings based on what field parameter is plotted
            preset = FIELD_PRESETS.get(field, DEFAULT_PRESET)

            # Plot level 0 data using pcolormesh
            pcm = ax.pcolormesh(theta_arr, r_arr, var, shading="auto",
                                cmap=preset["cmap"], norm=preset["norm"])

            # Overplot finer level data if available
            if ds.max_level > 0 and overplot_fine_levels:
                pcm = plot_finest_available(ax, ds, field)

            # Annotate optional contour lines based on contour_field -- this uses level 0 data for now
            if contour_field is not None:
                contour_var = cg0["boxlib", contour_field][:, :, 0].to_ndarray()
                contour_levels = CONTOUR_LEVEL_PRESETS.get(contour_field, 10)
                ax.contour(theta_arr, r_arr, contour_var, levels=contour_levels, linewidths=1, colors="k")

            # Annotate optional velocity or acceleration streamlines
            # note here x-velocity is in radial so its technically vertical
            # this uses level 0 data for now
            if annotate_velocity_streamlines:
                vt = cg0["boxlib", "y_velocity"][:, :, 0].to_ndarray()
                vr = cg0["boxlib", "x_velocity"][:, :, 0].to_ndarray()
                ax.streamplot(theta_arr, r_arr, vt, vr, color="tab:green", linewidth=2, density=1)

            if annotate_acceleration_streamlines and ("boxlib", "pressure") in ds.field_list:
                at = cg0["boxlib", "Dvt_Dt"][:, :, 0].to_ndarray()
                ar = cg0["boxlib", "Dvr_Dt"][:, :, 0].to_ndarray()
                ax.streamplot(theta_arr, r_arr, at, ar, color="tab:orange", linewidth=2, density=1)

            # Plot vertical line to indicate flame front and ash front
            if annotate_front:
                ash_front_theta = front_tracking_data["ash_theta"]
                flame_front_theta = front_tracking_data["flame_theta"]
                flame_tail_theta = front_tracking_data["flame_tail"]

                text_pos = ymin + 0.8*(ymax - ymin)
                ax.axvline(ash_front_theta, linestyle="-.", color="k", linewidth=1.5)
                ax.text(ash_front_theta, text_pos, "Ash\nFront",
                        color="k", fontsize=12, ha="center", va="center", rotation=90)

                ax.axvline(flame_front_theta, linestyle="-.", color="k", linewidth=1.5)
                ax.text(flame_front_theta, text_pos, "Flame\nFront",
                        color="k", fontsize=12, ha="center", va="center", rotation=90)

                ax.axvline(flame_tail_theta, linestyle="-.", color="k", linewidth=1.5)
                ax.text(flame_front_theta, text_pos, "Flame\nTail",
                        color="k", fontsize=12, ha="center", va="center", rotation=90)

            # Annotate optional AMR grids
            if annotate_grids:
                max_level = ds.max_level
                level_colors = plt.get_cmap("Set1", max_level + 1)
                for g in ds.index.grids:
                    r_lo  = g.LeftEdge[0].value  * 1e-5 # in km
                    r_hi  = g.RightEdge[0].value * 1e-5 # in km

                    th_lo = g.LeftEdge[1].value
                    th_hi = g.RightEdge[1].value

                    height = r_hi - r_lo
                    width  = th_hi - th_lo

                    rect = Rectangle((th_lo, r_lo), width, height, linewidth=0.6, edgecolor=level_colors(g.Level),
                                     facecolor="none", alpha=0.7, label=f"Level {g.Level}",
                                     transform=ax.transData,
                                     clip_on=True)
                    ax.add_patch(rect)

                # Deduplicated legend entries
                handles = [Rectangle((0,0), 1, 1, edgecolor=level_colors(l), facecolor="none", label=f"Level {l}")
                           for l in range(max_level + 1)]
                ax.legend(handles=handles, loc="upper right", framealpha=0.5)

            # Can I somehow annotate the slope line between unburnt fuel and ash interface -- and find its angle

            # Some formatting
            if i == 0:
                ax.set_title(f"t = {time:.2f}")

            ax.set_xlabel(r"$\theta$ [rad]")
            ax.set_ylabel("Radial [km]")
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

            cb_label = preset["label"]
            if not cb_label:
                cb_label = field
            fig.colorbar(pcm, cax=cax, label=cb_label)

    fig.set_size_inches(*figsize)
    fig.tight_layout()

    # Store to output, otherwise show plot
    if outName is not None:
        fig.savefig(outName, format="png", bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        A planar slice plot script for xrb_spherical problem.
        Given a list of plotfiles or a list of field parameters,
        it plots multiple slice plots.
        """)

    parser.add_argument('fnames', nargs='+', type=str,
                        help="""dataset file names for plotting. Accepts one or more datasets.
                        If multiple file names are given, a grid of slice plots of different
                        files will be plotted for a given field parameter.
                        Note that either fnames or field must be single valued.""")
    parser.add_argument('-f', '--fields', nargs='+', type=str,
                        help="""field parameters for plotting. Accepts one or more datasets.
                        If multiple parameters are given, a grid of slice plots of different
                        field parameters will be plotted for a given fname.
                        Note that either fnames or fields must be single valued.
                        """)
    parser.add_argument("--figsize", nargs=2, type=float, default=[16, 9],
                        metavar=("WIDTH", "HEIGHT"), help="Figure size in inches.")
    parser.add_argument("--xmin", type=float, default=None, metavar="THETA",
                        help="Minimum theta for plot xlim")
    parser.add_argument("--xmax", type=float, default=None, metavar="THETA",
                        help="Maximum theta for plot xlim")
    parser.add_argument("--ymin", type=float, default=None, metavar="R",
                        help="Minimum r [km] for plot ylim")
    parser.add_argument("--ymax", type=float, default=None, metavar="R",
                        help="Maximum r [km] for plot ylim")
    parser.add_argument("--contour-field", default=None,
                        help="Field variable to use for overplotting contour lines (e.g. 'pressure').")
    parser.add_argument("--overplot-fine-levels", action="store_true",
                        help="""Overplot finer AMR level data on top of level 0.
                        This can make plotting a lot slower""")
    parser.add_argument("--annotate-front", action="store_true",
                        help="Overlay vertical lines to indicate flame front and ash front.")
    parser.add_argument("--annotate-velocity-streamlines", action="store_true",
                        help="Overlay velocity streamlines.")
    parser.add_argument("--annotate-acceleration-streamlines", action="store_true",
                        help="Overlay acceleration streamlines.")
    parser.add_argument("--annotate-grids", action="store_true",
                        help="Overlay AMR grid boundaries colored by level.")
    parser.add_argument("-o", "--output", default=None, type=str, metavar="FILENAME",
                        help="Output filename (PNG). If not set, shows interactive plot.")

    args = parser.parse_args()

    if len(args.fnames) > 1 and len(args.fields) > 1:
        parser.error("Either fnames or fields must be single valued!")

    planar_slice(
        fnames                            = args.fnames,
        fields                            = args.fields,
        figsize                           = tuple(args.figsize),
        xmin                              = args.xmin,
        xmax                              = args.xmax,
        ymin                              = args.ymin,
        ymax                              = args.ymax,
        contour_field                     = args.contour_field,
        overplot_fine_levels              = args.overplot_fine_levels,
        annotate_front                    = args.annotate_front,
        annotate_velocity_streamlines     = args.annotate_velocity_streamlines,
        annotate_acceleration_streamlines = args.annotate_acceleration_streamlines,
        annotate_grids                    = args.annotate_grids,
        outName                           = args.output,
    )
