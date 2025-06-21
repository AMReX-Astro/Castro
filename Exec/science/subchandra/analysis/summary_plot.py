import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.ticker as ptick

import yt
from yt.units import cm
from yt.frontends.boxlib.api import CastroDataset

def doit(fns, f, prefix):

    fig = plt.figure()

    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    # These choices of keyword arguments produce a four panel plot with a single
    # shared narrow colorbar on the right hand side of the multipanel plot. Axes
    # labels are drawn for all plots since we're slicing along different directions
    # for each plot.
    grid = AxesGrid(
        fig,
        111,
        nrows_ncols=(2, 2),
        axes_pad=0.2,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_pad="3%"
    )

    xmin = 0*cm
    xmax = 4.e8*cm
    ymin = 1.0e9*cm
    ymax = 1.4e9*cm

    xctr = 0.5*(xmin + xmax)
    yctr = 0.5*(ymin + ymax)
    L_x = xmax - xmin
    L_y = ymax - ymin


    for i, fn in enumerate(fns):

        if fn is None:
            grid[i].remove()
            continue

        # Load the data and create a single plot
        ds = CastroDataset(fn)  # load data
        sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0*cm],
                          width=[L_x, L_y, 0.0*cm], fontsize="11")

        # Ensure the colorbar limits match for all plots
        if f == "Temp":
            sp.set_zlim(f, 5.e7, 4e9)
            sp.set_cmap(f, "magma")
        elif f == "enuc":
            sp.set_log(f, True, linthresh=1.e18)
            sp.set_zlim(f, -1.e22, 1.e22)
            sp.set_cmap(f, "bwr")
        elif f == "density":
            sp.set_zlim(f, 1.e-3, 5.e8)
        elif f == "z_velocity":
            sp.set_zlim(f, -2.e8, 2.e8)
            sp.set_log(f, False)
            sp.set_cmap(f, "bwr")
        elif f == "abar":
            sp.set_zlim(f, 4, 28)
            sp.set_log(f, False)
            sp.set_cmap(f, "plasma_r")

        sp.annotate_text((0.05, 0.05), fn.split("/")[0].replace("_tols", ""),
                         coord_system="axis", text_args={"color": "black"})

        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = sp.plots[f]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]


        # Finally, this actually redraws the plot.
        sp._setup_plots()

        grid[i].axes.ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True)
        grid[i].axes.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)


    fig.set_size_inches(10, 9)
    fig.tight_layout()
    plt.savefig(f"{prefix}_summary_plot.pdf")

def main():

    sdc_fns = [
        "sdc_subch2_cfl0.1_tols/subch_plt08239",
        "sdc_subch2_cfl0.2_tols/subch_plt04163",
        "sdc_subch2_cfl0.4_tols/subch_plt02135",
        "sdc_subch2_cfl0.2_tol1e-8/subch_plt04161"
    ]


    strang_fns = [
        "strang_subch2_cfl0.1_tols/subch_plt08240",
        "strang_subch2_cfl0.2_tols/subch_plt04182",
        "strang_subch2_cfl0.4_tols/subch_plt02172",
        "strang_subch2_cfl0.2_tol1.e-8/subch_plt04161"
    ]

    doit(sdc_fns, "enuc", "sdc_enuc")
    doit(sdc_fns, "abar", "sdc_abar")

    doit(strang_fns, "enuc", "strang_enuc")
    doit(strang_fns, "abar", "strang_abar")



if __name__ == "__main__":
    main()
