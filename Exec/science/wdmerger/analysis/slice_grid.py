import sys

import yt
from yt.units import cm


def doit(plotfile):

    ds = yt.load(plotfile)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    xctr = 0.0 * xmin
    L_x = xmax - xmin

    yctr = 0.5 * (ymin + ymax)
    L_y = ymax - ymin

    width_frac = 1./3.

    field = "density"

    slc = yt.SlicePlot(ds, "z", field,
                       center=[xctr, yctr, 0.0*cm],
                       width=[width_frac*L_x, width_frac*L_y, 0.0*cm],
                       fontsize=14)

    slc.set_zlim(field, 1.e-4, 5.e7)

    slc.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s",
                      coord_system="axis", text_args={"color": "white"})

    slc.set_buff_size((3072, 3072))
    slc.set_axes_unit("cm")

    slc.annotate_grids()
    slc.save(f"wdmerger_slice_grid_{plotfile}.png")


if __name__ == "__main__":

    pf = sys.argv[-1]
    doit(pf)
