# loop over all of the plotfiles in the current directory
# and render slices, output into a subdirectory {field}_slices

import argparse
from pathlib import Path

import yt


def doit(field, basename, is_small_plotfiles):

    curr_dir = Path.cwd()

    # find all of the plotfiles in the current directory
    plt = "smallplt" if is_small_plotfiles else "plt"
    plotfiles = [f for f in Path(".").glob(f"{basename}_{plt}*[0-9]") if f.is_dir()]

    out_dir = curr_dir / f"{field}_slices"
    if not out_dir.is_dir():
        out_dir.mkdir()

    # find any plotfiles we've already processed
    processed = [f for f in Path(".").glob(f"{basename}_{plt}*[0-9]") if Path(f"{out_dir}/{f}_{field}.png").is_file()]
    remaining = [f for f in plotfiles if f not in processed]

    ts_04 = yt.DatasetSeries(remaining)
    text_color = "white"


    for pfile in ts_04:
        ds = yt.load(str(pfile), hint="castro")

        sp = yt.SlicePlot(ds, 'z', field,
                          buff_size=(7680, 3840),
                          origin="domain", fontsize="14")

        if field == "magvel":
            sp.set_cmap(field, "inferno")
            sp.set_zlim(field, 1.0e5, 1.0e8)

        elif field == "magvort":
            sp.set_cmap(field, "viridis")
            sp.set_zlim(field, 1.e-4, 100)

        elif field == "enuc":
            sp.set_log(field, True, linthresh=1.e9)
            sp.set_zlim(field, 0., 1.e14)
            sp.set_cmap(field, "plasma")

        elif field == "Temp":
            sp.set_zlim(field, 1.e6, 1.e8)
            sp.set_cmap(field, "plasma")

        sp.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s",
                         coord_system="axis",
                         text_args={"color": text_color, "fontsize": "14"})
        sp.set_axes_unit("km")

        outfile = out_dir / f"{pfile}_{field}.png"
        sp.save(str(outfile))
        ds.index.clear_all_data()

if __name__ == "__main__":

    p = argparse.ArgumentParser(description="make slice plots of the suite of nova plotfiles")

    p.add_argument("--field", type=str, default="magvort",
                   help="field to visualize")
    p.add_argument("--basename", type=str, default="nova",
                   help="basename of the plotfiles")
    p.add_argument("--plt", action="store_true",
                   help="plotfiles are plt instead of smallplt?")
    args = p.parse_args()

    is_small_plotfiles = not args.plt
    field = args.field
    basename = args.basename

    doit(field, basename, is_small_plotfiles)
