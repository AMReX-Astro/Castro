import yt

ds = yt.load("plt00000")

slc = yt.SlicePlot(ds, "theta", "density")
slc.annotate_grids()
slc.set_figure_size(15)
slc.set_buff_size(1600)
slc.save("xrb_slice_grid.png")
