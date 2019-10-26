import yt

ds = yt.load("diffuse_plt00147")

slc = yt.SlicePlot(ds, "theta", "Temp")
slc.annotate_grids()
slc.set_figure_size(16)
slc.set_buff_size(1600)
slc.set_font_size(64)
slc.set_cmap("Temp", "plasma_r")
slc.set_log("Temp", False)
slc.set_axes_unit("cm")

slc.save("diffusion_temp_amr.pdf")
