import yt

ds = yt.load("smallplt10598")

field = "density"
slc = yt.SlicePlot(ds, "z", field)
slc.set_buff_size((2000, 2000))
#slc.set_cmap(field, "arbre")
#slc.set_cmap(field, "RdPu")
slc.set_cmap(field, "YlGn")
slc.set_zlim(field, 1.e4, 1.5e7)

#slc.annotate_grids()
slc.save("slice_grid.png")
