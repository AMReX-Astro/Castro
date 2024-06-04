import matplotlib.pyplot as plt
import numpy as np

import yt

ref_pf = "det_x_3km_plt00133"
compare_pf = "det_x_24km_plt00038"

field = "Temp"

fig, ax = plt.subplots()

legend_elem = []
legend_labels = []

for pf, label, color in [(ref_pf, "3 km", "C0"),
                         (compare_pf, "24 km", "C1")]:

    ds = yt.load(pf)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    ref = int(np.prod(ds.ref_factors[0:ds.index.max_level]))

    data = ds.covering_grid(ds.index.max_level,
                            left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions*ref, fields=field)

    print(data.shape)

    ct = ax.contour(data[field].d[:,:,0].T, levels=4, colors=color,
                    extent=[xmin, xmax, ymin, ymax])

    legend_elem.append(ct.legend_elements()[0][0])
    legend_labels.append(label)

ax.set_aspect("equal")

ax.legend(legend_elem, legend_labels)
fig.savefig("contour.png")

