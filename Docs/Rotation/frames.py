import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# from: https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

fig = plt.figure()

ax = fig.gca(projection="3d")

C = (0, 0, 0)
dx = 2.0

Cprime = (4, -2, 0)

p = (2, -3, 0)

# draw C axes

ax.plot([C[0], C[0]+dx], [C[1], C[1]], [C[2], C[2]], color="C0")
ax.plot([C[0], C[0]], [C[1], C[1]-dx], [C[2], C[2]], color="C0")
ax.plot([C[0], C[0]], [C[1], C[1]], [C[2], C[2]+dx], color="C0")

eps = 0.025
ax.text(C[0]+eps, C[1]+eps, C[2]+eps, "$C$", color="C0")

# draw Cprime

ax.plot([Cprime[0], Cprime[0]+dx], [Cprime[1], Cprime[1]], [Cprime[2], Cprime[2]], color="C1")
ax.plot([Cprime[0], Cprime[0]], [Cprime[1], Cprime[1]-dx], [Cprime[2], Cprime[2]], color="C1")
ax.plot([Cprime[0], Cprime[0]], [Cprime[1], Cprime[1]], [Cprime[2], Cprime[2]+dx], color="C1")

ax.text(Cprime[0]+eps, Cprime[1]+eps, Cprime[2]+eps, r"$\tilde{C}$", color="C1")

ax.scatter(p[0], p[1], p[2], color="r")

ax.text(p[0]+eps, p[1]+eps, p[2]+eps, r"$P$", color="r")


a = Arrow3D([C[0], p[0]], [C[1], p[1]], [C[2], p[2]], mutation_scale=20, 
            lw=1, arrowstyle="-|>", color="C0", ls=":")
ax.add_artist(a)

ax.text(0.5*(C[0] + p[0]), 0.5*(C[1] + p[1]), 0.5*(C[2] + p[2]), r"$\mathbf{r}$", color="C0")

a = Arrow3D([Cprime[0], p[0]], [Cprime[1], p[1]], [Cprime[2], p[2]], mutation_scale=20, 
            lw=1, arrowstyle="-|>", color="C1", ls=":")
ax.add_artist(a)

ax.text(0.5*(Cprime[0] + p[0]), 0.5*(Cprime[1] + p[1]), 0.5*(Cprime[2] + p[2]), r"$\tilde{\mathbf{r}}$", color="C1")

a = Arrow3D([C[0], Cprime[0]], [C[0], Cprime[1]], [C[0], Cprime[2]], mutation_scale=20, 
            lw=1, arrowstyle="-|>", color="k", ls=":")
ax.add_artist(a)

ax.text(0.5*(Cprime[0] + C[0]), 0.5*(Cprime[1] + C[1]), 0.5*(Cprime[2] + C[2]), r"$\mathbf{l}$", color="k")

ax.set_aspect('equal')
ax.auto_scale_xyz([-1,4], [-1,4], [-2,3])

fig.set_size_inches(8.0, 8.0)


plt.axis("off")
plt.tight_layout()
plt.savefig("frames.png", dpi=200)
