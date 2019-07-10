import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("summit_scaling_may19.txt")

for i in range (0,40,5):
    #print(i)
    nodes = data[i:i+5,0]
    avg_zone = data[i:i+5,2]

    if i in [0,5,20,25]:
        if data[i,1]==256 : #blue for 256^3 and orange for 512^3 when using gpu
            color = "C0"
        else:
            color = "C1"
    else:
        if data[i,1]==256 : #skyblue for 256^3 and red for 512^3 when using MPI+OpenMP
            color = "C9"
        else:
            color = "C3"



    if i in [5,15,25,35]: #amr triangle marker, no amr circle
        marker ="^"
    else:
        marker ="o"

    #plt.scatter(nodes, avg_zone)
    plt.plot(nodes,avg_zone, marker+color, ls=":")

plt.xlabel("number of nodes")
plt.ylabel(r"Avg # of zones advanced/ $\mu$s")

#legends
legs = []
legnames = []
legs.append(plt.Line2D((0,1),(0,0), color = "C0"))
legnames.append(r"$256^3$ gpu")
legs.append(plt.Line2D((0,1),(0,0), color = "C1"))
legnames.append(r"$512^3$ gpu")
legs.append(plt.Line2D((0,1),(0,0), color = "C9"))
legnames.append(r"$256^3$ MPI+OMP")
legs.append(plt.Line2D((0,1),(0,0), color = "C3"))
legnames.append(r"$512^3$ MPI+OMP")
legs.append(plt.Line2D((0,1),(0,0), color="k",
                       marker="o", markeredgecolor="k", markerfacecolor="k", linestyle="none"))
legnames.append("no AMR")
legs.append(plt.Line2D((0,1),(0,0), color="k",
                       marker="^", markeredgecolor="k", markerfacecolor="k",  linestyle="none"))
legnames.append("base + one 4x level")

plt.legend(legs, legnames, frameon=False, fontsize="8", numpoints=1, loc=0, ncol=3)

plt.savefig("summit_sedov.png", dpi=150)
