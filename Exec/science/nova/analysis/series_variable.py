import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt('grid_diag.out', skiprows=2)

time = file[: ,1] 
temperature = file[:, 20] 

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.plot(time, temperature)
ax.set_ylim([3.0e5, 2.0e8])
ax.set_xlabel(r"$time\quad[s]$")
ax.set_ylabel(r"$T\quad[K]$")

fig.savefig("temp_t.png")
