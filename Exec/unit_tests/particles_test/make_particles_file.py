"""
This script creates a particle file for the test problem.
Note it's not very clever - if you change problo and probhi in the inputs file,
you'll need to change them here as well.
"""

import numpy as np

outfile_name = "particle_file"

# number of particles
n_particles = 20

# copy these from the inputs file
problo = np.array([0, 0])
probhi = np.array([1, 1])

# maximum distance from center
max_R = np.max(0.5 * (probhi - problo))

dr_part = max_R / n_particles
dtheta_part = 2 * np.pi / n_particles

xs = np.zeros((n_particles, 2))

theta = np.linspace(0, n_particles, num=n_particles,
                    endpoint=False) * dtheta_part
r = (np.linspace(0, n_particles, num=n_particles, endpoint=False) + 0.5) * dr_part

xs[:, 0] = r * np.cos(theta)
xs[:, 1] = r * np.sin(theta)

xs[:, :] += 0.5 * (problo + probhi)[np.newaxis, :]

with open(outfile_name, 'w') as outfile:
    outfile.write(f"{n_particles}\n")

    for pos in xs:
        outfile.write(f"{pos[0]} {pos[1]}\n")
