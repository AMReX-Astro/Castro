"""
This script creates a particle file for the test problem.
Note it's not very clever - if you change problo and probhi in the inputs file,
you'll need to change them here as well.
"""

import numpy as np

outfile_name = "particle_file"

# number of particles
n_particles = 1000

# copy these from the inputs file
problo = np.array([0, 0])
probhi = np.array([1.84320e5, 3.072e4])

#creates a grid of particle positions at 20 different heights spanning the lower half of the y-domain
#and at 50 different lengths spanning the lower half of the x-domain
xs = np.zeros((n_particles, 2))
yps = np.linspace(0,n_particles,num=21,endpoint=True)
#I'm now setting custom heights, but could also use problo and probhi defined above
hs = np.linspace(0.02e5,0.1e5,num=20,endpoint=True)
ls = np.linspace(0,1.5e5,num=50,endpoint=True)

for i in range(0,20):
    xs[int(yps[i]):int(yps[i+1])+1,1] = hs[i]
for y in yps[0:20]:
    for i in range(0,50):
        xs[i + int(y),0] = ls[i]

with open(outfile_name, 'w') as outfile:
    outfile.write(f"{n_particles}\n")

    for pos in xs:
        outfile.write(f"{pos[0]} {pos[1]}\n")