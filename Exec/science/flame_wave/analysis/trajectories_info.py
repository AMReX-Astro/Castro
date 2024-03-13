"""
This script processes the Timestamp outputs from a Castro simulation with tracer particles enabled.
It loops through each Timestamp line-by-line and grabs unique particle identifiers (the particle ID and cpu number)
And stores the information at all timesteps for each particle into a unique txt file

After looping through all the Timestamps, it loops through all the new txt files and sorts the contents by timestep.
"""

import numpy as np
import glob
import pandas as pd
import subprocess

#trying to grab specific particle trajectories
Timestamps = sorted(glob.glob("/global/project/projectdirs/m3018/XRB/flame_wave_particles/fullrun1/particle_dir/Timestamp_*"))
#Timestamps = ["/global/project/projectdirs/m3018/XRB/flame_wave_particles/fullrun1/particle_dir/Timestamp_10","/global/project/projectdirs/m3018/XRB/flame_wave_particles/fullrun1/particle_dir/Timestamp_62"]

#create a set to track the number of unique particle identifiers
particles = set()
#loop through each timestamp
for t in Timestamps:
    print(t)
    #create a dictionary to store particle identifiers and info
    d = {}
    #loop through each line of a timestamp
    with open(t) as infile:
        for line in infile:
            #breaking up each line into an indexed array
            line = line.split()
            #store the pid information for this line
            pid = line[0]
            cpu = line[1]
            #line = [np.fromstring(line, sep = " ")]
            #check if a dict key for this unique particle identifier exists
            if f"{pid}_{cpu}" in d:
                #if it exists, append the current line of information to the key's value (numpy array)
                d[f"{pid}_{cpu}"].append(line[2:-1])
            else:
                #if the key/identifier doesn't yet exist, create the key and assign it a value
                d[f"{pid}_{cpu}"] = [line[2:-1]]
    #loop through the particle identifiers
    for key in d:
        #check if the pid key has already been encountered. If not, add it to the particles set to count how many particles have been encountered
        particles.add(key)
        #open each particle identifier text file, or create it if it doesn't exist yet, then append the key's array value
        with open(f"./{key}_info.txt",'a+') as outfile:
            np.savetxt(outfile, np.array(d[key]), delimiter=" ",fmt='%s')
            outfile.write('\n')
    #check how many unique particles have been identified after each timestamp
    print(len(particles))

files_to_sort = sorted(glob.glob("./*.txt"))
print("sorting files...")
for f in files_to_sort:
    #sorts lines in files: -g = sort scientific notation, -k 3 = use third column (time), -o = specify output filename
    subprocess.run(["sort","-g","-k 3",f,"-o",f], shell=False)