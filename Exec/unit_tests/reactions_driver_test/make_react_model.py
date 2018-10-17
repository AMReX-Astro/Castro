#!/usr/bin/env python
"""
Build a model file for the reactions unit test with
specified density, temperature, and composition values.
"""
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-dmin', '--dens_min', type=float, required=True,
                    help='Minimum density for grid in g/cm^3.')
parser.add_argument('-dmax', '--dens_max', type=float, required=True,
                    help='Maximum density for grid in g/cm^3.')
parser.add_argument('-nd', '--num_dens', type=int, required=True,
                    help='Number of points in density.')
parser.add_argument('-tmin', '--temp_min', type=float, required=True,
                    help='Minimum temperature for grid in K.')
parser.add_argument('-tmax', '--temp_max', type=float, required=True,
                    help='Maximum temperature for grid in K.')
parser.add_argument('-nt', '--num_temp', type=int, required=True,
                    help='Number of points in temperature.')
parser.add_argument('-xin', '--xin_file', type=str, required=True,
                    help='Name of input file containing composition grid.')
parser.add_argument('-o', '--outfile', type=str, default='model',
                    help='Name of output model file. Defaults to "model"')
args = parser.parse_args()

def get_composition():
    spec_re = re.compile('\A#\s*(\S*)\Z')
    
    species_names = []
    xin_values = []

    f = open(args.xin_file, 'r')
    lines = []
    for l in f:
        ls = l.strip()
        if ls:
            lines.append(ls)
    f.close()

    try:
        assert(len(lines) % 2 == 0)
    except:
        print('xin_file should have an even number of non-blank lines.')
        raise

    num_spec_zones = -1
    while(lines):
        sline = lines.pop(0)
        xline = lines.pop(0)
        match = spec_re.match(sline)
        try:
            assert(match)
        except:
            print('xin_file species line improperly formatted: {}'.format(sline))
            raise
        spec_name = match.group(1)
        xvals = [max(float(f.replace('d','e').replace('D','E')), 1.0e-10) for f in xline.split()]
        species_names.append(spec_name)
        xin_values.append(np.array(xvals))
        if num_spec_zones == -1:
            num_spec_zones = len(xvals)
        else:
            try:
                assert(num_spec_zones == len(xvals))
            except:
                print('Species must have the same number of mass fraction values in xin_file.')
                raise
    return species_names, xin_values, num_spec_zones

def get_logspace_grid(fmin, fmax, fnum):
    # Do this manually to make sure its the same
    # as in our Microphysics unit test.
    values_grid = []
    dlogf = (np.log10(fmax) - np.log10(fmin))/(fnum - 1)
    for i in range(fnum):
        values_grid.append(10.0**(np.log10(fmin) + float(i)*dlogf))
    return values_grid

def write_output(dens_grid, temp_grid, spec_names, spec_xin, num_spec_zones):
    f = open(args.outfile, 'w')

    # Write header
    num_zones = args.num_dens * args.num_temp * num_spec_zones
    f.write('# npts = {}\n'.format(num_zones))
    f.write('# num of variables = {}\n'.format(3+len(spec_names)))
    f.write('# density\n')
    f.write('# temperature\n')
    f.write('# pressure\n')
    for sn in spec_names:
        f.write('# {}\n'.format(sn))

    # Write zone data
    for k in range(num_spec_zones):
        for j in range(args.num_temp):
            for i in range(args.num_dens):
                line = '0.0 {} {} 0.0 '.format(dens_grid[i], temp_grid[j])
                xk = ['{}'.format(x[k]) for x in spec_xin]
                sline = ' '.join(xk)
                line = line + sline + '\n'
                f.write(line)
    f.close()

if __name__ == '__main__':
    spec_names, spec_xin, num_spec_zones = get_composition()
    dens_grid = get_logspace_grid(args.dens_min, args.dens_max, args.num_dens)
    temp_grid = get_logspace_grid(args.temp_min, args.temp_max, args.num_temp)
    write_output(dens_grid, temp_grid, spec_names, spec_xin, num_spec_zones)
