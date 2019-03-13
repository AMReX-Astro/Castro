#!/bin/bash
#SBATCH -J react_conv
#SBATCH -A m3018
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -t 12:00:00

# -N = number of nodes = # of MPI * number of threads / 24 (cores / node)

export OMP_NUM_THREADS=1


slack_job_start.py "starting reacting_convergence problem @ NERSC; ${restartFile}" @michael


# Edison has 2 Intel "IvyBridge" processors each with 12 cores.  Each processor
# connects to 1/2 of the node's memory, and is called a NUMA node.  So there
# are 2 NUMA nodes per Edison node.

# -n is the number of MPI tasks
# -c sets the number of OMP threads per MPI -- this should be the same as OMP_NUM_THREADS

# Strang
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.64 amr.plot_file=react_converge_64_strang_plt castro.time_integration_method=0
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.128 amr.plot_file=react_converge_128_strang_plt castro.time_integration_method=0
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.256 amr.plot_file=react_converge_256_strang_plt castro.time_integration_method=0
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.512 amr.plot_file=react_converge_512_strang_plt castro.time_integration_method=0

# SDC
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.64 amr.plot_file=react_converge_64_sdc_plt castro.sdc_order=2 castro.time_integration_method=2
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.128 amr.plot_file=react_converge_128_sdc_plt castro.sdc_order=2 castro.time_integration_method=2
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.256 amr.plot_file=react_converge_256_sdc_plt castro.sdc_order=2 castro.time_integration_method=2
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.512 amr.plot_file=react_converge_512_sdc_plt castro.sdc_order=2 castro.time_integration_method=2

# SDC-4
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.64 amr.plot_file=react_converge_64_sdc4_plt castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.128 amr.plot_file=react_converge_128_sdc4_plt castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.256 amr.plot_file=react_converge_256_sdc4_plt castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1
srun -n 24 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs.512 amr.plot_file=react_converge_512_sdc4_plt castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1





