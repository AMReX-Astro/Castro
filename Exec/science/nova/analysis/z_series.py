import yt
import os
import re
import numpy as np
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system
import matplotlib.pyplot as plt
from yt.units import dimensions

plt.rcParams['font.size'] = 21
plt.rcParams['font.family'] = 'serif'
plt.rcParams["axes.labelsize"] = 21
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['savefig.bbox']='tight'

yt.enable_parallelism()

def _zenv_a(field, data):

    fields = [ ('boxlib', 'X(B8)'),
                ('boxlib', 'X(Be7)'),
                ('boxlib', 'X(C12)'),
                ('boxlib', 'X(C13)'),
                ('boxlib', 'X(F17)'),
                ('boxlib', 'X(F18)'),
                ('boxlib', 'X(N13)'),
                ('boxlib', 'X(N14)'),
                ('boxlib', 'X(N15)'),
                ('boxlib', 'X(O14)'),
                ('boxlib', 'X(O15)'),
                ('boxlib', 'X(O16)'),
                ('boxlib', 'X(O17)') ]
    zenv_a = 0.0

    for f in fields:
        zenv_a += data[f]
    return zenv_a

def _zenv_b(field, data):

    fields = [ ('boxlib', 'X(b8)'),
                ('boxlib', 'X(be7)'),
                ('boxlib', 'X(c12)'),
                ('boxlib', 'X(c13)'),
                ('boxlib', 'X(f17)'),
                ('boxlib', 'X(f18)'),
                ('boxlib', 'X(n13)'),
                ('boxlib', 'X(n14)'),
                ('boxlib', 'X(n15)'),
                ('boxlib', 'X(o14)'),
                ('boxlib', 'X(o15)'),
                ('boxlib', 'X(o16)'),
                ('boxlib', 'X(o17)') ]
    zenv_b = 0.0

    for f in fields:
        zenv_b += data[f]
    return zenv_b

#=============================
curr_dir = os.getcwd()
plt_dir_a04 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_a04.sort(key=lambda x: int(x[3:]))

dir_paths_a04 = []
for dir in plt_dir_a04:
    dir_paths_a04.append(os.path.join(curr_dir, dir))

index_a04 = dict(zip(plt_dir_a04, [i for i in range(len(plt_dir_a04))]))

#=============================
curr_dir = "../nova_t7_0.8km_ext_pert"
plt_dir_a08 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_a08.sort(key=lambda x: int(x[3:]))

dir_paths_a08 = []
for dir in plt_dir_a08:
    dir_paths_a08.append(os.path.join(curr_dir, dir))

index_a08 = dict(zip(plt_dir_a08, [i for i in range(len(plt_dir_a08))]))

#=============================

curr_dir = "../nova_t7_0.2km"
plt_dir_b02 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_b02.sort(key=lambda x: int(x[3:]))

dir_paths_b02 = []
for dir in plt_dir_b02:
    dir_paths_b02.append(os.path.join(curr_dir, dir))

index_b02 = dict(zip(plt_dir_b02, [i for i in range(len(plt_dir_b02))]))

#=============================

curr_dir = "../nova_t7_0.4km"
plt_dir_b04 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_b04.sort(key=lambda x: int(x[3:]))

dir_paths_b04 = []
for dir in plt_dir_b04:
    dir_paths_b04.append(os.path.join(curr_dir, dir))

index_b04 = dict(zip(plt_dir_b04, [i for i in range(len(plt_dir_b04))]))

#=============================

ts_a08 = yt.DatasetSeries(dir_paths_a08)
ts_a04 = yt.DatasetSeries(dir_paths_a04)

ts_b04 = yt.DatasetSeries(dir_paths_b04[0:138])
ts_b02 = yt.DatasetSeries(dir_paths_b02[0:267])

N_a08 = len(plt_dir_a08)
time_a08 = np.zeros(N_a08)
z_value_a08 = np.zeros(N_a08)

N_a04 = len(plt_dir_a04)
time_a04 = np.zeros(N_a04)
z_value_a04 = np.zeros(N_a04)

N_b02 = len(plt_dir_b02[0:267])
time_b02 = np.zeros(N_b02)
z_value_b02 = np.zeros(N_b02)

N_b04 = len(plt_dir_b04[0:138])
time_b04 = np.zeros(N_b04)
z_value_b04 = np.zeros(N_b04)

comm = communication_system.communicators[-1]

# =========================
for ds in ts_a08.piter():

    ds.add_field(
        ("gas","zenv_a"),
        function=_zenv_a,
        sampling_type="local",
        units="auto",
        take_log=False,
        display_name=r"Z_{\mathrm{CNO}}",
        dimensions=dimensions.dimensionless,
    )

    id_i = index_a08[ds.basename]

    env = ds.region(center=[1.536e+08, 7.680e+07, 5.000e-01], left_edge=[0., 0.2259e+08, 0.], right_edge=[3.072e+08, 1.536e+08, 1.000e+00])
    zmean = env.mean(("gas", "zenv_a"), weight=("gas", "mass"))

    time_a08[id_i] = ds.current_time.value
    z_value_a08[id_i] = zmean

    ds.index.clear_all_data()

comm.barrier()

for ds in ts_a04.piter():

    ds.add_field(
        ("gas","zenv_a"),
        function=_zenv_a,
        sampling_type="local",
        units="auto",
        take_log=False,
        display_name=r"Z_{\mathrm{CNO}}",
        dimensions=dimensions.dimensionless,
    )

    id_i = index_a04[ds.basename]

    env = ds.region(center=[1.536e+08, 7.680e+07, 5.000e-01], left_edge=[0., 0.2259e+08, 0.], right_edge=[3.072e+08, 1.536e+08, 1.000e+00])
    zmean = env.mean(("gas", "zenv_a"), weight=("gas", "mass"))

    time_a04[id_i] = ds.current_time.value
    z_value_a04[id_i] = zmean

    ds.index.clear_all_data()

comm.barrier()

# ======================

for ds in ts_b02.piter():

    ds.add_field(
        ("gas","zenv_b"),
        function=_zenv_b,
        sampling_type="local",
        units="auto",
        take_log=False,
        display_name=r"Z_{\mathrm{CNO}}",
        dimensions=dimensions.dimensionless,
    )

    id_i = index_b02[ds.basename]

    env = ds.region(center=[1.024e+08, 5.120e+07, 5.000e-01], left_edge=[0., 0.2259e+08, 0.], right_edge=[2.048e+08, 1.024e+08, 1.000e+00])
    zmean = env.mean(("gas", "zenv_b"), weight=("gas", "mass"))

    time_b02[id_i] = ds.current_time.value
    z_value_b02[id_i] = zmean

    ds.index.clear_all_data()

comm.barrier()

#====================

for ds in ts_b04.piter():

    ds.add_field(
        ("gas","zenv_b"),
        function=_zenv_b,
        sampling_type="local",
        units="auto",
        take_log=False,
        display_name=r"Z_{\mathrm{CNO}}",
        dimensions=dimensions.dimensionless,
    )

    id_i = index_b04[ds.basename]

    env = ds.region(center=[1.024e+08, 5.120e+07, 5.000e-01], left_edge=[0., 0.2259e+08, 0.], right_edge=[2.048e+08, 1.024e+08, 1.000e+00])
    zmean = env.mean(("gas", "zenv_b"), weight=("gas", "mass"))

    time_b04[id_i] = ds.current_time.value
    z_value_b04[id_i] = zmean

    ds.index.clear_all_data()

comm.barrier()

time_a04[:] = comm.mpi_allreduce(time_a04[:], op="sum")
z_value_a04[:] = comm.mpi_allreduce(z_value_a04[:], op="sum")

time_a08[:] = comm.mpi_allreduce(time_a08[:], op="sum")
z_value_a08[:] = comm.mpi_allreduce(z_value_a08[:], op="sum")

time_b02[:] = comm.mpi_allreduce(time_b02[:], op="sum")
z_value_b02[:] = comm.mpi_allreduce(z_value_b02[:], op="sum")

time_b04[:] = comm.mpi_allreduce(time_b04[:], op="sum")
z_value_b04[:] = comm.mpi_allreduce(z_value_b04[:], op="sum")

if yt.is_root():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$t\,[\mathrm{s}]$')
    ax.set_ylabel(r'$Z_{\mathrm{env}}$')
    ax.plot(time_a04,z_value_a04, label="Model A4", linestyle='--')
    ax.plot(time_a08,z_value_a08, label="Model A8", linestyle='--')
    ax.plot(time_b02,z_value_b02, label="Model B2", linestyle='-', linewidth=3.0)
    ax.plot(time_b04,z_value_b04, label="Model B4", linestyle='-', linewidth=3.0)
    ax.set_ylim(0.0, 0.4)
    ax.legend()

    fig.set_size_inches(8, 6)
    ax.tick_params(labelsize=21)
    plt.tight_layout()

    plt.savefig('z_t_04_ext.pdf',  bbox_inches="tight")
