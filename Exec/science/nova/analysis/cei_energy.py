import yt
import os
import re
import numpy as np
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system
import matplotlib.pyplot as plt

yt.enable_parallelism()

plt.rcParams['font.size'] = 21
plt.rcParams['font.family'] = 'serif'
plt.rcParams["axes.labelsize"] = 21
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['savefig.bbox']='tight'

curr_dir = os.getcwd()
plt_dir_a04 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_a04.sort(key=lambda x: int(x[3:]))

dir_paths_a04 = []
for dir in plt_dir_a04:
    dir_paths_a04.append(os.path.join(curr_dir, dir))

index_a04 = dict(zip(plt_dir_a04, [i for i in range(len(plt_dir_a04))]))

curr_dir = "../nova_t7_0.8km_ext_pert"
plt_dir_a08 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_a08.sort(key=lambda x: int(x[3:]))

dir_paths_a08 = []
for dir in plt_dir_a08:
    dir_paths_a08.append(os.path.join(curr_dir, dir))

index_a08 = dict(zip(plt_dir_a08, [i for i in range(len(plt_dir_a08))]))

curr_dir = "../nova_t7_0.2km"
plt_dir_b02 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_b02.sort(key=lambda x: int(x[3:]))

dir_paths_b02 = []
for dir in plt_dir_b02:
    dir_paths_b02.append(os.path.join(curr_dir, dir))

index_b02 = dict(zip(plt_dir_b02, [i for i in range(len(plt_dir_b02))]))

curr_dir = "../nova_t7_0.4km"
plt_dir_b04 = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir_b04.sort(key=lambda x: int(x[3:]))

dir_paths_b04 = []
for dir in plt_dir_b04:
    dir_paths_b04.append(os.path.join(curr_dir, dir))

index_b04 = dict(zip(plt_dir_b04, [i for i in range(len(plt_dir_b04))]))

ts_a08 = yt.DatasetSeries(dir_paths_a08)
ts_a04 = yt.DatasetSeries(dir_paths_a04)

ts_b04 = yt.DatasetSeries(dir_paths_b04[0:138])
ts_b02 = yt.DatasetSeries(dir_paths_b02[0:267])

N_a08 = len(plt_dir_a08)
time_a08 = np.zeros(N_a08)
qvalue_a08 = np.zeros(N_a08)

N_a04 = len(plt_dir_a04)
time_a04 = np.zeros(N_a04)
qvalue_a04 = np.zeros(N_a04)

N_b04 = len(plt_dir_b04[0:138])
time_b04 = np.zeros(N_b04)
qvalue_b04 = np.zeros(N_b04)

N_b02 = len(plt_dir_b02[0:267])
time_b02 = np.zeros(N_b02)
qvalue_b02 = np.zeros(N_b02)

comm = communication_system.communicators[-1]

for ds in ts_a04.piter():

    id_i = index_a04[ds.basename]

    sp = ds.all_data()

    max_enuc = sp.max("enuc")
    time_a04[id_i] = ds.current_time.value
    qvalue_a04[id_i] = max_enuc

    ds.index.clear_all_data()

comm.barrier()

for ds in ts_a08.piter():

    id_i = index_a08[ds.basename]

    sp = ds.all_data()

    max_enuc = sp.max("enuc")
    time_a08[id_i] = ds.current_time.value
    qvalue_a08[id_i] = max_enuc

    ds.index.clear_all_data()

for ds in ts_b02.piter():

    id_i = index_b02[ds.basename]

    sp = ds.all_data()

    max_enuc = sp.max("enuc")
    time_b02[id_i] = ds.current_time.value
    qvalue_b02[id_i] = max_enuc

    ds.index.clear_all_data()

for ds in ts_b04.piter():

    id_i = index_b04[ds.basename]

    sp = ds.all_data()

    max_enuc = sp.max("enuc")
    time_b04[id_i] = ds.current_time.value
    qvalue_b04[id_i] = max_enuc

    ds.index.clear_all_data()

time_a04[:] = comm.mpi_allreduce(time_a04[:], op="sum")
qvalue_a04[:] = comm.mpi_allreduce(qvalue_a04[:], op="sum")

time_a08[:] = comm.mpi_allreduce(time_a08[:], op="sum")
qvalue_a08[:] = comm.mpi_allreduce(qvalue_a08[:], op="sum")

time_b02[:] = comm.mpi_allreduce(time_b02[:], op="sum")
qvalue_b02[:] = comm.mpi_allreduce(qvalue_b02[:], op="sum")

time_b04[:] = comm.mpi_allreduce(time_b04[:], op="sum")
qvalue_b04[:] = comm.mpi_allreduce(qvalue_b04[:], op="sum")

if yt.is_root():

    per = 10
    t_a04 = time_a04[::per]
    t_a08 = time_a08[::per]
    t_b02 = time_b02[::per]
    t_b04 = time_b04[::per]

    q_a04 = qvalue_a04[::per]
    q_a08 = qvalue_a08[::per]
    q_b02 = qvalue_b02[::per]
    q_b04 = qvalue_b04[::per]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t_a04, q_a04, label="Model A4", linestyle='--')
    ax.plot(t_a08, q_a08, label="Model A8", linestyle='--')
    ax.plot(t_b02, q_b02, label="Model B2", linestyle='-', linewidth=3.0)
    ax.plot(t_b04, q_b04, label="Model B4", linestyle='-', linewidth=3.0)

    ax.set_xlabel(r'$t\,[s]$')
    ax.set_ylabel(r'$\dot{e}_{\mathrm{nuc},\mathrm{max}}\,[\mathrm{erg}\cdot \mathrm{g}^{-1}\,\mathrm{s}^{-1}]$')
    ax.set_yscale('log')
    ax.set_ylim(1.0e13, 1.0e19)
    ax.legend(loc='upper left')

    fig.set_size_inches(8, 6)
    ax.tick_params(labelsize=21)
    plt.tight_layout()

    fig.savefig('q_series_04_ext.pdf', bbox_inches="tight")
