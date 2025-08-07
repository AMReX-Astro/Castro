import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'serif'
plt.rcParams["axes.labelsize"] = 16
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['savefig.bbox']='tight'

def ext_data(file):
    from diag_parser import deduplicate, read_diag_file
    temp = read_diag_file(file)
    data = deduplicate(temp)
    return data


def variables(data, per:int, end:int=None):
    if not end:
        time = data["TIME"][::per]
        mass = data["MASS"][::per]
        temp = data["MAXIMUM_TEMPERATURE"][::per]
    else:
        time = data["TIME"][:end:per]
        mass = data["MASS"][:end:per]
        temp = data["MAXIMUM_TEMPERATURE"][:end:per]
    return time, mass, temp

data_r04 = ext_data('grid_diag.out')
per = 10000
time_r04, mass_r04, temperature_r04 = variables(data_r04, per=per)


# i_max_peak_04 = temperature_r04.argmax()
# print(f"The TNR max temperature at 0.4km is: {temperature_r04[i_max_peak_04]:.3e}")
# print(f"The TNR max temperature time at 0.4km is: {time_r04[i_max_peak_04]:.2f}")
# i_max_peak_08 = temperature_r08.argmax()
# print(f"The TNR max temperature at 0.8km is:: {temperature_r08[i_max_peak_08]:.3e}")
# print(f"The CEI TNR time at 0.8km is: {time_r08[i_max_peak_08]:.2f}")

fig = plt.figure(figsize=(8,6))

# num_plots = 5
# colormap = plt.cm.gist_ncar
# plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, num_plots))))

ax = fig.add_subplot(111)
ax.plot(time_r04, temperature_r04/1.0e8, label=r"Model A4_drive", linestyle='--')

ax.set_ylim([0.5, 0.9])
#ax.set_title("Max Temperature evolution")
ax.set_xlabel(r"$time\quad[s]$")
ax.set_ylabel(r"$T_{\mathrm{peak}}\quad\times 10^8\,[K]$")
ax.legend(loc="upper left")

fig.set_size_inches(8, 6)
ax.tick_params(labelsize=16)
plt.tight_layout()

fig.savefig("temp_t_04_ext_drive.pdf", bbox_inches="tight")

# fig1 = plt.figure(figsize=(8,6))
# ax1 = fig1.add_subplot(111)
# ax1.plot(time_r04, mass_r04, label=r"Model B4")
# ax1.plot(time_r08, mass_r08, label=r"Model B8")
# ax1.set_ylim([2.90e20, 5.025e20])
# ax1.set_title("Total Mass evolution")
# ax1.set_xlabel(r"$time\quad[s]$")
# ax1.set_ylabel(r"$M_{\mathrm{total}}\quad[g]$")
# ax1.legend(loc="lower left")

# fig1.savefig("mass_t_04_ext.png")