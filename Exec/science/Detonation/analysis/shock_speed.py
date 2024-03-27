import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np


def get_T_profile(plotfile):

    ds = yt.load(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])

    return time, x_coord, temp


def find_x_for_T(x, T, T_0=2.e9):
    """ given a profile x(T), find the x_0 that corresponds to T_0 """

    # our strategy here assumes that the hot ash is in the early part
    # of the profile.  We then find the index of the first point where
    # T drops below T_0
    idx = np.where(T < T_0)[0][0]

    T1 = T[idx-1]
    x1 = x[idx-1]

    T2 = T[idx]
    x2 = x[idx]

    slope = (x2 - x1)/(T2 - T1)

    return x1 + slope*(T_0 - T1)


if __name__ == "__main__":

    try:
        plotfile = sys.argv[1:]
    except:
        sys.exit("must supply plotfile(s)")

    t = []
    v = []
    for n, p in enumerate(plotfile):
        time, x, T = get_T_profile(p)
        xpos = find_x_for_T(x, T)

        if n == 0:
            xpos_old = xpos
            time_old = time
        else:
            # difference with the previous file to find the det speed
            # note: the corresponding time is centered in the interval
            v.append((xpos - xpos_old)/(time - time_old))
            t.append(0.5*(time_old + time))

            xpos_old = xpos
            time_old = time

    plt.plot(t, v)
    ax = plt.gca()
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig("profile.png")

