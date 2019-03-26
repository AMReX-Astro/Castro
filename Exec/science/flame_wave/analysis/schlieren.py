import yt
import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__ == "__main__":

    if len(sys.argv) < 2:
        sys.exit("Please provide a plotfile!")

    ds = yt.load(sys.argv[1])

    level = 0
    dims = ds.domain_dimensions * ds.refine_by**level

    sq = ds.covering_grid(
        level, left_edge=[0, 0, 0], dims=dims, fields=["density"])

    def schlieren(data, cutoff=-16):
        """
        Calculates the Schlieren S = ln|(del^2 rho) / rho|. The cutoff is used to
        filter out low amplitude noise.
        """
        dens = data.field_data[('gas', 'density')]

        nr, nz, nl = dens.shape
        dens = dens[:, :, 0]

        coords = data.fcoords

        r = coords[:, 0].reshape((nr, nz, nl))[:, :, 0]
        z = coords[:, 1].reshape((nr, nz, nl))[:, :, 0]

        lapl = np.zeros_like(dens)

        # assume these are constant, which is at least true for this data
        dr = r[5, 5] - r[4, 5]
        dz = z[5, 5] - z[5, 4]

        # calculate the laplacian
        lapl[1:-1, :] = 1 / (2 * dr * r[1:-1]) * (dens[2:, :] - dens[:-2, :]) + \
            (dens[2:, :] + dens[:-2, :] - 2 * dens[1:-1, :]) / dr**2

        lapl[:, 1:-1] += 1 / dz**2 * \
            (dens[:, 2:] + dens[:, :-2] - 2 * dens[:, 1:-1])

        lapl[:, :] /= dens

        S = np.log(np.abs(lapl))
        S[S < cutoff] = cutoff

        return S

    fig, ax = plt.subplots(figsize=(18, 8))

    cax = ax.imshow(schlieren(sq, cutoff=-16).T,
                    origin='lower', cmap=plt.cm.bone_r)
    ax.set_xlabel(r'$r$')
    ax.set_ylabel(r'$z$')
    ax.set_title(r'$\mathcal{S} = \ln|(\nabla^2 \rho) / \rho|$')

    cbar = fig.colorbar(cax, orientation='horizontal')
    plt.show()

    fig.savefig("schlieren.png")
