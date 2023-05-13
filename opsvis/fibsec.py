import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Wedge

from .settings import *


# plot_fiber_section is inspired by plotSection matlab function
# written by D. Vamvatsikos available at
# http://users.ntua.gr/divamva/software.html (plotSection.zip)
def plot_fiber_section(fib_sec_list, fillflag=1,
                       matcolor=['y', 'b', 'r', 'g', 'm', 'k']):
    """Plot fiber cross-section.

    Args:
        fib_sec_list (list): list of lists in the format similar to the parameters
            for the section, layer, patch, fiber OpenSees commands

        fillflag (int): 1 - filled fibers with color specified in matcolor
            list, 0 - no color, only the outline of fibers

        matcolor (list): sequence of colors for various material tags
            assigned to fibers

    Examples:
        ::

            fib_sec_1 = [['section', 'Fiber', 1, '-GJ', 1.0e6],
                         ['patch', 'quad', 1, 4, 1,  0.032, 0.317, -0.311, 0.067, -0.266, 0.005, 0.077, 0.254],  # noqa: E501
                         ['patch', 'quad', 1, 1, 4,  -0.075, 0.144, -0.114, 0.116, 0.075, -0.144, 0.114, -0.116],  # noqa: E501
                         ['patch', 'quad', 1, 4, 1,  0.266, -0.005,  -0.077, -0.254,  -0.032, -0.317,  0.311, -0.067]  # noqa: E501
                         ]
            opsv.fib_sec_list_to_cmds(fib_sec_1)
            matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
            opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
            plt.axis('equal')
            # plt.savefig('fibsec_rc.png')
            plt.show()

    Notes:
        ``fib_sec_list`` can be reused by means of a python helper function
            ``opsvis.fib_sec_list_to_cmds(fib_sec_list_1)``

    See also:
        ``opsvis.fib_sec_list_to_cmds()``
    """

    fig, ax = plt.subplots()
    ax.set_xlabel('z')
    ax.invert_xaxis()  # To make z-axis positive to the left
    ax.set_ylabel('y')
    ax.grid(False)

    for item in fib_sec_list:

        if item[0] == 'layer':
            matTag = item[2]
            if item[1] == 'straight':
                n_bars = item[3]
                As = item[4]
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                r = np.sqrt(As / np.pi)
                Y = np.linspace(Iy, Jy, n_bars)
                Z = np.linspace(Iz, Jz, n_bars)
                for zi, yi in zip(Z, Y):
                    bar = Circle((zi, yi), r, ec='k', fc='k', zorder=10)
                    ax.add_patch(bar)
            if item[1] == 'circ':
                n_bars, As = item[3], item[4]
                yC, zC, arc_radius = item[5], item[6], item[7]
                if len(item) > 8:
                    a0_deg, a1_deg = item[8], item[9]
                    if ((a1_deg - a0_deg) >= 360. and n_bars > 0):
                        a1_deg = a0_deg + 360. - 360. / n_bars
                else:
                    a0_deg, a1_deg = 0., 360. - 360. / n_bars

                a0_rad, a1_rad = np.pi * a0_deg / 180., np.pi * a1_deg / 180.
                r_bar = np.sqrt(As / np.pi)
                thetas = np.linspace(a0_rad, a1_rad, n_bars)
                Y = yC + arc_radius * np.cos(thetas)
                Z = zC + arc_radius * np.sin(thetas)
                for zi, yi in zip(Z, Y):
                    bar = Circle((zi, yi), r_bar, ec='k', fc='k', zorder=10)
                    ax.add_patch(bar)

        if (item[0] == 'patch' and (item[1] == 'quad' or item[1] == 'quadr' or
                                  item[1] == 'rect')):
            matTag, nIJ, nJK = item[2], item[3], item[4]

            if item[1] == 'quad' or item[1] == 'quadr':
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                Ky, Kz, Ly, Lz = item[9], item[10], item[11], item[12]

            if item[1] == 'rect':
                Iy, Iz, Ky, Kz = item[5], item[6], item[7], item[8]
                Jy, Jz, Ly, Lz = Ky, Iz, Iy, Kz

            # check for convexity (vector products)
            outIJxIK = (Jy-Iy)*(Kz-Iz) - (Ky-Iy)*(Jz-Iz)
            outIKxIL = (Ky-Iy)*(Lz-Iz) - (Ly-Iy)*(Kz-Iz)
            # check if I, J, L points are colinear
            outIJxIL = (Jy-Iy)*(Lz-Iz) - (Ly-Iy)*(Jz-Iz)
            # outJKxJL = (Ky-Jy)*(Lz-Jz) - (Ly-Jy)*(Kz-Jz)

            if outIJxIK <= 0 or outIKxIL <= 0 or outIJxIL <= 0:
                print('\nWarning! Patch quad is non-convex or counter-clockwise defined or has at least 3 colinear points in line')  # noqa: E501

            IJz, IJy = np.linspace(Iz, Jz, nIJ+1), np.linspace(Iy, Jy, nIJ+1)
            JKz, JKy = np.linspace(Jz, Kz, nJK+1), np.linspace(Jy, Ky, nJK+1)
            LKz, LKy = np.linspace(Lz, Kz, nIJ+1), np.linspace(Ly, Ky, nIJ+1)
            ILz, ILy = np.linspace(Iz, Lz, nJK+1), np.linspace(Iy, Ly, nJK+1)

            if fillflag:
                Z = np.zeros((nIJ+1, nJK+1))
                Y = np.zeros((nIJ+1, nJK+1))

                for j in range(nIJ+1):
                    Z[j, :] = np.linspace(IJz[j], LKz[j], nJK+1)
                    Y[j, :] = np.linspace(IJy[j], LKy[j], nJK+1)

                for j in range(nIJ):
                    for k in range(nJK):
                        zy = np.array([[Z[j, k], Y[j, k]],
                                       [Z[j, k+1], Y[j, k+1]],
                                       [Z[j+1, k+1], Y[j+1, k+1]],
                                       [Z[j+1, k], Y[j+1, k]]])
                        poly = Polygon(zy, True, ec='k', fc=matcolor[matTag-1])
                        ax.add_patch(poly)

            else:
                # horizontal lines
                for az, bz, ay, by in zip(IJz, LKz, IJy, LKy):
                    plt.plot([az, bz], [ay, by], 'b-', zorder=1)

                # vertical lines
                for az, bz, ay, by in zip(JKz, ILz, JKy, ILy):
                    plt.plot([az, bz], [ay, by], 'b-', zorder=1)

        if item[0] == 'patch' and item[1] == 'circ':
            matTag, nc, nr = item[2], item[3], item[4]

            yC, zC, ri, re = item[5], item[6], item[7], item[8]
            a0, a1 = item[9], item[10]

            dr = (re - ri) / nr
            dth = (a1 - a0) / nc

            for j in range(nr):
                rj = ri + j * dr
                rj1 = rj + dr

                for i in range(nc):
                    thi = a0 + i * dth
                    thi1 = thi + dth
                    wedge = Wedge((zC, yC), rj1, thi, thi1, width=dr, ec='k',
                                  lw=1, fc=matcolor[matTag-1])
                    ax.add_patch(wedge)

            ax.axis('equal')


def fib_sec_list_to_cmds(fib_sec_list):
    """Reuses fib_sec_list to define fiber section in OpenSees.

    At present it is not possible to extract fiber section data from
    the OpenSees domain, this function is a workaround. The idea is to
    prepare data similar to the one the regular OpenSees commands
    (``section('Fiber', ...)``, ``fiber()``, ``patch()`` and/or
    ``layer()``) require.

    Args:
        fib_sec_list (list): is a list of fiber section data. First sub-list
        also defines the torsional stiffness (GJ).

    Warning:

    If you use this function, do not issue the regular OpenSees:
    section, Fiber, Patch or Layer commands.

    See also:

    ``opsvis.plot_fiber_section()``

    """
    for dat in fib_sec_list:
        if dat[0] == 'section':
            secTag, GJ = dat[2], dat[4]
            ops.section('Fiber', secTag, '-GJ', GJ)

        if dat[0] == 'layer':
            matTag = dat[2]
            n_bars = dat[3]
            As = dat[4]
            if dat[1] == 'straight':
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                ops.layer('straight', matTag, n_bars, As, Iy, Iz, Jy, Jz)

            elif dat[1] == 'circ':
                yC, zC, radius = dat[5], dat[6], dat[7]

                if len(dat) > 8:
                    begAng, endAng = dat[8], dat[9]
                    ops.layer('circ', matTag, n_bars, As, yC, zC, radius,
                              begAng, endAng)
                else:
                    ops.layer('circ', matTag, n_bars, As, yC, zC, radius)

        if dat[0] == 'patch':
            matTag = dat[2]
            nIJ = dat[3]
            nJK = dat[4]

            if dat[1] == 'quad' or dat[1] == 'quadr':
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                Ky, Kz, Ly, Lz = dat[9], dat[10], dat[11], dat[12]
                ops.patch('quad', matTag, nIJ, nJK, Iy, Iz, Jy, Jz, Ky, Kz,
                          Ly, Lz)

            elif dat[1] == 'rect':
                Iy, Iz, Ky, Kz = dat[5], dat[6], dat[7], dat[8]
                Jy, Jz, Ly, Lz = Ky, Iz, Iy, Kz
                ops.patch('rect', matTag, nIJ, nJK, Iy, Iz, Ky, Kz)

            elif dat[1] == 'circ':
                yC, zC, intRad, extRad = dat[5], dat[6], dat[7], dat[8]

                begAng, endAng = dat[9], dat[10]
                ops.patch('circ', matTag, nIJ, nJK, yC, zC, intRad, extRad,
                          begAng, endAng)
