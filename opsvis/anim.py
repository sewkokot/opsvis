import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from settings import *
from defo import *


def _anim_mode_2d(modeNo, sfac, nep, unDefoFlag, fmt_defo, fmt_undefo,
                  interpFlag, endDispFlag, fig_wi_he,
                  fig_lbrt, xlim, ylim, ax):

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he
            fig, ax = plt.subplots(figsize=(fig_wi/2.54, fig_he/2.54))
        else:
            fig, ax = plt.subplots()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

    ele_tags = ops.getEleTags()

    nen = np.shape(ops.eleNodes(ele_tags[0]))[0]

    # truss and beam/frame elements
    if nen == 2:

        ndf = ops.getNDF()[0]

        # truss element
        if ndf == 2:

            for ele_tag in ele_tags:
                nd1, nd2 = ops.eleNodes(ele_tag)

                # element x, y coordinates
                ex = np.array([ops.nodeCoord(nd1)[0],
                               ops.nodeCoord(nd2)[0]])
                ey = np.array([ops.nodeCoord(nd1)[1],
                               ops.nodeCoord(nd2)[1]])

                if modeNo:
                    eux = np.array([ops.nodeEigenvector(nd1, modeNo)[0],
                                    ops.nodeEigenvector(nd2, modeNo)[0]])
                    euy = np.array([ops.nodeEigenvector(nd1, modeNo)[1],
                                    ops.nodeEigenvector(nd2, modeNo)[1]])
                else:
                    eux = np.array([ops.nodeDisp(nd1)[0],
                                    ops.nodeDisp(nd2)[0]])
                    euy = np.array([ops.nodeDisp(nd1)[1],
                                    ops.nodeDisp(nd2)[1]])

                # displaced element coordinates (scaled by sfac factor)
                edx = np.array([ex[0] + sfac*eux[0], ex[1] + sfac*eux[1]])
                edy = np.array([ey[0] + sfac*euy[0], ey[1] + sfac*euy[1]])

                if unDefoFlag:
                    ax.plot(ex, ey, **fmt_undefo)

                ax.plot(edx, edy, **fmt_defo)

        # beam/frame element anim eigen
        elif ndf == 3:

            # fig, ax = plt.subplots(figsize=(fig_wi/2.54, fig_he/2.54))
            ax.axis('equal')
            ax.set_xlim(xlim[0], xlim[1])
            ax.set_ylim(ylim[0], ylim[1])

            nel = len(ele_tags)
            Ex = np.zeros((nel, 2))
            Ey = np.zeros((nel, 2))
            Ed = np.zeros((nel, 6))
            # time vector for one cycle (period)
            n_cycles = 10
            n_frames = n_cycles * 32 + 1
            t = np.linspace(0., n_cycles*2*np.pi, n_frames)
            lines = []

            for i, ele_tag in enumerate(ele_tags):
                nd1, nd2 = ops.eleNodes(ele_tag)

                # element x, y coordinates
                Ex[i, :] = np.array([ops.nodeCoord(nd1)[0],
                                     ops.nodeCoord(nd2)[0]])
                Ey[i, :] = np.array([ops.nodeCoord(nd1)[1],
                                     ops.nodeCoord(nd2)[1]])

                Ed[i, :] = np.array([ops.nodeEigenvector(nd1, modeNo)[0],
                                     ops.nodeEigenvector(nd1, modeNo)[1],
                                     ops.nodeEigenvector(nd1, modeNo)[2],
                                     ops.nodeEigenvector(nd2, modeNo)[0],
                                     ops.nodeEigenvector(nd2, modeNo)[1],
                                     ops.nodeEigenvector(nd2, modeNo)[2]])

                lines.append(ax.plot([], [], **fmt_defo)[0])

            def init():
                for j, ele_tag in enumerate(ele_tags):
                    lines[j].set_data([], [])
                return lines

            def animate(i):
                for j, ele_tag in enumerate(ele_tags):

                    if interpFlag:
                        xcdi, ycdi = beam_defo_interp_2d(Ex[j, :],
                                                         Ey[j, :],
                                                         Ed[j, :],
                                                         sfac*np.sin(t[i]),
                                                         nep)
                        lines[j].set_data(xcdi, ycdi)
                    else:
                        xdi, ydi = beam_disp_ends(Ex[j, :], Ey[j, :], Ed[j, :],
                                                  sfac*np.cos(t[i]))
                        lines[j].set_data(xdi, ydi)

                    # plt.plot(xcdi, ycdi, **fmt_defo)

                return lines

            anim = FuncAnimation(fig, animate, init_func=init,
                                 frames=n_frames, interval=50, blit=True,
                                 repeat=False)

        # plt.axis('equal')
        # plt.show()  # call this from main py file for more control

    # 2d triangular elements - todo
    # elif nen == 3:
    #     x = ex+sfac*ed[:, [0, 2, 4]]
    #     y = ex+sfac*ed[:, [1, 3, 5]]
    #     xc = [x, x[0, :]]
    #     yc = [x, x[0, :]]

    # 2d quadrilateral (quad) elements
    elif nen == 4:
        for ele_tag in ele_tags:
            nd1, nd2, nd3, nd4 = ops.eleNodes(ele_tag)

            # element x, y coordinates
            ex = np.array([ops.nodeCoord(nd1)[0],
                           ops.nodeCoord(nd2)[0],
                           ops.nodeCoord(nd3)[0],
                           ops.nodeCoord(nd4)[0]])
            ey = np.array([ops.nodeCoord(nd1)[1],
                           ops.nodeCoord(nd2)[1],
                           ops.nodeCoord(nd3)[1],
                           ops.nodeCoord(nd4)[1]])

            if modeNo:
                ed = np.array([ops.nodeEigenvector(nd1, modeNo)[0],
                               ops.nodeEigenvector(nd1, modeNo)[1],
                               ops.nodeEigenvector(nd2, modeNo)[0],
                               ops.nodeEigenvector(nd2, modeNo)[1],
                               ops.nodeEigenvector(nd3, modeNo)[0],
                               ops.nodeEigenvector(nd3, modeNo)[1],
                               ops.nodeEigenvector(nd4, modeNo)[0],
                               ops.nodeEigenvector(nd4, modeNo)[1]])
            else:
                ed = np.array([ops.nodeDisp(nd1)[0],
                               ops.nodeDisp(nd1)[1],
                               ops.nodeDisp(nd2)[0],
                               ops.nodeDisp(nd2)[1],
                               ops.nodeDisp(nd3)[0],
                               ops.nodeDisp(nd3)[1],
                               ops.nodeDisp(nd4)[0],
                               ops.nodeDisp(nd4)[1]])

            if unDefoFlag:
                ax.plot(np.append(ex, ex[0]), np.append(ey, ey[0]),
                        **fmt_undefo)

            # xcdi, ycdi = beam_defo_interp_2d(ex, ey, ed, sfac, nep)
            # xdi, ydi = beam_disp_ends(ex, ey, ed, sfac)
            # # interpolated displacement field
            # plt.plot(xcdi, ycdi, 'b.-')
            # # translations of ends only
            # plt.plot(xdi, ydi, 'ro')

            # test it with one element
            x = ex+sfac*ed[[0, 2, 4, 6]]
            y = ey+sfac*ed[[1, 3, 5, 7]]
            ax.plot(np.append(x, x[0]), np.append(y, y[0]), 'b.-')

        ax.axis('equal')
        ax.grid(False)

    # 2d 8-node quadratic elements
    # elif nen == 8:
    #     x = ex+sfac*ed[:, [0, 2, 4, 6, 8, 10, 12, 14]]
    #     y = ex+sfac*ed[:, [1, 3, 5, 7, 9, 11, 13, 15]]

    #     t = -1
    #     n = 0
    #     for s in range(-1, 1.4, 0.4):
    #         n += 1
    #     ...

    else:
        print(f'\nWarning! Elements not supported yet. nen: {nen}; must be: 2, 3, 4, 8.')  # noqa: E501

    return anim


def anim_mode(modeNo, sfac=False, nep=17, unDefoFlag=1, fmt_defo=fmt_defo,
              fmt_undefo=fmt_undefo, interpFlag=1, endDispFlag=1,
              Eo=0, az_el=az_el,
              fig_wi_he=False, fig_lbrt=False, xlim=[0, 1], ylim=[0, 1],
              ax=False):
    """Make animation of a mode shape obtained from eigenvalue solution.

    Args:
        modeNo (int): indicates which mode shape to animate.

        sfac (float): scale factor

        nep (integer): number of evaluation points inside the element and
            including both element ends

        unDefoFlag (integer): 1 - plot the undeformed model (mesh), 0 - do not
            plot the mesh

        interpFlag (integer): 1 - interpolate deformation inside element,
            0 - no interpolation

        endDispFlag (integer): 1 - plot marks at element ends, 0 - no marks

        fmt_defo (dict): format line string for interpolated (continuous)
            deformated shape. The format contains information on line color,
            style and marks as in the standard matplotlib plot function.

        az_el (tuple): a tuple containing the azimuth and elevation

        fig_lbrt (tuple): a tuple contating left, bottom, right and top offsets

        fig_wi_he (tuple): contains width and height of the figure

    Returns:
        anim: animation object

    Notes:

    See also:
        anim_mode()
    """


    node_tags = ops.getNodeTags()

    # calculate sfac
    # min_x, min_y, min_z = np.inf, np.inf, np.inf
    # max_x, max_y, max_z = -np.inf, -np.inf, -np.inf
    # max_ux, max_uy, max_uz = -np.inf, -np.inf, -np.inf
    min_x, min_y = np.inf, np.inf
    max_x, max_y = -np.inf, -np.inf
    max_ux, max_uy = -np.inf, -np.inf
    ratio = 0.1

    ndim = ops.getNDM()[0]

    if ndim == 2:
        if not sfac:
            for node_tag in node_tags:
                x_crd = ops.nodeCoord(node_tag)[0]
                y_crd = ops.nodeCoord(node_tag)[1]
                ux = ops.nodeEigenvector(node_tag, modeNo)[0]
                uy = ops.nodeEigenvector(node_tag, modeNo)[1]

                min_x = min(min_x, x_crd)
                min_y = min(min_y, y_crd)
                max_x = max(max_x, x_crd)
                max_y = max(max_y, y_crd)
                max_ux = max(max_ux, np.abs(ux))
                max_uy = max(max_uy, np.abs(uy))

            dxmax = max_x - min_x
            dymax = max_y - min_y
            dlmax = max(dxmax, dymax)
            edmax = max(max_ux, max_uy)
            sfac = ratio * dlmax/edmax

        anim = _anim_mode_2d(modeNo, sfac, nep, unDefoFlag, fmt_defo,
                             fmt_undefo, interpFlag, endDispFlag,
                             fig_wi_he, fig_lbrt,
                             xlim, ylim, ax)

    # elif ndim == 3:
    #     if not sfac:
    #         for node_tag in node_tags:
    #             x_crd = ops.nodeCoord(node_tag)[0]
    #             y_crd = ops.nodeCoord(node_tag)[1]
    #             z_crd = ops.nodeCoord(node_tag)[2]
    #             ux = ops.nodeEigenvector(node_tag, modeNo)[0]
    #             uy = ops.nodeEigenvector(node_tag, modeNo)[1]
    #             uz = ops.nodeEigenvector(node_tag, modeNo)[2]

    #             min_x = min(min_x, x_crd)
    #             min_y = min(min_y, y_crd)
    #             min_z = min(min_z, z_crd)
    #             max_x = max(max_x, x_crd)
    #             max_y = max(max_y, y_crd)
    #             max_z = max(max_z, z_crd)
    #             max_ux = max(max_ux, np.abs(ux))
    #             max_uy = max(max_uy, np.abs(uy))
    #             max_uz = max(max_uz, np.abs(uz))

    #         dxmax = max_x - min_x
    #         dymax = max_y - min_y
    #         dzmax = max_z - min_z
    #         dlmax = max(dxmax, dymax, dzmax)
    #         edmax = max(max_ux, max_uy, max_uz)
    #         sfac = ratio * dlmax/edmax

    #     _plot_defo_mode_3d(modeNo, sfac, nep, unDefoFlag, fmt_defo,
    #                        fmt_undefo, interpFlag, endDispFlag,
    #                        Eo, az_el, fig_wi_he, fig_lbrt)

    else:
        print(f'\nWarning! ndim: {ndim} not supported yet.')

    return anim


def _anim_defo_2d(Eds, timeV, sfac, nep, unDefoFlag, fmt_defo, fmt_undefo,
                  interpFlag, endDispFlag,
                  fig_wi_he, fig_lbrt, xlim, ylim, ax):

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he
            fig, ax = plt.subplots(figsize=(fig_wi/2.54, fig_he/2.54))
        else:
            fig, ax = plt.subplots()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

    ele_tags = ops.getEleTags()

    nen = np.shape(ops.eleNodes(ele_tags[0]))[0]

    # truss and beam/frame elements
    if nen == 2:

        ndf = ops.getNDF()[0]

        # truss element
        if ndf == 2:

            for ele_tag in ele_tags:
                nd1, nd2 = ops.eleNodes(ele_tag)

                # element x, y coordinates
                ex = np.array([ops.nodeCoord(nd1)[0],
                               ops.nodeCoord(nd2)[0]])
                ey = np.array([ops.nodeCoord(nd1)[1],
                               ops.nodeCoord(nd2)[1]])

                eux = np.array([ops.nodeDisp(nd1)[0],
                                ops.nodeDisp(nd2)[0]])
                euy = np.array([ops.nodeDisp(nd1)[1],
                                ops.nodeDisp(nd2)[1]])

                # displaced element coordinates (scaled by sfac factor)
                edx = np.array([ex[0] + sfac*eux[0], ex[1] + sfac*eux[1]])
                edy = np.array([ey[0] + sfac*euy[0], ey[1] + sfac*euy[1]])

                if unDefoFlag:
                    ax.plot(ex, ey, **fmt_undefo)

                ax.plot(edx, edy, **fmt_defo)

        # beam/frame element anim defo
        elif ndf == 3:

            ax.axis('equal')
            ax.set_xlim(xlim[0], xlim[1])
            ax.set_ylim(ylim[0], ylim[1])
            # ax.grid()

            nel = len(ele_tags)
            Ex = np.zeros((nel, 2))
            Ey = np.zeros((nel, 2))
            # no of frames equal to time intervals
            n_frames, _, _ = np.shape(Eds)
            lines = []

            # time_text = ax.set_title('')  # does not work
            time_text = ax.text(.05, .95, '', transform=ax.transAxes)
            for i, ele_tag in enumerate(ele_tags):
                nd1, nd2 = ops.eleNodes(ele_tag)

                # element x, y coordinates
                Ex[i, :] = np.array([ops.nodeCoord(nd1)[0],
                                     ops.nodeCoord(nd2)[0]])
                Ey[i, :] = np.array([ops.nodeCoord(nd1)[1],
                                     ops.nodeCoord(nd2)[1]])

                lines.append(ax.plot([], [], **fmt_defo)[0])

            def init():
                for j, ele_tag in enumerate(ele_tags):
                    lines[j].set_data([], [])

                time_text.set_text('')

                return tuple(lines) + (time_text,)

            def animate(i):
                for j, ele_tag in enumerate(ele_tags):

                    if interpFlag:
                        xcdi, ycdi = beam_defo_interp_2d(Ex[j, :],
                                                         Ey[j, :],
                                                         Eds[i, j, :],
                                                         sfac,
                                                         nep)
                        lines[j].set_data(xcdi, ycdi)
                    else:
                        xdi, ydi = beam_disp_ends(Ex[j, :], Ey[j, :],
                                                  Eds[i, j, :], sfac)
                        lines[j].set_data(xdi, ydi)

                    # plt.plot(xcdi, ycdi, **fmt_defo)

                # time_text.set_text(f'f')
                time_text.set_text(f'frame: {i+1}/{n_frames}, \
time: {timeV[i]:.3f} s')

                return tuple(lines) + (time_text,)

            anim = FuncAnimation(fig, animate, init_func=init, frames=n_frames,
                                 interval=50, blit=True, repeat=False)

        # plt.axis('equal')
        # plt.show()  # call this from main py file for more control

    # 2d triangular elements
    # elif nen == 3:
    #     x = ex+sfac*ed[:, [0, 2, 4]]
    #     y = ex+sfac*ed[:, [1, 3, 5]]
    #     xc = [x, x[0, :]]
    #     yc = [x, x[0, :]]

    # 2d quadrilateral (quad) elements
    elif nen == 4:
        for ele_tag in ele_tags:
            nd1, nd2, nd3, nd4 = ops.eleNodes(ele_tag)

            # element x, y coordinates
            ex = np.array([ops.nodeCoord(nd1)[0],
                           ops.nodeCoord(nd2)[0],
                           ops.nodeCoord(nd3)[0],
                           ops.nodeCoord(nd4)[0]])
            ey = np.array([ops.nodeCoord(nd1)[1],
                           ops.nodeCoord(nd2)[1],
                           ops.nodeCoord(nd3)[1],
                           ops.nodeCoord(nd4)[1]])

            # if modeNo:
            #     ed = np.array([ops.nodeEigenvector(nd1, modeNo)[0],
            #                    ops.nodeEigenvector(nd1, modeNo)[1],
            #                    ops.nodeEigenvector(nd2, modeNo)[0],
            #                    ops.nodeEigenvector(nd2, modeNo)[1],
            #                    ops.nodeEigenvector(nd3, modeNo)[0],
            #                    ops.nodeEigenvector(nd3, modeNo)[1],
            #                    ops.nodeEigenvector(nd4, modeNo)[0],
            #                    ops.nodeEigenvector(nd4, modeNo)[1]])
            # else:
            ed = np.array([ops.nodeDisp(nd1)[0],
                           ops.nodeDisp(nd1)[1],
                           ops.nodeDisp(nd2)[0],
                           ops.nodeDisp(nd2)[1],
                           ops.nodeDisp(nd3)[0],
                           ops.nodeDisp(nd3)[1],
                           ops.nodeDisp(nd4)[0],
                           ops.nodeDisp(nd4)[1]])

            if unDefoFlag:
                ax.plot(np.append(ex, ex[0]), np.append(ey, ey[0]),
                        **fmt_undefo)

            # xcdi, ycdi = beam_defo_interp_2d(ex, ey, ed, sfac, nep)
            # xdi, ydi = beam_disp_ends(ex, ey, ed, sfac)
            # # interpolated displacement field
            # plt.plot(xcdi, ycdi, 'b.-')
            # # translations of ends only
            # plt.plot(xdi, ydi, 'ro')

            # test it with one element
            x = ex+sfac*ed[[0, 2, 4, 6]]
            y = ey+sfac*ed[[1, 3, 5, 7]]
            ax.plot(np.append(x, x[0]), np.append(y, y[0]), 'b.-')

        ax.axis('equal')
        ax.grid(False)

    # 2d 8-node quadratic elements
    # elif nen == 8:
    #     x = ex+sfac*ed[:, [0, 2, 4, 6, 8, 10, 12, 14]]
    #     y = ex+sfac*ed[:, [1, 3, 5, 7, 9, 11, 13, 15]]

    #     t = -1
    #     n = 0
    #     for s in range(-1, 1.4, 0.4):
    #         n += 1
    #     ...

    else:
        print(f'\nWarning! Elements not supported yet. nen: {nen}; must be: 2, 3, 4, 8.')  # noqa: E501

    return anim


def anim_defo(Eds, timeV, sfac, nep=17, unDefoFlag=1, fmt_defo=fmt_defo,
              fmt_undefo=fmt_undefo, interpFlag=1, endDispFlag=1,
              az_el=az_el,
              fig_wi_he=False, fig_lbrt=False, xlim=[0, 1],
              ylim=[0, 1], ax=False):
    """Make animation of the deformed shape computed by transient analysis

    Args:
        Eds (ndarray): A 3d array (n_steps x n_eles x n_dof_per_element)
            containing the collected displacements per element for all
            time steps.

        timeV (1darray): vector of discretized time values

        sfac (float): scale factor

        nep (integer): number of evaluation points inside the element and
            including both element ends

        unDefoFlag (integer): 1 - plot the undeformed model (mesh), 0 - do not
            plot the mesh

        interpFlag (integer): 1 - interpolate deformation inside element,
            0 - no interpolation

        endDispFlag (integer): 1 - plot marks at element ends, 0 - no marks

        fmt_defo (dict): format line string for interpolated (continuous)
            deformated shape. The format contains information on line color,
            style and marks as in the standard matplotlib plot function.

        az_el (tuple): a tuple containing the azimuth and elevation

        fig_lbrt (tuple): a tuple contating left, bottom, right and top offsets

        fig_wi_he (tuple): contains width and height of the figure

    Examples:

    Notes:

    See also:
    """

    node_tags = ops.getNodeTags()

    ndim = ops.getNDM()[0]

    if ndim == 2:
        anim = _anim_defo_2d(Eds, timeV, sfac, nep, unDefoFlag, fmt_defo,
                             fmt_undefo, interpFlag, endDispFlag,
                             fig_wi_he, fig_lbrt, xlim, ylim, ax)

    else:
        print(f'\nWarning! ndim: {ndim} not supported yet.')

    return anim
