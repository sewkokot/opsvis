import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Polygon, Wedge
from matplotlib.animation import FuncAnimation
import matplotlib.tri as tri

from .settings import *
from . import model as opsvmodel


def _plot_defo_mode_2d(modeNo, sfac, nep, unDefoFlag, fmt_defo, fmt_undefo,
                       fmt_defo_faces, fmt_undefo_faces,
                       interpFlag, endDispFlag, fmt_nodes,
                       fig_wi_he, fig_lbrt, node_supports, ax):

    node_tags = ops.getNodeTags()
    ele_tags = ops.getEleTags()

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he
            fig, ax = plt.subplots(figsize=(fig_wi/2.54, fig_he/2.54))
        else:
            fig, ax = plt.subplots()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]
        nen = np.shape(ops.eleNodes(ele_tag))[0]

        if (ele_classtag == EleClassTag.truss
            or ele_classtag == EleClassTag.trussSection):

            nen, ndf = 2, 2
            ele_node_tags = ops.eleNodes(ele_tag)
            # nd1, nd2 = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)[:2]
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)[:2]

            # displaced element coordinates (scaled by sfac factor)
            xy = ecrd + sfac * ed

            if unDefoFlag:
                ax.plot(ecrd[:, 0], ecrd[:, 1], **fmt_undefo)

            ax.plot(xy[:, 0], xy[:, 1], **fmt_defo)

        elif (ele_classtag == EleClassTag.ZeroLength
              or ele_classtag == EleClassTag.ZeroLengthSection
              or ele_classtag == EleClassTag.CoupledZeroLength
              or ele_classtag == EleClassTag.TwoNodeLink):

            # nen, ndf = 2, 3
            ndf = ops.getNDF()[0]
            ele_node_tags = ops.eleNodes(ele_tag)
            nen = len(ele_node_tags)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)


            # fix for a sdof system, or specify sfac explicitly
            if np.isclose(sfac, 0.0):
                sfac = 1.0

            # displaced element coordinates (scaled by sfac factor)
            xy = ecrd + sfac * ed[:, :2]

            if unDefoFlag:
                ax.plot(ecrd[:, 0], ecrd[:, 1], **fmt_undefo)

            ax.plot(xy[:, 0], xy[:, 1], **fmt_defo_zeroLenght)

        # beam/frame element plot_defo
        elif (ele_classtag == EleClassTag.ElasticBeam2d
              or ele_classtag == EleClassTag.ForceBeamColumn2d
              or ele_classtag == EleClassTag.DispBeamColumn2d
              or ele_classtag == EleClassTag.TimoshenkoBeamColumn2d
              or ele_classtag == EleClassTag.ElasticTimoshenkoBeam2d):

            nen, ndf = 2, 3

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd_nodes = np.zeros((nen, 2))
            # ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                # ecrd[i, :] = ops.nodeCoord(ele_node_tag)
                ecrd_nodes[i, :] = ops.nodeCoord(ele_node_tag)

            ecrd_eles0 = np.copy(ecrd_nodes)
            ecrd_eles = np.copy(ecrd_nodes)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            ele_offsets = np.array(ops.eleResponse(ele_tag, 'offsets'))
            nz_offsets = np.nonzero(ele_offsets)[0]  # tuple of arrays
            if np.any(nz_offsets):
                Lro_l, fi0_l = ro_length_init_rot(ele_offsets[0], ele_offsets[1])
                Lro_r, fi0_r = ro_length_init_rot(ele_offsets[3], ele_offsets[4])
                # ed[0, 2] - rotations
                x_ro_rot_l, y_ro_rot_l = ro_rotated(Lro_l, fi0_l, ed[0, 2]*sfac)
                x_ro_rot_r, y_ro_rot_r = ro_rotated(Lro_r, fi0_r, ed[1, 2]*sfac)
                # modify ecrd_eles
                ecrd_eles0[:, 0] += ele_offsets[[0, 3]]
                ecrd_eles0[:, 1] += ele_offsets[[1, 4]]
                ecrd_eles[:, 0] += [x_ro_rot_l, x_ro_rot_r]
                ecrd_eles[:, 1] += [y_ro_rot_l, y_ro_rot_r]
            else:
                pass

            if unDefoFlag:
                ax.plot(ecrd_eles0[:, 0], ecrd_eles0[:, 1], **fmt_undefo)

            # interpolated displacement field
            if interpFlag:
                xcdi, ycdi = beam_defo_interp_2d(ecrd_eles, ed.flatten(), sfac, nep)
            else:
                xcdi, ycdi = beam_disp_ends(ecrd_eles, ed.flatten(), sfac)

            ax.plot(xcdi, ycdi, **fmt_defo)

            # plot rigid offsets
            # xdi, ydi = beam_disp_ro(ecrd_eles[:, 0],
            #                         ecrd_eles[:, 1],
            #                         ed.flatten(), sfac)
            ax.plot([ecrd_nodes[0, 0] + sfac * ed[0, 0],
                     ecrd_eles[0, 0] + sfac * ed[0, 0]] ,
                    [ecrd_nodes[0, 1] + sfac * ed[0, 1],
                     ecrd_eles[0, 1] + sfac * ed[0, 1]] ,
                    **fmt_model_rigid_offset)
            ax.plot([ecrd_nodes[1, 0] + sfac * ed[1, 0],
                     ecrd_eles[1, 0] + sfac * ed[1, 0]] ,
                    [ecrd_nodes[1, 1] + sfac * ed[1, 1],
                     ecrd_eles[1, 1] + sfac * ed[1, 1]],
                    **fmt_model_rigid_offset)

            # ax.plot([ecrd_nodes[1, 0], ecrd_eles[1, 0]],
            #         [ecrd_nodes[1, 1], ecrd_eles[1, 1]],
            #         **fmt_model_rigid_offset)

            # ax.plot([ecrd[0, 0] + sfac*ed[0], ecrd[1, 0] + sfac*ed[3]],
            #         [ecrd[0, 1] + sfac*ed[1], ecrd[1, 1] + sfac*ed[4]], **fmt_model_rigid_offset)

            # translations of ends
            if endDispFlag:
                xdi, ydi = beam_disp_ends(ecrd_eles, ed.flatten(), sfac)
                ax.plot(xdi, ydi, **fmt_nodes)

        # 2d triangular (tri31) elements plot_defo
        elif (ele_classtag == EleClassTag.tri3n):

            nen, ndf = 3, 2
            nodes_geo_order = [0, 1, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        **fmt_undefo)

            xy = ecrd + sfac * ed

            ax.plot(xy[nodes_geo_order, 0],
                    xy[nodes_geo_order, 1],
                    **fmt_defo)

        # 2d quadrilateral (quad) elements plot_defo
        elif (ele_classtag == EleClassTag.quad4n):

            nen, ndf = 4, 2
            nodes_geo_order = [0, 1, 2, 3, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        **fmt_undefo)
                # verts = [ecrd[nodes_geo_order]]
                # ax.add_collection(PolyCollection(verts, **fmt_undefo_faces))

            xy = ecrd + sfac * ed

            verts = [xy[nodes_geo_order]]
            ax.add_collection(PolyCollection(verts, **fmt_defo_faces))

        # 2d quadrilateral (quad8n) elements plot_defo
        elif (ele_classtag == EleClassTag.quad8n):

            nen, ndf = 8, 2
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 3, 7, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        **fmt_undefo)

            xy = ecrd + sfac * ed

            ax.plot(xy[nodes_geo_order, 0],
                    xy[nodes_geo_order, 1],
                    **fmt_defo)

        # 2d quadrilateral (quad9n) elements plot_defo
        elif (ele_classtag == EleClassTag.quad9n):

            nen, ndf = 9, 2
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 3, 7, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        **fmt_undefo)

            xy = ecrd + sfac * ed

            ax.plot(xy[nodes_geo_order, 0],
                    xy[nodes_geo_order, 1],
                    **fmt_defo)
            ax.plot([xy[8, 0]], [xy[8, 1]], 'b.-')

        # 2d triangle (tri6n) elements plot_defo
        elif (ele_classtag == EleClassTag.tri6n):

            nen, ndf = 6, 2
            nodes_geo_order = [0, 3, 1, 4, 2, 5, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        **fmt_undefo)

            xy = ecrd + sfac * ed

            ax.plot(xy[nodes_geo_order, 0],
                    xy[nodes_geo_order, 1],
                    **fmt_defo)

        # 2d joint2d element plot_defo
        elif (ele_classtag == EleClassTag.Joint2D):

            nen, ndf = 5, 3

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags[:4]):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags[:4]):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags[:4]):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            xy = ecrd + sfac * ed[:, :2]

            ax.plot([xy[0, 0], xy[2, 0], xy[2, 0], xy[0, 0],
                     xy[0, 0]],
                    [xy[1, 1], xy[1, 1], xy[3, 1], xy[3, 1],
                     xy[1, 1]], **fmt_model_joint2d)

        else:
            print(f'\nWarning! Elements not supported yet. nen: {nen}; must be: 2, 3, 4, 8.')  # noqa: E501

    if node_supports:
        opsvmodel._plot_supports(node_tags, ax)

    ax.axis('equal')
    ax.grid(False)

    return ax


def ro_length_init_rot(xo, yo):
    """Return lenght and initial rotation of the rigid offset in 2d
    """

    L = np.sqrt(xo**2. + yo**2.)
    # fi0 = np.arctan(yo/xo)
    fi0 = np.arctan2(yo, xo)

    return L, fi0


def ro_length_init_rot_3d(xo, yo, zo, along='x'):
    """Return lenght and initial rotation of the rigid offset in 3d
    """
    L = np.sqrt(xo**2. + yo**2. + zo**2.)
    if along == 'x':
        fi0 = np.arctan2(zo, xo)
    elif along == 'y':
        fi0 = np.arctan2(zo, xo)
    elif along == 'z':
        fi0 = 0.


    return L, fi0


def ro_rotated(L, fi0, fi):
    dfi = fi0 + fi
    x = L * np.cos(dfi)
    y = L * np.sin(dfi)
    return x, y


def _plot_defo_mode_3d(modeNo, sfac, nep, unDefoFlag, fmt_defo, fmt_undefo,
                       fmt_defo_faces, fmt_undefo_faces,
                       interpFlag, endDispFlag, fmt_nodes, az_el,
                       fig_wi_he, fig_lbrt, node_supports, ax):

    node_tags = ops.getNodeTags()
    ele_tags = ops.getEleTags()

    azim, elev = az_el

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he

            fig = plt.figure(figsize=(fig_wi/2.54, fig_he/2.54))
            # fig.subplots_adjust(left=.08, bottom=.08, right=.985, top=.94)

        else:
            fig = plt.figure()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

        ax = fig.add_subplot(111, projection=Axes3D.name)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.view_init(azim=azim, elev=elev)

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]
        nen = np.shape(ops.eleNodes(ele_tag))[0]

        # if (ele_classtag == EleClassTag.truss):
        if (ele_classtag == EleClassTag.ElasticBeam3d
            or ele_classtag == EleClassTag.ForceBeamColumn3d
            or ele_classtag == EleClassTag.DispBeamColumn3d
            or ele_classtag == EleClassTag.ElasticTimoshenkoBeam3d):

            nen, ndf = 2, 6

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            # eo = Eo[i, :]
            g = np.vstack((ops.eleResponse(ele_tag, 'xlocal'),
                           ops.eleResponse(ele_tag, 'ylocal'),
                           ops.eleResponse(ele_tag, 'zlocal')))

            if unDefoFlag:
                plt.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], **fmt_undefo)

            if interpFlag:
                # xcd, ycd, zcd = beam_defo_interp_3d(ex, ey, ez, g,
                #                                     ed, sfac, nep)
                xcd, ycd, zcd = beam_defo_interp_3d(ecrd, g,
                                                    ed.flatten(), sfac, nep)
            else:
                xcd, ycd, zcd = beam_disp_ends3d(ecrd, ed.flatten(), sfac)

            ax.plot(xcd, ycd, zcd, **fmt_defo)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            # translations of ends
            if endDispFlag:
                xd, yd, zd = beam_disp_ends3d(ecrd, ed.flatten(), sfac)
                ax.plot(xd, yd, zd, **fmt_nodes)

        elif (ele_classtag == EleClassTag.truss
              or ele_classtag == EleClassTag.trussSection):

            nen, ndf = 2, 3

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)[:3]
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)[:3]

            if unDefoFlag:
                plt.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], **fmt_undefo)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            # displaced element coordinates (scaled by sfac factor)
            xyz = ecrd + sfac * ed
            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], **fmt_defo)

        elif (ele_classtag == EleClassTag.ZeroLength
              or ele_classtag == EleClassTag.ZeroLengthSection
              or ele_classtag == EleClassTag.CoupledZeroLength
              or ele_classtag == EleClassTag.TwoNodeLink):

            # nen, ndf = 2, 6
            ndf = ops.getNDF()[0]
            ele_node_tags = ops.eleNodes(ele_tag)
            nen = len(ele_node_tags)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                plt.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], **fmt_undefo)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            # fix for a sdof system, or specify sfac explicitly
            if np.isclose(sfac, 0.0):
                sfac = 1.0

            # displaced element coordinates (scaled by sfac factor)
            xyz = ecrd + sfac * ed[:, :3]

            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], **fmt_defo_zeroLenght)

        # plot: four node shell in 3d
        elif ele_classtag in [EleClassTag.ShellMITC4,
                               EleClassTag.ASDShellQ4,
                               EleClassTag.ShellDKGQ,EleClassTag.ShellNLDKGQ]:

            # ndf is limited to translations x, y, z (no rotations)
            nen, ndf = 4, 3
            nodes_geo_order = [0, 1, 2, 3, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)[:3]
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)[:3]

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        ecrd[nodes_geo_order, 2],
                        **fmt_undefo)

            xyz = ecrd + sfac * ed
            # x = ex+sfac*ed[[0, 3, 6, 9]]
            # y = ey+sfac*ed[[1, 4, 7, 10]]
            # z = ez+sfac*ed[[2, 5, 8, 11]]

            verts = [[xyz[0], xyz[1], xyz[2], xyz[3]]]
            ax.add_collection3d(Poly3DCollection(verts, **fmt_defo_faces))

            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=0)

        # plot: three node shell in 3d
        elif ele_classtag in [EleClassTag.ShellDKGT,EleClassTag.ShellNLDKGT]:

            # ndf is limited to translations x, y, z (no rotations)
            nen, ndf = 3, 3
            nodes_geo_order = [0, 1, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))  # without rotations

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)[:3]
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)[:3]

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        ecrd[nodes_geo_order, 2],
                        **fmt_undefo)

            xyz = ecrd + sfac * ed
            # x = ex+sfac*ed[[0, 3, 6, 9]]
            # y = ey+sfac*ed[[1, 4, 7, 10]]
            # z = ez+sfac*ed[[2, 5, 8, 11]]

            verts = [[xyz[0], xyz[1], xyz[2]]]
            ax.add_collection3d(Poly3DCollection(verts, **fmt_defo_faces))

            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=0)

        # 4-node tetrahedron, 3d defo
        elif (ele_classtag == EleClassTag.FourNodeTetrahedron):

            nen, ndf = 4, 3
            nodes_geo_order = [0, 1, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        ecrd[nodes_geo_order, 2],
                        **fmt_undefo)

                for j in range(3):
                    ax.plot(ecrd[[j, 3], 0],
                            ecrd[[j, 3], 1],
                            ecrd[[j, 3], 2], **fmt_undefo)

            xyz = ecrd + sfac * ed

            ax.plot(xyz[nodes_geo_order, 0],
                    xyz[nodes_geo_order, 1],
                    xyz[nodes_geo_order, 2],
                    'b.-', lw=0.4, ms=2, mfc='g', mec='g')

            for j in range(3):
                ax.plot(xyz[[j, 3], 0],
                         xyz[[j, 3], 1],
                         xyz[[j, 3], 2], 'b.-',
                         lw=0.4, ms=2, mfc='g', mec='g')

        # tetra10n, 3d defo
        elif (ele_classtag in {EleClassTag.TenNodeTetrahedron,
                               EleClassTag.TenNodeTetrahedronSK}):

            nen, ndf = 10, 3
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        ecrd[nodes_geo_order, 2],
                        **fmt_undefo)

                for j in range(3):
                    ax.plot(ecrd[[j, 7+j, 3], 0],
                            ecrd[[j, 7+j, 3], 1],
                            ecrd[[j, 7+j, 3], 2], **fmt_undefo)

            xyz = ecrd + sfac * ed

            ax.plot(xyz[nodes_geo_order, 0],
                    xyz[nodes_geo_order, 1],
                    xyz[nodes_geo_order, 2],
                    'b.-', lw=0.4, ms=2, mfc='g', mec='g')

            for j in range(3):
                ax.plot(xyz[[j, 7+j, 3], 0],
                        xyz[[j, 7+j, 3], 1],
                        xyz[[j, 7+j, 3], 2], 'b.-',
                        lw=0.4, ms=2, mfc='g', mec='g')

        # 8-node brick, 3d defo
        elif (ele_classtag == EleClassTag.brick8n):

            nen, ndf = 8, 3
            nodes_geo_order = np.array([0, 1, 2, 3, 0], dtype=int)  # bottom face nodes loop

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                ax.plot(ecrd[nodes_geo_order, 0],
                        ecrd[nodes_geo_order, 1],
                        ecrd[nodes_geo_order, 2],
                        **fmt_undefo)
                ax.plot(ecrd[nodes_geo_order + 4, 0],
                        ecrd[nodes_geo_order + 4, 1],
                        ecrd[nodes_geo_order + 4, 2],
                        **fmt_undefo)

                for j in range(4):
                    ax.plot(ecrd[[j, j+4], 0],
                            ecrd[[j, j+4], 1],
                            ecrd[[j, j+4], 2], **fmt_undefo)

            xyz = ecrd + sfac * ed

            ax.plot(xyz[nodes_geo_order, 0],
                    xyz[nodes_geo_order, 1],
                    xyz[nodes_geo_order, 2],
                    **fmt_defo)
            ax.plot(xyz[nodes_geo_order + 4, 0],
                    xyz[nodes_geo_order + 4, 1],
                    xyz[nodes_geo_order + 4, 2],
                    **fmt_defo)

            for j in range(4):
                ax.plot(xyz[[j, j+4], 0],
                        xyz[[j, j+4], 1],
                        xyz[[j, j+4], 2], **fmt_defo)

        # plot: quad in 3d defo, but also shell which case more dofs
        # elif (ele_classtag == EleClassTag.quad4n):

        # 20-node brick20n, 3d defo
        elif (ele_classtag == EleClassTag.brick20n):

            nen, ndf = 20, 3
            nodes_geo_order = [0, 1, 2, 3, 0]  # bottom face nodes loop
            nodes_geo_order2 = [4, 5, 6, 7, 4]  # top face nodes loop

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            if modeNo:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeEigenvector(ele_node_tag, modeNo)
            else:
                for i, ele_node_tag in enumerate(ele_node_tags):
                    ed[i, :] = ops.nodeDisp(ele_node_tag)

            if unDefoFlag:
                for j in range(2):
                    ax.plot(ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 0],
                            ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 1],
                            ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 2],
                            **fmt_undefo)

                for j in range(4):
                    ax.plot(ecrd[[j, 16+j, 4+j], 0],
                            ecrd[[j, 16+j, 4+j], 1],
                            ecrd[[j, 16+j, 4+j], 2], **fmt_undefo)

            xyz = ecrd + sfac * ed

            for j in range(2):
                ax.plot(xyz[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 0],
                        xyz[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 1],
                        xyz[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 2],
                        'b.-', lw=0.4, ms=2, mfc='g', mec='g')

            for j in range(4):
                ax.plot(xyz[[j, 16+j, 4+j], 0],
                        xyz[[j, 16+j, 4+j], 1],
                        xyz[[j, 16+j, 4+j], 2], 'b.-',
                        lw=0.4, ms=2, mfc='g', mec='g')

    if node_supports:
        opsvmodel._plot_supports(node_tags, ax)

    ax.set_box_aspect((np.ptp(ax.get_xlim3d()),
                       np.ptp(ax.get_ylim3d()),
                       np.ptp(ax.get_zlim3d())))

    return ax


def plot_defo(sfac=False, nep=17, unDefoFlag=1, fmt_defo=fmt_defo,
              fmt_undefo=fmt_undefo, fmt_defo_faces=fmt_defo_faces,
              fmt_undefo_faces=fmt_undefo_faces, interpFlag=1, endDispFlag=0,
              fmt_nodes=fmt_nodes, Eo=0, az_el=az_el,
              fig_wi_he=False, fig_lbrt=False, node_supports=True, ax=False):
    """Plot deformed shape of the structure.

    Args:
        sfac (float): scale factor to increase/decrease displacements obtained
            from FE analysis. If not specified (False), sfac is automatically
            calculated based on the maximum overall displacement and this
            maximum displacement is plotted as 20 percent (hardcoded) of
            the maximum model dimension.

        interpFlag (int): 1 - use interpolated deformation using shape
            function, 0 - do not use interpolation, just show displaced element
            nodes (default is 1)

        nep (int): number of evaluation points for shape function interpolation
            (default: 17)

        node_supports (bool): True - show the supports. Default: True.

    Returns:
        sfac (float): the automatically calculated scale factor can be
            returned.

    Usage:

    ``sfac = plot_defo()`` - plot deformed shape with default parameters and
    automatically calcutated scale factor.

    ``plot_defo(interpFlag=0)`` - plot simplified deformation by
    displacing the nodes connected with straight lines (shape function
    interpolation)

    ``plot_defo(sfac=1.5)`` - plot with specified scale factor

    ``plot_defo(unDefoFlag=0, endDispFlag=0)`` - plot without showing
    undeformed (original) mesh and without showing markers at the
    element ends.
    """
    node_tags = ops.getNodeTags()

    # calculate sfac
    min_x, min_y, min_z = np.inf, np.inf, np.inf
    max_x, max_y, max_z = -np.inf, -np.inf, -np.inf
    max_ux, max_uy, max_uz = -np.inf, -np.inf, -np.inf
    ratio = 0.1

    ndim = ops.getNDM()[0]

    if ndim == 2:
        if not sfac:
            for node_tag in node_tags:
                x_crd = ops.nodeCoord(node_tag)[0]
                y_crd = ops.nodeCoord(node_tag)[1]
                ux = ops.nodeDisp(node_tag)[0]
                uy = ops.nodeDisp(node_tag)[1]

                min_x = min(min_x, x_crd)
                min_y = min(min_y, y_crd)
                max_x = max(max_x, x_crd)
                max_y = max(max_y, y_crd)
                max_ux = max(max_ux, np.abs(ux))
                max_uy = max(max_uy, np.abs(uy))

            dxmax = max_x - min_x
            dymax = max_y - min_y
            dlmax = max(dxmax, dymax)
            max_u_abs = max(max_ux, max_uy)

            max_u_abs = max_u_abs_from_beam_defo_interp_2d(max_u_abs, nep)

            sfac = ratio * dlmax/max_u_abs

        ax = _plot_defo_mode_2d(0, sfac, nep, unDefoFlag, fmt_defo,
                                fmt_undefo, fmt_defo_faces,
                                fmt_undefo_faces, interpFlag,
                                endDispFlag, fmt_nodes, fig_wi_he,
                                fig_lbrt, node_supports, ax)

    elif ndim == 3:
        if not sfac:
            for node_tag in node_tags:
                x_crd = ops.nodeCoord(node_tag)[0]
                y_crd = ops.nodeCoord(node_tag)[1]
                z_crd = ops.nodeCoord(node_tag)[2]
                ux = ops.nodeDisp(node_tag)[0]
                uy = ops.nodeDisp(node_tag)[1]
                uz = ops.nodeDisp(node_tag)[2]

                min_x = min(min_x, x_crd)
                min_y = min(min_y, y_crd)
                min_z = min(min_z, z_crd)
                max_x = max(max_x, x_crd)
                max_y = max(max_y, y_crd)
                max_z = max(max_z, z_crd)
                max_ux = max(max_ux, np.abs(ux))
                max_uy = max(max_uy, np.abs(uy))
                max_uz = max(max_uz, np.abs(uz))

            dxmax = max_x - min_x
            dymax = max_y - min_y
            dzmax = max_z - min_z
            dlmax = max(dxmax, dymax, dzmax)
            max_u_abs = max(max_ux, max_uy, max_uz)

            max_u_abs = max_u_abs_from_beam_defo_interp_3d(max_u_abs, nep)

            sfac = ratio * dlmax/max_u_abs

        ax = _plot_defo_mode_3d(0, sfac, nep, unDefoFlag, fmt_defo,
                                fmt_undefo, fmt_defo_faces,
                                fmt_undefo_faces, interpFlag,
                                endDispFlag, fmt_nodes, az_el,
                                fig_wi_he, fig_lbrt, node_supports,
                                ax)

    else:
        print(f'\nWarning! ndim: {ndim} not supported yet.')

    return sfac


def plot_mode_shape(modeNo, sfac=False, nep=17, unDefoFlag=1,
                    fmt_defo=fmt_defo,
                    fmt_undefo=fmt_undefo, fmt_defo_faces=fmt_defo_faces,
                    fmt_undefo_faces=fmt_undefo_faces,
                    interpFlag=1, endDispFlag=1,
                    fmt_nodes=fmt_nodes, Eo=0,
                    az_el=az_el, fig_wi_he=False,
                    fig_lbrt=False, node_supports=True, ax=False):
    """Plot mode shape of the structure obtained from eigenvalue analysis.

    Args:
        modeNo (int): indicates which mode shape to plot

        sfac (float): scale factor to increase/decrease displacements obtained
            from FE analysis. If not specified (False), sfac is automatically
            calculated based on the maximum overall displacement and this
            maximum displacement is plotted as 20 percent (hardcoded) of
            the maximum model dimension.

        interpFlag (int): 1 - use interpolated deformation using shape
            function, 0 - do not use interpolation, just show displaced element
            nodes (default is 1)

        nep (int): number of evaluation points for shape function interpolation
            (default: 17)

    Usage:

    ``plot_mode_shape(1)`` - plot the first mode shape with default parameters
    and automatically calcutated scale factor.

    ``plot_mode_shape(2, interpFlag=0)`` - plot the 2nd mode shape by
    displacing the nodes connected with straight lines (shape function
    interpolation)

    ``plot_mode_shape(3, sfac=1.5)`` - plot the 3rd mode shape with specified
    scale factor

    ``plot_mode_shape(4, unDefoFlag=0, endDispFlag=0)`` - plot the 4th mode
    shape without showing undeformed (original) mesh and without showing
    markers at the element ends.

    Examples:

    Notes:

    See also:
        plot_defo()
    """

    node_tags = ops.getNodeTags()

    # calculate sfac
    min_x, min_y, min_z = np.inf, np.inf, np.inf
    max_x, max_y, max_z = -np.inf, -np.inf, -np.inf
    max_ux, max_uy, max_uz = -np.inf, -np.inf, -np.inf
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
            max_u_abs = max(max_ux, max_uy)
            sfac = ratio * dlmax/max_u_abs

        ax = _plot_defo_mode_2d(modeNo, sfac, nep, unDefoFlag,
                                fmt_defo, fmt_undefo, fmt_defo_faces,
                                fmt_undefo_faces, interpFlag,
                                endDispFlag, fmt_nodes, fig_wi_he,
                                fig_lbrt, node_supports, ax)

    elif ndim == 3:
        if not sfac:
            for node_tag in node_tags:
                x_crd = ops.nodeCoord(node_tag)[0]
                y_crd = ops.nodeCoord(node_tag)[1]
                z_crd = ops.nodeCoord(node_tag)[2]
                ux = ops.nodeEigenvector(node_tag, modeNo)[0]
                uy = ops.nodeEigenvector(node_tag, modeNo)[1]
                uz = ops.nodeEigenvector(node_tag, modeNo)[2]

                min_x = min(min_x, x_crd)
                min_y = min(min_y, y_crd)
                min_z = min(min_z, z_crd)
                max_x = max(max_x, x_crd)
                max_y = max(max_y, y_crd)
                max_z = max(max_z, z_crd)
                max_ux = max(max_ux, np.abs(ux))
                max_uy = max(max_uy, np.abs(uy))
                max_uz = max(max_uz, np.abs(uz))

            dxmax = max_x - min_x
            dymax = max_y - min_y
            dzmax = max_z - min_z
            dlmax = max(dxmax, dymax, dzmax)
            max_u_abs = max(max_ux, max_uy, max_uz)
            sfac = ratio * dlmax/max_u_abs

        ax = _plot_defo_mode_3d(modeNo, sfac, nep, unDefoFlag,
                                fmt_defo, fmt_undefo, fmt_defo_faces,
                                fmt_undefo_faces, interpFlag,
                                endDispFlag, fmt_nodes, az_el,
                                fig_wi_he, fig_lbrt, node_supports,
                                ax)

    else:
        print(f'\nWarning! ndim: {ndim} not supported yet.')


def beam_disp_ends(ecrd, d, sfac):
    """
    Calculate the element deformation at element ends only.
    """

    #  indx: 0   1   2   3   4   5
    # Ed = ux1 uy1 ur1 ux2 uy2 ur2
    exd = np.array([ecrd[0, 0] + sfac * d[0], ecrd[1, 0] + sfac * d[3]])
    eyd = np.array([ecrd[0, 1] + sfac * d[1], ecrd[1, 1] + sfac * d[4]])

    return exd, eyd


def beam_disp_ro(ecrd, d, sfac):
    """
    Calculate the element deformation at element ends only.
    """

    #  indx: 0   1   2   3   4   5
    # Ed = ux1 uy1 ur1 ux2 uy2 ur2
    exd = np.array([ecrd[0, 0] + sfac * d[0], ecrd[1, 0] + sfac * d[3]])
    eyd = np.array([ecrd[0, 1] + sfac * d[1], ecrd[1, 1] + sfac * d[4]])

    return exd, eyd


def beam_defo_interp_2d(ecrd, u, sfac, nep=17):
    """
    Interpolate element displacements at nep points.

    Parametrs:
    ecrd : element x, y coordinates (2 x ndm),
    u : element nodal displacements
    sfac : scale factor for deformation plot
    nep : number of evaluation points (including end nodes)

    Returns:
    crd_xc, crd_yc : x, y coordinates of interpolated (at nep points)
    beam deformation required for plot_defo() function
    """


    G, L, cosa, cosb = opsvmodel.rot_transf_2d(ecrd)

    u_l = G @ u

    # longitudinal deformation (1)
    N_a = beam_axial_shape_functions(L, nep)
    u_ac = N_a @ np.array([u_l[0], u_l[3]])

    # transverse deformation (2)
    N_t = beam_transverse_shape_functions(L, nep)
    u_tc = N_t @ np.array([u_l[1], u_l[2], u_l[4], u_l[5]])
    # for u_tci in u_tc:

    # combined two row vectors
    # 1-st vector longitudinal deformation (1)
    # 2-nd vector transverse deformation (2)
    u_atc = np.vstack((u_ac, u_tc))

    # project longitudinal (u_ac) and transverse deformation
    # (local u and v) to (global u and v)
    G1 = np.array([[cosa, -cosb],
                   [cosb, cosa]])

    u_xyc = G1 @ u_atc

    # discretize element coordinates
    # first  row = X + [0 dx 2dx ... 4-dx 4]
    # second row = Y + [0 dy 2dy ... 3-dy 3]
    xy_c = np.vstack((np.linspace(ecrd[0, 0], ecrd[1, 0], num=nep),
                      np.linspace(ecrd[0, 1], ecrd[1, 1], num=nep)))

    # Continuous x, y displacement coordinates
    crd_xc = xy_c[0, :] + sfac * u_xyc[0, :]
    crd_yc = xy_c[1, :] + sfac * u_xyc[1, :]
    # latex_array(ecrd_xc)
    # latex_array(ecrd_yc)

    return crd_xc, crd_yc


def beam_disp_ends3d(ecrd, d, sfac):
    """
    Calculate the element deformation at element ends only.
    """

    #  indx: 0   1   2   3   4   5   6   7   8   9  10  11
    # Ed = ux1 uy1 uz1 rx1 ry1 rz1 ux2 uy2 uz2 rx2 ry2 rz2
    exd = np.array([ecrd[0, 0] + sfac * d[0], ecrd[1, 0] + sfac * d[6]])
    eyd = np.array([ecrd[0, 1] + sfac * d[1], ecrd[1, 1] + sfac * d[7]])
    ezd = np.array([ecrd[0, 2] + sfac * d[2], ecrd[1, 2] + sfac * d[8]])

    return exd, eyd, ezd


def beam_defo_interp_3d(ecrd, g, u, sfac, nep=17):
    """
    3d beam version of beam_defo_interp_2d.
    """
    G, L = opsvmodel.rot_transf_3d(ecrd, g)
    ul = G @ u

    _, crd_yc = beam_defo_interp_2d(np.array([[0., 0.], [L, 0.]]),
                                    np.array([ul[0], ul[1], ul[5], ul[6],
                                              ul[7], ul[11]]), sfac, nep)
    crd_xc, crd_zc = beam_defo_interp_2d(np.array([[0., 0.], [L, 0.]]),
                                         np.array([ul[0], ul[2], -ul[4], ul[6],
                                                   ul[8], -ul[10]]), sfac, nep)

    xl = np.linspace(0., L, num=nep)
    crd_xc = crd_xc - xl

    crd_xyzc = np.vstack([crd_xc, crd_yc, crd_zc])

    u_xyzc = np.transpose(g) @ crd_xyzc

    xyz_c = np.vstack((np.linspace(ecrd[0, 0], ecrd[1, 0], num=nep),
                       np.linspace(ecrd[0, 1], ecrd[1, 1], num=nep),
                       np.linspace(ecrd[0, 2], ecrd[1, 2], num=nep)))

    crd_xc = xyz_c[0, :] + u_xyzc[0, :]
    crd_yc = xyz_c[1, :] + u_xyzc[1, :]
    crd_zc = xyz_c[2, :] + u_xyzc[2, :]

    return crd_xc, crd_yc, crd_zc


def max_u_abs_from_beam_defo_interp_2d(max_u_abs, nep):
    """
    Find maximum displacement from the interpolated tranverse deformation.

    Useful if the translational nodal displacement are small and transverse
    deformation is due to nodal rotations.
    """

    ele_tags = ops.getEleTags()

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]

        if (ele_classtag == EleClassTag.ElasticBeam2d or
            ele_classtag == EleClassTag.ForceBeamColumn2d or
            ele_classtag == EleClassTag.DispBeamColumn2d):

            nen, ndf = 2, 3
            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            for i, ele_node_tag in enumerate(ele_node_tags):
                ed[i, :] = ops.nodeDisp(ele_node_tag)

            u = ed.flatten()

            G, L, *_ = opsvmodel.rot_transf_2d(ecrd)

            u_l = G @ u

            N_t = beam_transverse_shape_functions(L, nep)

            u_tc = N_t @ np.array([u_l[1], u_l[2], u_l[4], u_l[5]])

            max_u_tc = max(np.abs(u_tc))

            max_u_abs = max(max_u_abs, np.abs(max_u_tc))

    return max_u_abs


def max_u_abs_from_beam_defo_interp_3d(max_u_abs, nep):
    """
    Find maximum displacement from the interpolated tranverse deformation.

    Useful if the translational nodal displacement are small and transverse
    deformation is due to nodal rotations.
    """

    nen, ndf = 2, 6

    ele_tags = ops.getEleTags()

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]

        if (ele_classtag == EleClassTag.ElasticBeam3d or
            ele_classtag == EleClassTag.ForceBeamColumn3d or
            ele_classtag == EleClassTag.DispBeamColumn3d):

            ele_node_tags = ops.eleNodes(ele_tag)

            g = np.vstack((ops.eleResponse(ele_tag, 'xlocal'),
                           ops.eleResponse(ele_tag, 'ylocal'),
                           ops.eleResponse(ele_tag, 'zlocal')))

            ecrd = np.zeros((nen, 3))
            ed = np.zeros((nen, ndf))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            for i, ele_node_tag in enumerate(ele_node_tags):
                ed[i, :] = ops.nodeDisp(ele_node_tag)

            u = ed.flatten()

            G, L = opsvmodel.rot_transf_3d(ecrd, g)
            ul = G @ u

            N_t = beam_transverse_shape_functions(L, nep)

            # account for two cross section axis
            # (1) first cross section axis
            u_tc = N_t @ np.array([ul[1], ul[5], ul[7], ul[11]])
            max_u_tc = max(np.abs(u_tc))
            max_u_abs = max(max_u_abs, np.abs(max_u_tc))


            # (2) second cross section axis
            u_tc = N_t @ np.array([ul[2], -ul[4], ul[8], -ul[10]])
            max_u_tc = max(np.abs(u_tc))
            max_u_abs = max(max_u_abs, np.abs(max_u_tc))


    return max_u_abs


def beam_transverse_shape_functions(L, nep):

    xl = np.linspace(0., L, num=nep)
    one = np.ones(xl.shape)
    N_t = np.column_stack((one - 3 * xl**2 / L**2 + 2 * xl**3 / L**3,
                           xl - 2 * xl**2 / L + xl**3 / L**2,
                           3 * xl**2 / L**2 - 2 * xl**3 / L**3,
                           -xl**2 / L + xl**3 / L**2))

    return N_t


def beam_axial_shape_functions(L, nep):

    xl = np.linspace(0., L, num=nep)
    one = np.ones(xl.shape)
    N_a = np.column_stack((one - xl / L, xl / L))

    return N_a
