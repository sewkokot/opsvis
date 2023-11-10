import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle, Polygon, Wedge
from matplotlib.animation import FuncAnimation
from matplotlib.path import Path
import matplotlib.tri as tri

from .settings import *
from .defo import *


def _plot_model_2d(node_labels, element_labels, offset_nd_label, axis_off,
                   fig_wi_he, fig_lbrt, nodes_only, fmt_model,
                   fmt_model_nodes_only, node_supports, gauss_points, fmt_gauss_points,
                   fmt_model_truss, truss_node_offset, ax):

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he
            fig, ax = plt.subplots(figsize=(fig_wi / 2.54, fig_he / 2.54))
        else:
            fig, ax = plt.subplots()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

    max_x_crd, max_y_crd, max_crd = -np.inf, -np.inf, -np.inf

    node_tags = ops.getNodeTags()
    ele_tags = ops.getEleTags()

    # nodes only in the OpenSees domain
    if len(ele_tags) == 0 or nodes_only:
        for node_tag in node_tags:
            ax.plot(ops.nodeCoord(node_tag)[0],
                    ops.nodeCoord(node_tag)[1],
                    **fmt_model_nodes_only)
            if node_labels:
                ax.text(ops.nodeCoord(node_tag)[0],
                        ops.nodeCoord(node_tag)[1],
                        f'{node_tag}', va='bottom', ha='left', color='blue')

        ax.axis('equal')
        ax.grid(False)

        return ax

    for node_tag in node_tags:
        x_crd = ops.nodeCoord(node_tag)[0]
        y_crd = ops.nodeCoord(node_tag)[1]
        if x_crd > max_x_crd:
            max_x_crd = x_crd
        if y_crd > max_y_crd:
            max_y_crd = y_crd

    max_crd = np.amax([max_x_crd, max_y_crd])
    _offset = 0.005 * max_crd
    _offnl = 0.003 * max_crd

    # 0. first print node labels
    if node_labels:
        for j, node_tag in enumerate(node_tags):
            # xycrd = ops.nodeCoord(node_tag)
            if not offset_nd_label == 'above':
                offset_nd_label_x, offset_nd_label_y = _offnl, _offnl
                va = 'bottom'
                # va = 'center'
                ha = 'left'
            else:
                offset_nd_label_x, offset_nd_label_y = 0.0, _offnl
                va = 'bottom'
                ha = 'center'

            ax.text(ops.nodeCoord(node_tag)[0]+offset_nd_label_x,
                    ops.nodeCoord(node_tag)[1]+offset_nd_label_y,
                    f'{node_tag}', va=va, ha=ha, color='blue')

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]

        if (ele_classtag == EleClassTag.ElasticBeam2d
                or ele_classtag == EleClassTag.ForceBeamColumn2d
                or ele_classtag == EleClassTag.DispBeamColumn2d
                or ele_classtag == EleClassTag.TimoshenkoBeamColumn2d
                or ele_classtag == EleClassTag.ElasticTimoshenkoBeam2d
                or ele_classtag == EleClassTag.MVLEM
                or ele_classtag == EleClassTag.SFI_MVLEM):

            nen = 2
            ele_node_tags = ops.eleNodes(ele_tag)
            ecrd_nodes = np.zeros((nen, 2))
            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd_nodes[i, :] = ops.nodeCoord(ele_node_tag)

            ecrd_eles = np.copy(ecrd_nodes)

            ele_offsets = np.array(ops.eleResponse(ele_tag, 'offsets'))
            nz_offsets = np.nonzero(ele_offsets)[0]  # tuple of arrays
            if np.any(nz_offsets):
                # modify ecrd_eles
                ecrd_eles[:, 0] += ele_offsets[[0, 3]]
                ecrd_eles[:, 1] += ele_offsets[[1, 4]]
                ax.plot([ecrd_nodes[0, 0], ecrd_eles[0, 0]],
                        [ecrd_nodes[0, 1], ecrd_eles[0, 1]],
                        **fmt_model_rigid_offset)
                ax.plot([ecrd_nodes[1, 0], ecrd_eles[1, 0]],
                        [ecrd_nodes[1, 1], ecrd_eles[1, 1]],
                        **fmt_model_rigid_offset)

            else:
                pass

            ax.plot(ecrd_eles[:, 0], ecrd_eles[:, 1], **fmt_model)

            if gauss_points:
                Lgps = ops.eleResponse(ele_tag, 'integrationPoints')
                for Lgpi in Lgps:
                    dxi = (ecrd_eles[1, 0] - ecrd_eles[0, 0]) * Lgpi / bar_length(ecrd_eles)
                    dyi = (ecrd_eles[1, 1] - ecrd_eles[0, 1]) * Lgpi / bar_length(ecrd_eles)
                    ax.plot(ecrd_eles[0, 0] + dxi, ecrd_eles[0, 1] + dyi, **fmt_gauss_points)

            # location of label
            xt = sum(ecrd_eles[:, 0]) / nen
            yt = sum(ecrd_eles[:, 1]) / nen

            if element_labels:
                if ecrd_eles[1, 0] - ecrd_eles[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y = _offset, 0.0
                elif ecrd_eles[1, 1] - ecrd_eles[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y = 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y = 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, f'{ele_tag}', va=va, ha=ha,
                        color='red')

        if (ele_classtag == EleClassTag.truss
                or ele_classtag == EleClassTag.trussSection):

            nen = 2
            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)


            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            if truss_node_offset:
                Lx = ecrd[1, 0] - ecrd[0, 0]
                Ly = ecrd[1, 1] - ecrd[0, 1]
                dLx = truss_node_offset * Lx
                dLy = truss_node_offset * Ly

                ecrd[0, 0] += dLx
                ecrd[1, 0] -= dLx
                ecrd[0, 1] += dLy
                ecrd[1, 1] -= dLy

            ax.plot(ecrd[:, 0], ecrd[:, 1], **fmt_model_truss)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y = _offset, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y = 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y = 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, f'{ele_tag}', va=va, ha=ha,
                        color='red')

        # 2d triangular (tri31) elements plot_model
        elif (ele_classtag == EleClassTag.tri3n):

            nen = 3
            nodes_geo_order = [0, 1, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1], **fmt_model)

            if element_labels:
                va = 'center'
                ha = 'center'
                ax.text(xt, yt, f'{ele_tag}', va=va, ha=ha, color='red')

        # 2d quadrilateral (quad) elements plot_model
        elif (ele_classtag == EleClassTag.quad4n):

            nen = 4
            nodes_geo_order = [0, 1, 2, 3, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1], **fmt_model)

            if element_labels:
                va = 'center'
                ha = 'center'
                ax.text(xt, yt, f'{ele_tag}', va=va, ha=ha, color='red')

        # 2d quadrilateral (quad8n) elements plot_model
        elif (ele_classtag == EleClassTag.quad8n):

            nen = 8
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 3, 7, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1], **fmt_model)

            if element_labels:
                va = 'center'
                ha = 'center'
                ax.text(xt, yt, f'{ele_tag}', va=va, ha=ha, color='red')

        # 2d quadrilateral (quad9n) elements plot_model
        elif (ele_classtag == EleClassTag.quad9n):

            nen = 9
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 3, 7, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1], **fmt_model)

            if element_labels:
                va = 'center'
                ha = 'center'
                ax.text(xt, yt, f'{ele_tag}', va=va, ha=ha, color='red')

            ax.scatter([ecrd[8, 0]], [ecrd[8, 1]], s=2, color='g')

        # 2d triangle (tri6n) elements plot_model
        elif (ele_classtag == EleClassTag.tri6n):

            nen = 6
            nodes_geo_order = [0, 3, 1, 4, 2, 5, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1], **fmt_model)

            if element_labels:
                va = 'center'
                ha = 'center'
                ax.text(xt, yt, f'{ele_tag}', va=va, ha=ha, color='red')

        # 2d joint2d element plot_model
        elif (ele_classtag == EleClassTag.Joint2D):

            nen = 5

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label - central node
            xt = ecrd[4, 0]
            yt = ecrd[4, 1]

            ax.plot([ecrd[0, 0], ecrd[2, 0], ecrd[2, 0], ecrd[0, 0],
                     ecrd[0, 0]],
                    [ecrd[1, 1], ecrd[1, 1], ecrd[3, 1], ecrd[3, 1],
                     ecrd[1, 1]], **fmt_model_joint2d)

            if element_labels:
                ax.text(xt, yt, f'{ele_tag}', va='center', ha='center', color='red')

        # 2d zeroLength, two node link element plot_model
        elif (ele_classtag == EleClassTag.ZeroLength or
              ele_classtag == EleClassTag.TwoNodeLink):

            nen = 2

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 2))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = ecrd[0, 0]
            yt = ecrd[0, 1]

            # ax.plot([ecrd[0, 0], ecrd[2, 0], ecrd[2, 0], ecrd[0, 0],
            #          ecrd[0, 0]],
            #         [ecrd[1, 1], ecrd[1, 1], ecrd[3, 1], ecrd[3, 1],
            #          ecrd[1, 1]], **fmt_model_joint2d)

            if element_labels:
                ax.text(xt, yt, f'{ele_tag}', va='top', ha='left',
                        color='red')

    if node_supports:
        _plot_supports(node_tags, ax)

    ax.axis('equal')
    ax.grid(False)

    return ax


def _plot_supports(node_tags, ax):

    ndim = ops.getNDM()[0]

    # fix 2d support: square
    verts = [(-1., 0.),
             (1., 0.),
             (1., -2.),
             (-1., -2.),
             (-1., 0.)]

    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY]

    path_fix = Path(verts, codes)

    # hinge/pin 2d support: ^
    verts = [(0., 0.),
             (-0.2, -0.25),
             (0.2, -0.25),
             (0., 0.)]

    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY]

    path_pin = Path(verts, codes)

    # roller horiz 2d support: o_
    verts = [(0., 0.),
             (1., 0.),
             (1., -1.),
             (1., -2.),
             (0., -2.),
             (-1., -2.),
             (-1., -1.),
             (-1., 0.),
             (0., 0.),
             (-1.5, -2.),
             (1.5, -2.)]

    codes = [Path.MOVETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.MOVETO,
             Path.LINETO]

    path_roller_horiz = Path(verts, codes)

    # roller vert 2d support: o|
    verts = [(0., 0.),
             (0., 1.),
             (1., 1.),
             (2., 1.),
             (2., 0.),
             (2., -1.),
             (1., -1.),
             (0., -1.),
             (0., 0.),
             (2, -1.5),
             (2, 1.5)]

    codes = [Path.MOVETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.CURVE3,
             Path.LINETO,
             Path.MOVETO,
             Path.LINETO]

    path_roller_vert = Path(verts, codes)

    m_color = 'm'
    m_fstyle = 'full'
    m_size = 16

    fixed_nodes = ops.getFixedNodes()

    for node_tag in fixed_nodes:
        nd_crd = ops.nodeCoord(node_tag)

        fixed_dofs = ops.getFixedDOFs(node_tag)

        if ndim == 2:
            if (fixed_dofs == [1, 2, 3]):
                m_type = path_fix
            elif (fixed_dofs == [1, 2]):
                m_type = path_pin
            elif (fixed_dofs == [1]):
                m_type = path_roller_vert
            elif (fixed_dofs == [2]):
                m_type = path_roller_horiz

            ax.plot(nd_crd[0], nd_crd[1], marker=m_type, markersize=m_size,
                    color=m_color, fillstyle=m_fstyle)

        elif ndim == 3:
            if (fixed_dofs == [1, 2, 3, 4, 5, 6]):
                m_type = path_fix
            elif (fixed_dofs == [1, 2, 3]):
                m_type = path_pin
            else:
                m_type = None

            ax.plot(nd_crd[0], nd_crd[1], nd_crd[2], marker=m_type, markersize=m_size,
                    color=m_color, fillstyle=m_fstyle)

    return ax


def _plot_model_3d(node_labels, element_labels, offset_nd_label,
                   axis_off, az_el, fig_wi_he, fig_lbrt, local_axes,
                   nodes_only, fmt_model, node_supports, gauss_points,
                   fmt_gauss_points, fmt_model_truss,
                   truss_node_offset, ax):

    node_tags = ops.getNodeTags()
    ele_tags = ops.getEleTags()

    azim, elev = az_el

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he

            fig = plt.figure(figsize=(fig_wi / 2.54, fig_he / 2.54))
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

    max_x_crd, max_y_crd, max_z_crd, max_crd = -np.inf, -np.inf, \
        -np.inf, -np.inf

    # nodes only in the OpenSees domain
    if len(ele_tags) == 0 or nodes_only:
        for node_tag in node_tags:
            ax.plot(ops.nodeCoord(node_tag)[0],
                    ops.nodeCoord(node_tag)[1],
                    ops.nodeCoord(node_tag)[2],
                    color='green', marker='o')
            if node_labels:
                ax.text(ops.nodeCoord(node_tag)[0],
                        ops.nodeCoord(node_tag)[1],
                        ops.nodeCoord(node_tag)[2],
                        f'{node_tag}', va='bottom', ha='left', color='blue')

        return ax

    for node_tag in node_tags:
        x_crd = ops.nodeCoord(node_tag)[0]
        y_crd = ops.nodeCoord(node_tag)[1]
        z_crd = ops.nodeCoord(node_tag)[2]
        if x_crd > max_x_crd:
            max_x_crd = x_crd
        if y_crd > max_y_crd:
            max_y_crd = y_crd
        if z_crd > max_z_crd:
            max_z_crd = z_crd

    if offset_nd_label == 0 or offset_nd_label == 0.:
        _offset = 0.
    else:
        max_crd = np.amax([max_x_crd, max_y_crd, max_z_crd])
        _offset = 0.005 * max_crd

    if node_labels:
        for node_tag in node_tags:
            ax.text(ops.nodeCoord(node_tag)[0] + _offset,
                    ops.nodeCoord(node_tag)[1] + _offset,
                    ops.nodeCoord(node_tag)[2] + _offset,
                    f'{node_tag}', va='bottom', ha='left', color='blue')

    for i, ele_tag in enumerate(ele_tags):
        ele_classtag = ops.getEleClassTags(ele_tag)[0]
        nen = np.shape(ops.eleNodes(ele_tag))[0]

        if (ele_classtag == EleClassTag.ElasticBeam3d or
            ele_classtag == EleClassTag.ForceBeamColumn3d or
            ele_classtag == EleClassTag.DispBeamColumn3d or
            ele_classtag == EleClassTag.ElasticTimoshenkoBeam3d):

            nen = 2
            ele_node_tags = ops.eleNodes(ele_tag)
            ecrd_nodes = np.zeros((nen, 3))
            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd_nodes[i, :] = ops.nodeCoord(ele_node_tag)

            ecrd_eles = np.copy(ecrd_nodes)

            ele_offsets = np.array(ops.eleResponse(ele_tag, 'offsets'))
            nz_offsets = np.nonzero(ele_offsets)[0]  # tuple of arrays
            if np.any(nz_offsets):
                # modify ecrd_eles
                ecrd_eles[:, 0] += ele_offsets[[0, 3]]
                ecrd_eles[:, 1] += ele_offsets[[1, 4]]
                ecrd_eles[:, 2] += ele_offsets[[2, 5]]
                ax.plot([ecrd_nodes[0, 0], ecrd_eles[0, 0]],
                        [ecrd_nodes[0, 1], ecrd_eles[0, 1]],
                        [ecrd_nodes[0, 2], ecrd_eles[0, 2]],
                        **fmt_model_rigid_offset)
                ax.plot([ecrd_nodes[1, 0], ecrd_eles[1, 0]],
                        [ecrd_nodes[1, 1], ecrd_eles[1, 1]],
                        [ecrd_nodes[1, 2], ecrd_eles[1, 2]],
                        **fmt_model_rigid_offset)

            else:
                pass

            ax.plot(ecrd_eles[:, 0],
                    ecrd_eles[:, 1],
                    ecrd_eles[:, 2], **fmt_model)
            # ax.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], **fmt_model)

            # location of label
            xt = sum(ecrd_eles[:, 0]) / nen
            yt = sum(ecrd_eles[:, 1]) / nen
            zt = sum(ecrd_eles[:, 2]) / nen

            # fixme: placement of node_tag labels
            if element_labels:
                if ecrd_eles[1, 0] - ecrd_eles[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd_eles[1, 1] - ecrd_eles[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd_eles[1, 2] - ecrd_eles[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt + offset_x,
                        yt + offset_y,
                        zt + offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            if local_axes:
                # eo = Eo[i, :]
                xloc = ops.eleResponse(ele_tag, 'xlocal')
                yloc = ops.eleResponse(ele_tag, 'ylocal')
                zloc = ops.eleResponse(ele_tag, 'zlocal')
                g = np.vstack((xloc, yloc, zloc))
                L = bar_length(ecrd_eles)
                alen = 0.1 * L

                # plot local axis directional vectors: workaround quiver = arrow
                ax.quiver(xt, yt, zt, g[0, 0], g[0, 1], g[0, 2], color='g',
                          lw=2, length=alen, alpha=.8, normalize=True)
                ax.quiver(xt, yt, zt, g[1, 0], g[1, 1], g[1, 2], color='b',
                          lw=2, length=alen, alpha=.8, normalize=True)
                ax.quiver(xt, yt, zt, g[2, 0], g[2, 1], g[2, 2], color='r',
                          lw=2, length=alen, alpha=.8, normalize=True)

            if gauss_points:
                Lgps = ops.eleResponse(ele_tag, 'integrationPoints')
                for Lgpi in Lgps:
                    dxi = (ecrd_eles[1, 0] - ecrd_eles[0, 0]) * Lgpi / bar_length(ecrd_eles)
                    dyi = (ecrd_eles[1, 1] - ecrd_eles[0, 1]) * Lgpi / bar_length(ecrd_eles)
                    dzi = (ecrd_eles[1, 2] - ecrd_eles[0, 2]) * Lgpi / bar_length(ecrd_eles)
                    ax.plot(ecrd_eles[0, 0] + dxi, ecrd_eles[0, 1] + dyi, ecrd_eles[0, 2] + dzi, **fmt_gauss_points)

        elif (ele_classtag == EleClassTag.truss
              or ele_classtag == EleClassTag.trussSection):

            nen = 2
            ele_node_tags = ops.eleNodes(ele_tag)
            ecrd = np.zeros((nen, 3))
            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            if truss_node_offset:
                Lx = ecrd[1, 0] - ecrd[0, 0]
                Ly = ecrd[1, 1] - ecrd[0, 1]
                Lz = ecrd[1, 2] - ecrd[0, 2]
                dLx = truss_node_offset * Lx
                dLy = truss_node_offset * Ly
                dLz = truss_node_offset * Lz

                ecrd[0, 0] += dLx
                ecrd[1, 0] -= dLx
                ecrd[0, 1] += dLy
                ecrd[1, 1] -= dLy
                ecrd[0, 2] += dLz
                ecrd[1, 2] -= dLz

            ax.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], **fmt_model_truss)

            # fixme: placement of node_tag labels
            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

        # quad in 3d
        # elif nen == 4:
        elif ele_classtag == EleClassTag.FourNodeTetrahedron:

            nen = 4
            nodes_geo_order = [0, 1, 2, 0]
            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model)

            for j in range(3):
                ax.plot(ecrd[[j, 3], 0],
                        ecrd[[j, 3], 1],
                        ecrd[[j, 3], 2], **fmt_model)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        elif (ele_classtag in {EleClassTag.TenNodeTetrahedron,
                               EleClassTag.TenNodeTetrahedronSK}):

            nen = 10
            nodes_geo_order = [0, 4, 1, 5, 2, 6, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    'b.-')

            for j in range(3):
                ax.plot(ecrd[[j, 7+j, 3], 0],
                        ecrd[[j, 7+j, 3], 1],
                        ecrd[[j, 7+j, 3], 2], 'b.-')

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        elif (ele_classtag == EleClassTag.ShellMITC4 or
              ele_classtag == EleClassTag.ASDShellQ4 or
              ele_classtag == EleClassTag.ShellNLDKGQ or
              ele_classtag == EleClassTag.ShellDKGQ or                            
              ele_classtag == EleClassTag.quad4n):

            nen = 4
            nodes_geo_order = [0, 1, 2, 3, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        elif (ele_classtag == EleClassTag.ShellNLDKGT or
              ele_classtag == EleClassTag.ShellDKGT):
            
            nen = 3
            nodes_geo_order = [0, 1, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')
            
        # mvlem 3d plot model
        elif (ele_classtag == EleClassTag.MVLEM_3D
              or ele_classtag == EleClassTag.SFI_MVLEM_3D):

            nen = 4
            nodes_geo_order = [0, 1, 3, 2, 0]

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        # 8-node brick, 3d model
        # elif nen == 8:
        elif (ele_classtag == EleClassTag.brick8n or
              ele_classtag == EleClassTag.SSPbrick):

            nen = 8
            nodes_geo_order = np.array([0, 1, 2, 3, 0], dtype=int)  # bottom face nodes loop

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model)
            ax.plot(ecrd[nodes_geo_order + 4, 0],
                    ecrd[nodes_geo_order + 4, 1],
                    ecrd[nodes_geo_order + 4, 2],
                    **fmt_model)

            for j in range(4):
                ax.plot(ecrd[[j, j+4], 0],
                        ecrd[[j, j+4], 1],
                        ecrd[[j, j+4], 2], **fmt_model)

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        # 20-node brick, 3d model
        # elif nen == 20:
        elif (ele_classtag == EleClassTag.brick20n):

            nen = 20
            nodes_geo_order = np.array([0, 8, 1, 9, 2, 10, 3, 11, 0], dtype=int)  # bottom face nodes loop

            ele_node_tags = ops.eleNodes(ele_tag)
            ecrd = np.zeros((nen, 3))
            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            ax.plot(ecrd[nodes_geo_order, 0],
                    ecrd[nodes_geo_order, 1],
                    ecrd[nodes_geo_order, 2],
                    **fmt_model, mfc='g', mec='g')
            ax.plot(ecrd[nodes_geo_order + 4, 0],
                    ecrd[nodes_geo_order + 4, 1],
                    ecrd[nodes_geo_order + 4, 2],
                    **fmt_model, mfc='g', mec='g')

            # for j in range(2):
            #     ax.plot(ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 0],
            #             ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 1],
            #             ecrd[[0+j*4, 8+j*4, 1+j*4, 9+j*4, 2+j*4, 10+j*4, 3+j*4, 11+j*4, 0+j*4], 2],
            #             **fmt_model, mfc='g', mec='g')

            for j in range(4):
                ax.plot(ecrd[[j, 16+j, 4+j], 0],
                        ecrd[[j, 16+j, 4+j], 1],
                        ecrd[[j, 16+j, 4+j], 2], **fmt_model,
                        mfc='g', mec='g')

            # location of label
            xt = sum(ecrd[:, 0]) / nen
            yt = sum(ecrd[:, 1]) / nen
            zt = sum(ecrd[:, 2]) / nen

            if element_labels:
                if ecrd[1, 0] - ecrd[0, 0] == 0:
                    va = 'center'
                    ha = 'left'
                    offset_x, offset_y, offset_z = _offset, 0.0, 0.0
                elif ecrd[1, 1] - ecrd[0, 1] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, _offset, 0.0
                elif ecrd[1, 2] - ecrd[0, 2] == 0:
                    va = 'bottom'
                    ha = 'center'
                    offset_x, offset_y, offset_z = 0.0, 0.0, _offset
                else:
                    va = 'bottom'
                    ha = 'left'
                    offset_x, offset_y, offset_z = 0.03, 0.03, 0.03

                ax.text(xt+offset_x, yt+offset_y, zt+offset_z, f'{ele_tag}',
                        va=va, ha=ha, color='red')

            # if node_labels:
            #     for node_tag in node_tags:
            #         ax.text(ops.nodeCoord(node_tag)[0]+_offset,
            #                 ops.nodeCoord(node_tag)[1]+_offset,
            #                 ops.nodeCoord(node_tag)[2]+_offset,
            #                 f'{node_tag}', va='bottom', ha='left', color='blue')

        # 3d zeroLength, two node link element plot_model
        elif (ele_classtag == EleClassTag.ZeroLength or
              ele_classtag == EleClassTag.TwoNodeLink):

            nen = 2

            ele_node_tags = ops.eleNodes(ele_tag)

            ecrd = np.zeros((nen, 3))

            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd[i, :] = ops.nodeCoord(ele_node_tag)

            # location of label
            xt = ecrd[0, 0]
            yt = ecrd[0, 1]
            zt = ecrd[0, 2]

            # ax.plot([ecrd[0, 0], ecrd[2, 0], ecrd[2, 0], ecrd[0, 0],
            #          ecrd[0, 0]],
            #         [ecrd[1, 1], ecrd[1, 1], ecrd[3, 1], ecrd[3, 1],
            #          ecrd[1, 1]], **fmt_model_joint2d)

            if element_labels:
                ax.text(xt, yt, zt, f'{ele_tag}', va='top', ha='left',
                        color='red')

    if node_supports:
        _plot_supports(node_tags, ax)

    ax.set_box_aspect((np.ptp(ax.get_xlim3d()),
                       np.ptp(ax.get_ylim3d()),
                       np.ptp(ax.get_zlim3d())))

    return ax


def plot_model(node_labels=1, element_labels=1, offset_nd_label=False,
               axis_off=0, az_el=az_el, fig_wi_he=False,
               fig_lbrt=False, local_axes=True, nodes_only=False,
               fmt_model=fmt_model,
               fmt_model_nodes_only=fmt_model_nodes_only,
               node_supports=True, gauss_points=True,
               fmt_gauss_points=fmt_gauss_points,
               fmt_model_truss=fmt_model_truss,
               truss_node_offset=0.96, ax=False):
    """Plot defined model of the structure.

    Args:
        node_labels (int): 1 - plot node labels, 0 - do not plot them;
            (default: 1)

        element_labels (int): 1 - plot element labels, 0 - do not plot
            them; (default: 1)

        offset_nd_label (bool): False - do not offset node labels from the
            actual node location. This option can enhance visibility.

        axis_off (int): 0 - turn off axes, 1 - display axes; (default: 0)

        az_el (tuple): contains azimuth and elevation for 3d plots. For 2d
            plots this parameter is neglected.

        fig_wi_he (tuple): contains width and height of the figure

        fig_lbrt (tuple): a tuple contating left, bottom, right and top offsets

        local_axes (bool): True - show cross section local axes or False.
            The green, red and blue arrows denote the element axis direction,
            the z-local axis and the y-local axis.

        nodes_only (bool): True - show the nodes only, although the elements
            are defined. Default: False.

        fmt_model (dict): A dictionary containing formatting the line and markers
            of the model elements. The formatting options can be: linewidth, color,
            marker, markersize. See matplotlib.plot documentation for more details,
            if necessary.

        node_supports (bool): True - show the supports. Default: True.

        gauss_points (bool): True - show the integration (Gauss) points. Default: True.

        fmt_gauss_points (dict): A dictionary containing formatting the marker
            of the gauss point.

        truss_node_offset (float): If non-zero, the nodes are offset to show pin markers.
            The number should preferably be between 0.90-0.97.
            Zero (0) or False means no offset. Default: 0.96.

        ax: axis object.

    Usage:

    ``plot_model()`` - plot model with node and element labels.

    ``plot_model(node_labels=0, element_labels=0)`` - plot model without node
    element labels

    ``plot_model(fig_wi_he=(20., 14.))`` - plot model in a window 20 cm long,
    and 14 cm high.

    ``plot_model(nodes_only=True)`` - plot only the nodes even thought the elements
    are defined.

    ``plot_model(node_supports=False)`` - plot the model without the supports.

    """

    # az_el - azimut, elevation used for 3d plots only
    node_tags = ops.getNodeTags()

    ndim = ops.getNDM()[0]

    if ndim == 2:
        ax = _plot_model_2d(node_labels, element_labels,
                            offset_nd_label, axis_off, fig_wi_he,
                            fig_lbrt, nodes_only, fmt_model,
                            fmt_model_nodes_only, node_supports,
                            gauss_points, fmt_gauss_points,
                            fmt_model_truss, truss_node_offset, ax)
        if axis_off:
            ax.axis('off')

    elif ndim == 3:
        ax = _plot_model_3d(node_labels, element_labels,
                            offset_nd_label, axis_off, az_el,
                            fig_wi_he, fig_lbrt, local_axes,
                            nodes_only, fmt_model, node_supports,
                            gauss_points, fmt_gauss_points,
                            fmt_model_truss, truss_node_offset, ax)
        if axis_off:
            ax.axis('off')

    else:
        print(f'\nWarning! ndim: {ndim} not supported yet.')

    # plt.show()  # call this from main py file for more control
    return ax


def plot_supports_and_loads_2d(nep=17):
    """This function is removed. Use plot_model() and
    plot_loads_2d() for showing support conditions and loads,
    respectively.
    """

    print('\nWarning! plot_supports_and_loads_2d has been removed.')  # noqa: E501
    print('\nWarning! Use plot_model(node_supports=True) for showing supports and')  # noqa: E501
    print('\nWarning! plot_loads_2d() for showing loads')  # noqa: E501


def plot_loads_2d(nep=17, sfac=False, fig_wi_he=False,
                  fig_lbrt=False, fmt_model_loads=fmt_model_loads,
                  node_supports=True, truss_node_offset=0, ax=False):
    """Display the nodal and element loads applied to the 2d models.

    Args:
        nep (int): number of arrows when displacing element distributed loads
            (default: 17)

        node_supports (bool): True - show the node support conditions.
            Default: False.

    """

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he
            fig, ax = plt.subplots(figsize=(fig_wi / 2.54, fig_he / 2.54))
        else:
            fig, ax = plt.subplots()

        if fig_lbrt:
            fleft, fbottom, fright, ftop = fig_lbrt
            fig.subplots_adjust(left=fleft, bottom=fbottom, right=fright, top=ftop)

    ax = plot_model(node_labels=0, element_labels=0, fmt_model=fmt_model_loads,
                    node_supports=node_supports,
                    truss_node_offset=truss_node_offset, ax=ax)
    # ax.axis('equal')

    if not sfac:
        ratio = 0.1
        min_x, max_x = ax.get_xlim()
        min_y, max_y = ax.get_ylim()
        xsfac = ratio * abs(max_x - min_x)
        ysfac = ratio * abs(max_y - min_y)
        sfac = max(xsfac, ysfac)

    node_tags = ops.getNodeTags()
    ele_tags = ops.getEleTags()

    Wx = False
    waa = False
    wab = False

    ### get Ew data
    Ew = get_Ew_data_from_ops_domain()

    for ele_tag in ele_tags:

        ele_class_tag = ops.getEleClassTags(ele_tag)[0]

        if (ele_class_tag == EleClassTag.ElasticBeam2d
            or ele_class_tag == EleClassTag.ForceBeamColumn2d
            or ele_class_tag == EleClassTag.DispBeamColumn2d
            or ele_class_tag == EleClassTag.TimoshenkoBeamColumn2d
            or ele_class_tag == EleClassTag.ElasticTimoshenkoBeam2d):

            # by default no element load
            ele_load_data = []
            if ele_tag in Ew:
                ele_load_data = Ew[ele_tag]

            nen = 2
            ele_node_tags = ops.eleNodes(ele_tag)
            ecrd_nodes = np.zeros((nen, 2))
            for i, ele_node_tag in enumerate(ele_node_tags):
                ecrd_nodes[i, :] = ops.nodeCoord(ele_node_tag)

            ecrd_eles = np.copy(ecrd_nodes)

            ele_offsets = np.array(ops.eleResponse(ele_tag, 'offsets'))
            nz_offsets = np.nonzero(ele_offsets)[0]  # tuple of arrays
            if np.any(nz_offsets):
                # modify ecrd_eles
                ecrd_eles[:, 0] += ele_offsets[[0, 3]]
                ecrd_eles[:, 1] += ele_offsets[[1, 4]]
            else:
                pass

            # step 1: first plot model itself
            # ax.plot(ecrd[:, 0], ecrd[:, 1], 'k-', solid_capstyle='round', solid_joinstyle='round',
            #         dash_capstyle='butt', dash_joinstyle='round')

            # step 2
            Lxy = ecrd_eles[1, :] - ecrd_eles[0, :]
            L = np.sqrt(Lxy @ Lxy)
            cosa, cosb = Lxy / L

            # if not sfac:
            #     sfac = ratio * L

            xl = np.linspace(0., L, nep)
            xl2 = np.linspace(0., L, 5)
            # xl2 = np.logspace(0., 1., nep) * L / nep
            dxs = np.array([0.1, 0.2, 0.3, 0.4])
            # dxs = np.array([0.15, 0.2, 0.275, 0.375])
            xl3 = np.append(0., np.cumsum(dxs)) * L

            one = np.ones(nep)
            # pl = ops.eleResponse(ele_tag, 'localForces')

            # s_all, xl, nep = section_force_distribution_2d(ex, ey, pl, nep, ele_load_data)

            for ele_load_data_i in ele_load_data:
                ele_load_type = ele_load_data_i[0]

                if ele_load_type == '-beamPoint':
                    Pt, aL, Pa = ele_load_data_i[1:4]
                    a = aL * L
                    s = sfac * np.sign(Pt)

                    s_0 = np.zeros(2)
                    # s_0 = [ecrd[0, 0], ecrd[0, 1]]

                    s_0[0] = ecrd_eles[0, 0] + a * cosa
                    s_0[1] = ecrd_eles[0, 1] + a * cosb

                    s_p = np.copy(s_0)

                    s_p[0] -= s * cosb
                    s_p[1] += s * cosa


                    # plot arrows
                    ax.arrow(s_0[0], s_0[1],
                             s_p[0] - s_0[0], s_p[1] - s_0[1],
                             # width = 0.01,
                             head_starts_at_zero=True,  # default False
                             # overhang=0.2,
                             # lw = 1,
                             head_width=0.1 * sfac, head_length=0.2 * sfac,
                             fc='b', ec='b',
                             length_includes_head=False, shape='full')
                    ax.text(sum(ecrd_eles[:, 0])/2, sum(ecrd_eles[:, 1])/2, f'{abs(Pt)}', color='b')
                    # ax.annotate("", xy=(s_p[0], s_p[1]),
                    #             xytext=(s_0[0], s_0[1]),
                    #             arrowprops=dict(arrowstyle="->", color='r', lw=3,
                    #                             connectionstyle="arc3"))

                elif ele_load_type == '-beamUniform':

                    n_ele_load_data = len(ele_load_data_i)

                    # constant uniform element load
                    if n_ele_load_data == 3:
                        # eload_type, Wy, Wx = ele_load_data[0], ele_load_data[1], ele_load_data[2]
                        Wy, Wx = ele_load_data_i[1], ele_load_data_i[2]
                        text_string = f'q = {abs(Wy)}, {abs(Wx)}'

                        # s = sfac * Wy * one
                        s = sfac * one * np.sign(Wy)
                        # plt.text(sum(ecrd[:, 0])/2, sum(ecrd[:, 1])/2, f'q = {Wy}, {Wx}', va='bottom', ha='center', color='r')

                    # triangular or trapezoidal element load
                    elif n_ele_load_data == 7:
                        wta, waa, aL, bL, wtb, wab = ele_load_data_i[1:7]
                        text_string = f'q = {abs(wta)}, {abs(waa)}, {aL}, {bL}, {abs(wtb)}, {abs(wab)}'
                        a, b = aL * L, bL * L
                        bma = b - a
                        s = np.zeros(nep)

                        indx = 0
                        for x in np.nditer(xl):
                            xma = x - a
                            wtx = wta + (wtb - wta) * xma / bma

                            if x < a:
                                pass
                            elif x >= a and x <= b:
                                s[indx] = wtx
                            elif x > b:
                                pass

                            indx += 1

                    s = s * sfac

                    s_0 = np.zeros((nep, 2))
                    # s_0[0, :] = [ex[0], ey[0]]
                    s_0[0, :] = [ecrd_eles[0, 0], ecrd_eles[0, 1]]

                    s_0[1:, 0] = s_0[0, 0] + xl[1:] * cosa
                    s_0[1:, 1] = s_0[0, 1] + xl[1:] * cosb

                    s_p = np.copy(s_0)

                    s_p[:, 0] -= s * cosb
                    s_p[:, 1] += s * cosa
                    # plt.axis('equal')

                    # reference perpendicular lines
                    for i in np.arange(nep):
                        ax.arrow(s_0[i, 0], s_0[i, 1],
                                 s_p[i, 0] - s_0[i, 0], s_p[i, 1] - s_0[i, 1],
                                 # width = 0.005,
                                 # lw = 2,
                                 head_width=0.1*sfac, head_length=0.2*sfac,
                                 head_starts_at_zero=True,  # default False
                                 # overhang=0.5,
                                 fc='r', ec='r',
                                 length_includes_head=True, shape='full')
                        # ax.annotate("", xy=(s_p[i, 0], s_p[i, 1]),
                        #             xytext=(s_0[i, 0], s_0[i, 1]),
                        #             arrowprops=dict(arrowstyle="->", color='r', lw=1,
                        #                             connectionstyle="arc3"))

                    # connecting beg-end line - redundant ?
                    # plt.plot([s_p[0, 0], s_p[-1, 0]],[s_p[0, 1], s_p[-1, 1]], 'r')
                    ax.text(sum(ecrd_eles[:, 0])/2, sum(ecrd_eles[:, 1])/2, text_string, va='bottom', ha='center', color='r')

                    if Wx != 0:
                        # for i, xl in enumerate(xl2[:-1]):
                        #     plt.arrow(sa[i, 0], sa[i, 1],
                        #               sa[i+1, 0]-sa[i, 0], sa[i+1, 1]-sa[i, 1],
                        #               width = 0.01,
                        #               # lw = 1,
                        #               # head_width=0.02, head_length=0.05,
                        #               # head_starts_at_zero=True,  # default False
                        #               # overhang=0.5,
                        #               fc='m', ec='m',
                        #               length_includes_head=True, shape='full')

                        ax.quiver(s_0[:-1, 0], s_0[:-1, 1],
                                  s_0[1:, 0]-s_0[:-1, 0], s_0[1:, 1]-s_0[:-1, 1],
                                  scale_units='xy', angles='xy', scale=0.8, color='m')

                    if waa != 0 or wab != 0:
                        sa = np.zeros((5, 2))
                        # sa[0, :] = [ex[0], ey[0]]
                        sa[0, :] = [ecrd_eles[0, 0], ecrd_eles[0, 1]]
                        sa[1:, 0] = sa[0, 0] + xl3[1:] * cosa
                        sa[1:, 1] = sa[0, 1] + xl3[1:] * cosb

                        for i, xl in enumerate(xl3[:-1]):
                            ax.arrow(sa[i, 0], sa[i, 1],
                                     sa[i+1, 0]-sa[i, 0], sa[i+1, 1]-sa[i, 1],
                                     # width = 0.05,
                                     width = 0.01*L,
                                     # lw = 1,
                                     # head_width=0.02, head_length=0.05,
                                     # head_starts_at_zero=True,  # default False
                                     # overhang=0.5,
                                     fc='m', ec='m',
                                     length_includes_head=True, shape='full')

                        # ax.quiver(sa[:-1, 0], sa[:-1, 1],
                        #            sa[1:, 0]-sa[:-1, 0], sa[1:, 1]-sa[:-1, 1],
                        #            scale_units='xy', angles='xy', scale=0.8, color='g')

                        # ax.plot([s_0[i, 0], s_p[i, 0]], [s_0[i, 1], s_p[i, 1]],
                        #          fmt_secforce, solid_capstyle='round',
                        #          solid_joinstyle='round', dash_capstyle='butt',
                        #          dash_joinstyle='round')
                        # plot arrows
                        # ax.annotate("",
                        #              xy=(s_p[i, 0], s_p[i, 1]), xycoords='data',
                        #              xytext=(s_0[i, 0], s_0[i, 1]), textcoords='data',
                        #              arrowprops=dict(arrowstyle="->", color='r', lw=2,
                        #                              connectionstyle="arc3"))

    for node_tag in node_tags:
        nd_crd = ops.nodeCoord(node_tag)

        # 2. nodal loads 2d
        nodal_loads = ops.nodeUnbalance(node_tag)
        nodal_loads_idx = np.nonzero(nodal_loads)

        if nodal_loads_idx[0].size:
            for kier in nodal_loads_idx[0]:
                # horizontal or vertical nodal force (load)
                if kier == 0 or kier == 1:
                    if kier == 0:
                        kier2 = np.sign(nodal_loads[kier])
                        if kier2 > 0:
                            pos_or_neg = '+'
                        elif kier2 < 0:
                            pos_or_neg = '-'
                        dx = sfac*np.sign(kier2)
                        dy = 0.

                    elif kier == 1:
                        kier2 = np.sign(nodal_loads[kier])
                        if kier2 > 0:
                            pos_or_neg = '+'
                        elif kier2 < 0:
                            pos_or_neg = '-'
                        dx = 0.
                        dy = sfac*np.sign(kier2)

                    ax.arrow(nd_crd[0], nd_crd[1],
                             dx, dy,
                             # width = 0.01,
                             # head_starts_at_zero=True,  # default False
                             # overhang=0.2,
                             lw=3,
                             head_width=0.1*sfac, head_length=0.2*sfac,
                             fc='b', ec='b',
                             length_includes_head=True, shape='full')
                    ax.text(nd_crd[0]+dx, nd_crd[1]+dy, f' {abs(nodal_loads[kier]):.5g}', color='b')

                # concentrated bending moment
                elif kier == 2:
                    kier2 = np.sign(nodal_loads[kier])
                    if kier2 > 0:
                        pos_or_neg = 'anti-clockwise'
                        # marker_type=r'$\circlearrowleft$'
                        marker_type=r'$\curvearrowleft$'
                    elif kier2 < 0:
                        pos_or_neg = 'clockwise'
                        # marker_type=r'$\circlearrowright$'
                        marker_type=r'$\curvearrowright$'

                    ax.plot(nd_crd[0], nd_crd[1], marker=marker_type, markersize=35, color='b')
                    ax.text(nd_crd[0], nd_crd[1], f'{abs(nodal_loads[kier]):.5g}', color='b', va='bottom', ha='left')

    ax.axis('equal')
    ax.grid(False)

    return ax


def get_nodal_loads_from_ops_domain():
    load_at_nodes = {}
    ibeg, iend = 0, 0

    node_load_tags_all_patterns = ops.getNodeLoadTags()
    node_load_data_all_patterns = ops.getNodeLoadData()

    for node_load_tag in node_load_tags_all_patterns:
        load_at_nodes[node_load_tag] = []

    for node_load_tag in node_load_tags_all_patterns:
        ndf = ops.getNDF(node_load_tag)[0]
        iend = ibeg+ndf
        node_load_data = node_load_data_all_patterns[ibeg:iend]
        ibeg = iend
        load_at_nodes[node_load_tag].append(node_load_data)

    return load_at_nodes


def get_Ew_data_from_ops_domain():

    Ew = {}
    ibeg = 0
    iend = 0
    ele_load_tags_all_patterns = ops.getEleLoadTags()
    ele_load_types_all_patterns = ops.getEleLoadClassTags()
    ele_load_data_all_patterns = ops.getEleLoadData()

    for ele_load_tag in ele_load_tags_all_patterns:
        Ew[ele_load_tag] = []

    # iterate over ele_load_tags/classtags
    for ele_load_type, ele_load_tag in zip(ele_load_types_all_patterns,
                                           ele_load_tags_all_patterns):

        # -beamUniform
        if ele_load_type == LoadTag.Beam2dUniformLoad:
            iend = ibeg+LoadTag.Beam2dUniformLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            wy, wx = ele_load_data[0], ele_load_data[1]
            Ew[ele_load_tag].append(['-beamUniform', wy, wx])
        # -beamUniform partial
        elif ele_load_type == LoadTag.Beam2dPartialUniformLoad:
            iend = ibeg+LoadTag.Beam2dPartialUniformLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            wta, wtb, waa, wab, aL, bL = ele_load_data
            Ew[ele_load_tag].append(['-beamUniform', wta, waa, aL, bL, wtb, wab])
        # -beamPoint
        elif ele_load_type == LoadTag.Beam2dPointLoad:
            iend = ibeg+LoadTag.Beam2dPointLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            Pt, Pa, aL = ele_load_data
            Ew[ele_load_tag].append(['-beamPoint', Pt, aL, Pa])
        else:
            print(f'\nWarning! ele_load_type:\n{ele_load_type} - Unknown element load Error')


    return Ew


def get_Ew_data_from_ops_domain_3d():

    Ew = {}
    ibeg = 0
    iend = 0
    ele_load_tags_all_patterns = ops.getEleLoadTags()
    ele_load_types_all_patterns = ops.getEleLoadClassTags()
    ele_load_data_all_patterns = ops.getEleLoadData()

    for ele_load_tag in ele_load_tags_all_patterns:
        Ew[ele_load_tag] = []

    # iterate over ele_load_tags/classtags
    for ele_load_type, ele_load_tag in zip(ele_load_types_all_patterns,
                                           ele_load_tags_all_patterns):

        # -beamUniform
        if ele_load_type == LoadTag.Beam3dUniformLoad:
            iend = ibeg+LoadTag.Beam3dUniformLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            wy, wz, wx = ele_load_data[0], ele_load_data[1], ele_load_data[2]
            Ew[ele_load_tag].append(['-beamUniform', wy, wz, wx])
        # -beamUniform partial
        elif ele_load_type == LoadTag.Beam3dPartialUniformLoad:
            iend = ibeg+LoadTag.Beam3dPartialUniformLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            wy, wz, wa, aL, bL = ele_load_data
            Ew[ele_load_tag].append(['-beamUniform', wy, wz, wa, aL, bL])
        # -beamPoint
        elif ele_load_type == LoadTag.Beam3dPointLoad:
            iend = ibeg+LoadTag.Beam3dPointLoad_ndata
            ele_load_data = ele_load_data_all_patterns[ibeg:iend]
            ibeg = iend
            Py, Pz, Px, aL = ele_load_data
            Ew[ele_load_tag].append(['-beamPoint', Py, Pz, aL, Px])
        else:
            print(f'\nWarning! ele_load_type:\n{ele_load_type} - Unknown element load Error')


    return Ew


def plot_extruded_shapes_3d(ele_shapes, az_el=az_el,
                            fig_wi_he=False,
                            fig_lbrt=False, ax=False):
    """Plot an extruded 3d model based on cross-section dimenions.

    Three arrows present local section axes: green - local x-axis,
    red - local z-axis, blue - local y-axis.

    Args:
        ele_shapes (dict): keys are ele_tags and values are lists of:
            shape_type (str): 'rect' - rectangular shape, 'I' - double T shape
            and shape_args (list): list of floats, which necessary section
            dimensions.
            For 'rect' the list is [b d]; width and depth,
            for 'I' shape - [bf d tw tf]; flange width, section depth, web
            and flange thicknesses
            Example: ele_shapes = {1: ['rect', [b, d]],
            2: ['I', [bf, d, tw, tf]]}

        az_el (tuple): azimuth and elevation

        fig_wi_he: figure width and height in centimeters

        fig_lbrt: figure left, bottom, right, top boundaries

    Usage:
        ::

            ele_shapes = {1: ['circ', [b]],
                          2: ['rect', [b, h]],
                          3: ['I', [b, h, b/10., h/6.]]}
            opsv.plot_extruded_shapes_3d(ele_shapes)

    Notes:

    - For now only rectangular, circular and double T sections are supported.

    - This function can be a source of inconsistency because OpenSees lacks
      functions to return section dimensions. A workaround is to have own
      Python helper functions to reuse data specified once
    """

    ele_tags = ops.getEleTags()

    azim, elev = az_el

    if not ax:
        if fig_wi_he:
            fig_wi, fig_he = fig_wi_he

            fig = plt.figure(figsize=(fig_wi / 2.54, fig_he / 2.54))
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

        ele_node_tags = ops.eleNodes(ele_tag)
        nen = 2
        ecrd = np.zeros((nen, 3))
        for i, ele_node_tag in enumerate(ele_node_tags):
            ecrd[i, :] = ops.nodeCoord(ele_node_tag)

        # mesh outline
        ax.plot(ecrd[:, 0], ecrd[:, 1], ecrd[:, 2], 'k--', solid_capstyle='round',
                solid_joinstyle='round', dash_capstyle='butt',
                dash_joinstyle='round')

        # eo = Eo[i, :]
        xloc = ops.eleResponse(ele_tag, 'xlocal')
        yloc = ops.eleResponse(ele_tag, 'ylocal')
        zloc = ops.eleResponse(ele_tag, 'zlocal')
        g = np.vstack((xloc, yloc, zloc))

        G, L = rot_transf_3d(ecrd, g)

        # by default empty
        shape_type, shape_args = ['circ', [0.]]
        if ele_tag in ele_shapes:
            shape_type, shape_args = ele_shapes[ele_tag]

        if shape_type == 'rect' or shape_type == 'I':
            if shape_type == 'rect':
                verts = _plot_extruded_shapes_3d_rect(ecrd, g,
                                                      shape_args)

            elif shape_type == 'I':
                verts = _plot_extruded_shapes_3d_double_T(ecrd, g,
                                                          shape_args)

            # plot 3d sides
            ax.add_collection3d(Poly3DCollection(verts, linewidths=0.6,
                                                 edgecolors='k', alpha=.25))

        elif shape_type == 'circ':
            X, Y, Z = _plot_extruded_shapes_3d_circ(ecrd, g,
                                                    shape_args)
            ax.plot_surface(X, Y, Z, linewidths=0.4, edgecolors='k',
                            alpha=.25)

        # common for all members
        Xm, Ym, Zm = sum(ecrd[:, 0]) / 2, sum(ecrd[:, 1]) / 2, sum(ecrd[:, 2]) / 2

        alen = 0.1*L

        # plot local axis directional vectors: workaround quiver = arrow
        ax.quiver(Xm, Ym, Zm, g[0, 0], g[0, 1], g[0, 2], color='g',
                  lw=2, length=alen, alpha=.8, normalize=True)
        ax.quiver(Xm, Ym, Zm, g[1, 0], g[1, 1], g[1, 2], color='b',
                  lw=2, length=alen, alpha=.8, normalize=True)
        ax.quiver(Xm, Ym, Zm, g[2, 0], g[2, 1], g[2, 2], color='r',
                  lw=2, length=alen, alpha=.8, normalize=True)

    ax.set_box_aspect((np.ptp(ax.get_xlim3d()),
                       np.ptp(ax.get_ylim3d()),
                       np.ptp(ax.get_zlim3d())))


def _plot_extruded_shapes_3d_double_T(ecrd, g, shape_args):
    bf, d, tw, tf = shape_args

    za, zb = tw / 2, bf / 2
    ya, yb = d / 2 - tf, d / 2

    g10a, g11a, g12a = g[1, 0] * ya, g[1, 1] * ya, g[1, 2] * ya
    g10b, g11b, g12b = g[1, 0] * yb, g[1, 1] * yb, g[1, 2] * yb

    g20a, g21a, g22a = g[2, 0] * za, g[2, 1] * za, g[2, 2] * za
    g20b, g21b, g22b = g[2, 0] * zb, g[2, 1] * zb, g[2, 2] * zb

    pts = np.zeros((24, 3))
    # beg node (I) cross-section vertices crds, counter-clockwise order
    pts[0] = [ecrd[0, 0] - g10b - g20b, ecrd[0, 1] - g11b - g21b, ecrd[0, 2] - g12b - g22b]
    pts[1] = [ecrd[0, 0] - g10a - g20b, ecrd[0, 1] - g11a - g21b, ecrd[0, 2] - g12a - g22b]
    pts[2] = [ecrd[0, 0] - g10a - g20a, ecrd[0, 1] - g11a - g21a, ecrd[0, 2] - g12a - g22a]

    pts[3] = [ecrd[0, 0] + g10a - g20a, ecrd[0, 1] + g11a - g21a, ecrd[0, 2] + g12a - g22a]
    pts[4] = [ecrd[0, 0] + g10a - g20b, ecrd[0, 1] + g11a - g21b, ecrd[0, 2] + g12a - g22b]
    pts[5] = [ecrd[0, 0] + g10b - g20b, ecrd[0, 1] + g11b - g21b, ecrd[0, 2] + g12b - g22b]

    pts[6] = [ecrd[0, 0] + g10b + g20b, ecrd[0, 1] + g11b + g21b, ecrd[0, 2] + g12b + g22b]
    pts[7] = [ecrd[0, 0] + g10a + g20b, ecrd[0, 1] + g11a + g21b, ecrd[0, 2] + g12a + g22b]
    pts[8] = [ecrd[0, 0] + g10a + g20a, ecrd[0, 1] + g11a + g21a, ecrd[0, 2] + g12a + g22a]

    pts[9] = [ecrd[0, 0] - g10a + g20a, ecrd[0, 1] - g11a + g21a, ecrd[0, 2] - g12a + g22a]
    pts[10] = [ecrd[0, 0] - g10a + g20b, ecrd[0, 1] - g11a + g21b, ecrd[0, 2] - g12a + g22b]
    pts[11] = [ecrd[0, 0] - g10b + g20b, ecrd[0, 1] - g11b + g21b, ecrd[0, 2] - g12b + g22b]

    # end node (J) cross-section vertices
    pts[12] = [ecrd[1, 0] - g10b - g20b, ecrd[1, 1] - g11b - g21b, ecrd[1, 2] - g12b - g22b]
    pts[13] = [ecrd[1, 0] - g10a - g20b, ecrd[1, 1] - g11a - g21b, ecrd[1, 2] - g12a - g22b]
    pts[14] = [ecrd[1, 0] - g10a - g20a, ecrd[1, 1] - g11a - g21a, ecrd[1, 2] - g12a - g22a]

    pts[15] = [ecrd[1, 0] + g10a - g20a, ecrd[1, 1] + g11a - g21a, ecrd[1, 2] + g12a - g22a]
    pts[16] = [ecrd[1, 0] + g10a - g20b, ecrd[1, 1] + g11a - g21b, ecrd[1, 2] + g12a - g22b]
    pts[17] = [ecrd[1, 0] + g10b - g20b, ecrd[1, 1] + g11b - g21b, ecrd[1, 2] + g12b - g22b]

    pts[18] = [ecrd[1, 0] + g10b + g20b, ecrd[1, 1] + g11b + g21b, ecrd[1, 2] + g12b + g22b]
    pts[19] = [ecrd[1, 0] + g10a + g20b, ecrd[1, 1] + g11a + g21b, ecrd[1, 2] + g12a + g22b]
    pts[20] = [ecrd[1, 0] + g10a + g20a, ecrd[1, 1] + g11a + g21a, ecrd[1, 2] + g12a + g22a]

    pts[21] = [ecrd[1, 0] - g10a + g20a, ecrd[1, 1] - g11a + g21a, ecrd[1, 2] - g12a + g22a]
    pts[22] = [ecrd[1, 0] - g10a + g20b, ecrd[1, 1] - g11a + g21b, ecrd[1, 2] - g12a + g22b]
    pts[23] = [ecrd[1, 0] - g10b + g20b, ecrd[1, 1] - g11b + g21b, ecrd[1, 2] - g12b + g22b]

    # list of sides
    verts = [[pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6],
              pts[7], pts[8], pts[9], pts[10], pts[11]],  # beg
             [pts[12], pts[13], pts[14], pts[15], pts[16], pts[17],
              pts[18], pts[19], pts[20], pts[21], pts[22], pts[23]],  # end
             [pts[0], pts[12], pts[23], pts[11]],  # bot 1
             [pts[5], pts[17], pts[18], pts[6]],  # top 2
             [pts[9], pts[8], pts[20], pts[21]],  # 3
             [pts[2], pts[3], pts[15], pts[14]],  # 4
             [pts[11], pts[10], pts[22], pts[23]],  # 5
             [pts[0], pts[1], pts[13], pts[12]],  # 6
             [pts[7], pts[6], pts[18], pts[19]],  # 7
             [pts[4], pts[5], pts[17], pts[16]],  # 8
             [pts[10], pts[9], pts[21], pts[22]],  # 9
             [pts[2], pts[1], pts[13], pts[14]],  # 10
             [pts[7], pts[8], pts[20], pts[19]],  # 11
             [pts[3], pts[4], pts[16], pts[15]]]  # 12

    return verts


def _plot_extruded_shapes_3d_rect(ecrd, g, shape_args):
    b, d = shape_args
    b2, d2 = b / 2, d / 2

    g10, g11, g12 = g[1, 0] * d2, g[1, 1] * d2, g[1, 2] * d2
    g20, g21, g22 = g[2, 0] * b2, g[2, 1] * b2, g[2, 2] * b2

    pts = np.zeros((8, 3))
    # beg node cross-section vertices
    # collected i-beg, j-end node coordinates, counter-clockwise order
    pts[0] = [ecrd[0, 0] - g10 - g20, ecrd[0, 1] - g11 - g21, ecrd[0, 2] - g12 - g22]
    pts[1] = [ecrd[0, 0] + g10 - g20, ecrd[0, 1] + g11 - g21, ecrd[0, 2] + g12 - g22]
    pts[2] = [ecrd[0, 0] + g10 + g20, ecrd[0, 1] + g11 + g21, ecrd[0, 2] + g12 + g22]
    pts[3] = [ecrd[0, 0] - g10 + g20, ecrd[0, 1] - g11 + g21, ecrd[0, 2] - g12 + g22]
    # end node cross-section vertices
    pts[4] = [ecrd[1, 0] - g10 - g20, ecrd[1, 1] - g11 - g21, ecrd[1, 2] - g12 - g22]
    pts[5] = [ecrd[1, 0] + g10 - g20, ecrd[1, 1] + g11 - g21, ecrd[1, 2] + g12 - g22]
    pts[6] = [ecrd[1, 0] + g10 + g20, ecrd[1, 1] + g11 + g21, ecrd[1, 2] + g12 + g22]
    pts[7] = [ecrd[1, 0] - g10 + g20, ecrd[1, 1] - g11 + g21, ecrd[1, 2] - g12 + g22]

    # list of 4-node sides
    verts = [[pts[0], pts[1], pts[2], pts[3]],  # beg
             [pts[4], pts[5], pts[6], pts[7]],  # end
             [pts[0], pts[4], pts[5], pts[1]],  # bottom
             [pts[3], pts[7], pts[6], pts[2]],  # top
             [pts[0], pts[4], pts[7], pts[3]],  # front
             [pts[1], pts[5], pts[6], pts[2]]]  # back

    return verts


def _plot_extruded_shapes_3d_circ(ecrd, g, shape_args):
    d = shape_args[0]
    r = d / 2

    Lxyz = ecrd[1, :] - ecrd[0, :]
    L = np.sqrt(Lxyz @ Lxyz)

    xl = np.linspace(0, L, 3)  # subdivde member length
    alpha = np.linspace(0, 2 * np.pi, 10)  # subdivide circle

    xl, alpha = np.meshgrid(xl, alpha)
    X = (ecrd[0, 0] + g[0, 0] * xl + r * np.sin(alpha) * g[1, 0]
         + r * np.cos(alpha) * g[2, 0])
    Y = (ecrd[0, 1] + g[0, 1] * xl + r * np.sin(alpha) * g[1, 1]
         + r * np.cos(alpha) * g[2, 1])
    Z = (ecrd[0, 2] + g[0, 2] * xl + r * np.sin(alpha) * g[1, 2]
         + r * np.cos(alpha) * g[2, 2])

    return X, Y, Z


def bar_length(ecrd):
    Lxyz = ecrd[1, :] - ecrd[0, :]
    L = np.sqrt(Lxyz @ Lxyz)

    return L


def rot_transf_2d(ecrd):

    Lxy = ecrd[1, :] - ecrd[0, :]
    L = np.sqrt(Lxy @ Lxy)
    cosa, cosb = Lxy / L
    G = np.array([[cosa, cosb, 0., 0., 0., 0.],
                  [-cosb, cosa, 0., 0., 0., 0.],
                  [0., 0., 1., 0., 0., 0.],
                  [0., 0., 0., cosa, cosb, 0.],
                  [0., 0., 0., -cosb, cosa, 0.],
                  [0., 0., 0., 0., 0., 1.]])

    return G, L, cosa, cosb


def rot_transf_3d(ecrd, g):

    Lxyz = ecrd[1, :] - ecrd[0, :]
    L = np.sqrt(Lxyz @ Lxyz)

    z = np.zeros((3, 3))
    G = np.block([[g, z, z, z],
                  [z, g, z, z],
                  [z, z, g, z],
                  [z, z, z, g]])

    return G, L
