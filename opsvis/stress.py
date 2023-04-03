import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from .settings import *


def stress_2d_ele_tags_only(ele_tags):

    Stress2dEleClasstags = set([EleClassTag.tri3n,
                                EleClassTag.tri6n,
                                EleClassTag.quad4n,
                                EleClassTag.quad8n,
                                EleClassTag.quad9n])

    ele_classtags = ops.getEleClassTags()

    idx = [i for i, e in enumerate(ele_classtags) if e in Stress2dEleClasstags]

    ele_tags_2d_tris_quads_only = [ele_tags[i] for i in idx]

    return ele_tags_2d_tris_quads_only


def sig_out_per_node(how_many='all'):
    """Return a 2d numpy array of stress components per OpenSees node.

    Three first stress components (sxx, syy, sxy) are calculated and
    extracted from OpenSees, while the rest svm (Huber-Mieses-Hencky),
    two principal stresses (s1, s2) and directional angle are calculated
    as postprocessed quantities.

    Args:
        how_many (str): supported options are: 'all' - all components,
            'sxx', 'syy', 'sxy', 'svm' (or 'vmis'), 's1', 's2', 'angle'.

    Returns:
        sig_out (ndarray): a 2d array of stress components per node with
            the following components: sxx, syy, sxy, svm, s1, s2, angle.
            Size (n_nodes x 7).

    Examples:
        sig_out = opsv.sig_out_per_node()

    Notes:
       s1, s2: principal stresses
       angle: angle of the principal stress s1
    """
    ele_tags = ops.getEleTags()
    node_tags = ops.getNodeTags()
    n_nodes = len(node_tags)

    ele_classtag = ops.getEleClassTags(ele_tags[0])[0]

    # initialize helper arrays
    sig_out = np.zeros((n_nodes, 7))

    nodes_tag_count = np.zeros((n_nodes, 2), dtype=int)
    nodes_tag_count[:, 0] = node_tags

    nen = np.shape(ops.eleNodes(ele_tags[0]))[0]

    for i, ele_tag in enumerate(ele_tags):
        ele_node_tags = ops.eleNodes(ele_tag)

        tmp_list = [0]*nen
        for j, ele_node_tag in enumerate(ele_node_tags):
            tmp_list[j] = node_tags.index(ele_node_tag)

        nodes_tag_count[tmp_list, 1] += 1

        if ele_classtag == EleClassTag.SSPquad:
            sig_ip_el = ops.eleResponse(ele_tag, 'stress')
            # sigM_nd = sigM_ip
            sigM_np = np.tile(sig_ip_el, (nen, 1))

        else:
            sig_nd_el = ops.eleResponse(ele_tag, 'stressAtNodes')
            sigM_nd = np.reshape(sig_nd_el, (-1, 3))

        # sxx yy xy components
        sig_out[tmp_list, 0] += sigM_nd[:nen, 0]
        sig_out[tmp_list, 1] += sigM_nd[:nen, 1]
        sig_out[tmp_list, 2] += sigM_nd[:nen, 2]

    indxs, = np.where(nodes_tag_count[:, 1] > 1)

    # n_indxs < n_nodes: e.g. 21<25 (bous), 2<6 (2el) etc.
    n_indxs = np.shape(indxs)[0]

    # divide summed stresses by the number of common nodes
    sig_out[indxs, :] = \
        sig_out[indxs, :]/nodes_tag_count[indxs, 1].reshape(n_indxs, 1)

    if how_many == 'all' or how_many == 'svm' or how_many == 'vmis':
        # warning reshape from (pts,ncomp) to (ncomp,pts)
        vm_out = vm_stress(np.transpose(sig_out[:, :3]))
        sig_out[:, 3] = vm_out

    if (how_many == 'all' or how_many == 's1' or how_many == 's2' or how_many == 'angle'):  # noqa: E501
        princ_sig_out = princ_stress(np.transpose(sig_out[:, :3]))

        sig_out[:, 4:7] = np.transpose(princ_sig_out)

    print('-- WARNING!!! full sig_out matrix calculated here --')

    return sig_out


def sig_component_per_node(stress_str):
    """Return a 2d numpy array of stress components per OpenSees node.

    Three first stress components (sxx, syy, sxy) are calculated and
    extracted from OpenSees, while the rest svm (Huber-Mieses-Hencky),
    two principal stresses (s1, s2) and directional angle are calculated
    as postprocessed quantities.

    Args:
        how_many (str): supported options are: 'all' - all components,
            'sxx', 'syy', 'sxy', 'svm' (or 'vmis'), 's1', 's2', 'angle'.

    Returns:
        sig_out (ndarray): a 2d array of stress components per node with
            the following components: sxx, syy, sxy, svm, s1, s2, angle.
            Size (n_nodes x 7).

    Examples:
        sig_out = opsv.sig_out_per_node()

    Notes:
       s1, s2: principal stresses
       angle: angle of the principal stress s1
    """

    ele_tags_all = ops.getEleTags()
    ele_tags = stress_2d_ele_tags_only(ele_tags_all)

    ele_classtag = ops.getEleClassTags(ele_tags[0])[0]

    node_tags = ops.getNodeTags()
    n_nodes = len(node_tags)

    # initialize helper arrays
    sig_out = np.zeros((n_nodes, 4))

    nodes_tag_count = np.zeros((n_nodes, 2), dtype=int)
    nodes_tag_count[:, 0] = node_tags

    nen = np.shape(ops.eleNodes(ele_tags[0]))[0]

    for i, ele_tag in enumerate(ele_tags):
        ele_node_tags = ops.eleNodes(ele_tag)

        tmp_list = [0]*nen
        for j, ele_node_tag in enumerate(ele_node_tags):
            tmp_list[j] = node_tags.index(ele_node_tag)

        nodes_tag_count[tmp_list, 1] += 1

        if ele_classtag == EleClassTag.SSPquad:
            sig_ip_el = ops.eleResponse(ele_tag, 'stress')
            # sigM_nd = sigM_ip
            sigM_np = np.tile(sig_ip_el, (nen, 1))

        else:
            sig_nd_el = ops.eleResponse(ele_tag, 'stressAtNodes')
            sigM_nd = np.reshape(sig_nd_el, (-1, 3))

        # sxx yy xy components
        sig_out[tmp_list, 0] += sigM_nd[:nen, 0]
        sig_out[tmp_list, 1] += sigM_nd[:nen, 1]
        sig_out[tmp_list, 2] += sigM_nd[:nen, 2]

    indxs, = np.where(nodes_tag_count[:, 1] > 1)

    # n_indxs < n_nodes: e.g. 21<25 (bous), 2<6 (2el) etc.
    n_indxs = np.shape(indxs)[0]

    # divide summed stresses by the number of common nodes
    sig_out[indxs, :] = \
        sig_out[indxs, :]/nodes_tag_count[indxs, 1].reshape(n_indxs, 1)

    if stress_str == 'sxx':
        sig_out_vec = sig_out[:, 0]
    elif stress_str == 'syy':
        sig_out_vec = sig_out[:, 1]
    elif stress_str == 'sxy':
        sig_out_vec = sig_out[:, 2]
    elif stress_str == 'svm' or stress_str == 'vmis':
        # warning reshape from (pts,ncomp) to (ncomp,pts)
        sig_out_vec = vm_stress(np.transpose(sig_out[:, :3]))
    elif (stress_str == 's1' or stress_str == 's2' or stress_str == 'angle'):
        princ_sig_out = princ_stress(np.transpose(sig_out[:, :3]))
        if stress_str == 's1':
            # sig_out_vec = np.transpose(princ_sig_out)[:, 0]
            sig_out_vec = princ_sig_out[0, :]
        elif stress_str == 's2':
            sig_out_vec = princ_sig_out[1, :]
        elif stress_str == 'angle':
            sig_out_vec = princ_sig_out[2, :]

    return sig_out_vec


def princ_stress(sig):
    """Return a tuple (s1, s2, angle): principal stresses (plane stress) and angle
    Args:
        sig (ndarray): input array of stresses at nodes: sxx, syy, sxy (tau)

    Returns:
        out (ndarray): 1st row is first principal stress s1, 2nd row is second
           principal stress s2, 3rd row is the angle of s1
    """
    sx, sy, tau = sig[0], sig[1], sig[2]

    ds = (sx-sy)/2
    R = np.sqrt(ds**2 + tau**2)

    s1 = (sx+sy)/2. + R
    s2 = (sx+sy)/2. - R
    angle = np.arctan2(tau, ds)/2

    out = np.vstack((s1, s2, angle))

    return out


def vm_stress(sig):
    n_sig_comp, n_pts = np.shape(sig)
    if n_sig_comp > 3:
        x, y, z, xy, xz, yz = sig
    else:
        x, y, xy = sig
        z, xz, yz = 0., 0., 0.

    _a = 0.5*((x-y)**2 + (y-z)**2 + (z-x)**2 + 6.*(xy**2 + xz**2 + yz**2))
    return np.sqrt(_a)


def quads_to_4tris(quads_conn, nds_crd, nds_val):
    """
    Get triangles connectivity, coordinates and new values at quad centroids.

    Args:
        quads_conn (ndarray):

        nds_crd (ndarray):

        nds_val (ndarray):

    Returns:
        tris_conn, nds_c_crd, nds_c_val (tuple):

    Notes:
        Triangles connectivity array is based on
        quadrilaterals connectivity.
        Each quad is split into four triangles.
        New nodes are created at the quad centroid.

    See also:
        function: quads_to_8tris_9n, quads_to_8tris_8n
    """
    n_quads, _ = quads_conn.shape
    n_nds, _ = nds_crd.shape

    # coordinates and values at quad centroids _c_
    nds_c_crd = np.zeros((n_quads, 2))
    nds_c_val = np.zeros(n_quads)

    tris_conn = np.zeros((4*n_quads, 3), dtype=int)

    for i, quad_conn in enumerate(quads_conn):
        j = 4*i
        n0, n1, n2, n3 = quad_conn

        # quad centroids
        nds_c_crd[i] = np.array([np.sum(nds_crd[[n0, n1, n2, n3], 0])/4.,
                                 np.sum(nds_crd[[n0, n1, n2, n3], 1])/4.])
        nds_c_val[i] = np.sum(nds_val[[n0, n1, n2, n3]])/4.

        # triangles connectivity
        tris_conn[j] = np.array([n0, n1, n_nds+i])
        tris_conn[j+1] = np.array([n1, n2, n_nds+i])
        tris_conn[j+2] = np.array([n2, n3, n_nds+i])
        tris_conn[j+3] = np.array([n3, n0, n_nds+i])

    return tris_conn, nds_c_crd, nds_c_val


def bricks_to_24tris(bricks_conn, nds_crd, nds_val, disps=None):
    """
    Get triangles connectivity, coordinates and new vals at brick face centroids.

    Args:
        bricks_conn (ndarray):

        nds_crd (ndarray):

        nds_val (ndarray):

    Returns:
        tris_conn, nds_c_crd, nds_c_val (tuple):

    Notes:
        Triangles connectivity array is based on
        stdBricks connectivity.
        Each of 6 brick faces is split into four triangles.
        New nodes are created at the face centroid.

    See also:
        function: bricks_to_8tris_9n, bricks_to_8tris_8n
    """
    n_bricks, _ = bricks_conn.shape
    n_nds, _ = nds_crd.shape

    # coordinates and values at brick centroids _c_
    nds_c_crd = np.zeros((n_bricks*6, 3))
    nds_c_val = np.zeros(n_bricks*6)

    disps_c = None
    if disps is not None:
        disps_c = np.zeros((n_bricks*6, 3))

    tris_conn = np.zeros((24*n_bricks, 3), dtype=int)

    for i, brick_conn in enumerate(bricks_conn):
        j = 24*i
        n0, n1, n2, n3, n4, n5, n6, n7 = brick_conn

        # brick centroids
        nds_c_crd[i*6] = np.array([np.sum(nds_crd[[n0, n1, n5, n4], 0])/4.,
                                   np.sum(nds_crd[[n0, n1, n5, n4], 1])/4.,
                                   np.sum(nds_crd[[n0, n1, n5, n5], 2])/4.])
        nds_c_crd[i*6+1] = np.array([np.sum(nds_crd[[n1, n2, n6, n5], 0])/4.,
                                     np.sum(nds_crd[[n1, n2, n6, n5], 1])/4.,
                                     np.sum(nds_crd[[n1, n2, n6, n5], 2])/4.])
        nds_c_crd[i*6+2] = np.array([np.sum(nds_crd[[n2, n3, n7, n6], 0])/4.,
                                     np.sum(nds_crd[[n2, n3, n7, n6], 1])/4.,
                                     np.sum(nds_crd[[n2, n3, n7, n6], 2])/4.])
        nds_c_crd[i*6+3] = np.array([np.sum(nds_crd[[n3, n0, n4, n7], 0])/4.,
                                     np.sum(nds_crd[[n3, n0, n4, n7], 1])/4.,
                                     np.sum(nds_crd[[n3, n0, n4, n7], 2])/4.])
        nds_c_crd[i*6+4] = np.array([np.sum(nds_crd[[n4, n5, n6, n7], 0])/4.,
                                     np.sum(nds_crd[[n4, n5, n6, n7], 1])/4.,
                                     np.sum(nds_crd[[n4, n5, n6, n7], 2])/4.])
        nds_c_crd[i*6+5] = np.array([np.sum(nds_crd[[n0, n1, n2, n3], 0])/4.,
                                     np.sum(nds_crd[[n0, n1, n2, n3], 1])/4.,
                                     np.sum(nds_crd[[n0, n1, n2, n3], 2])/4.])

        nds_c_val[6*i] = np.sum(nds_val[[n0, n1, n5, n4]])/4.
        nds_c_val[6*i+1] = np.sum(nds_val[[n1, n2, n6, n5]])/4.
        nds_c_val[6*i+2] = np.sum(nds_val[[n2, n3, n7, n6]])/4.
        nds_c_val[6*i+3] = np.sum(nds_val[[n3, n0, n4, n7]])/4.
        nds_c_val[6*i+4] = np.sum(nds_val[[n4, n5, n6, n7]])/4.
        nds_c_val[6*i+5] = np.sum(nds_val[[n0, n1, n2, n3]])/4.

        # triangles connectivity
        tris_conn[j] = np.array([n0, n1, n_nds+i*6])
        tris_conn[j+1] = np.array([n1, n5, n_nds+i*6])
        tris_conn[j+2] = np.array([n5, n4, n_nds+i*6])
        tris_conn[j+3] = np.array([n4, n0, n_nds+i*6])

        tris_conn[j+4] = np.array([n1, n2, n_nds+i*6+1])
        tris_conn[j+5] = np.array([n2, n6, n_nds+i*6+1])
        tris_conn[j+6] = np.array([n6, n5, n_nds+i*6+1])
        tris_conn[j+7] = np.array([n5, n1, n_nds+i*6+1])

        tris_conn[j+8] = np.array([n2, n3, n_nds+i*6+2])
        tris_conn[j+9] = np.array([n3, n7, n_nds+i*6+2])
        tris_conn[j+10] = np.array([n7, n6, n_nds+i*6+2])
        tris_conn[j+11] = np.array([n6, n2, n_nds+i*6+2])

        tris_conn[j+12] = np.array([n3, n0, n_nds+i*6+3])
        tris_conn[j+13] = np.array([n0, n4, n_nds+i*6+3])
        tris_conn[j+14] = np.array([n4, n7, n_nds+i*6+3])
        tris_conn[j+15] = np.array([n7, n3, n_nds+i*6+3])

        tris_conn[j+16] = np.array([n4, n5, n_nds+i*6+4])
        tris_conn[j+17] = np.array([n5, n6, n_nds+i*6+4])
        tris_conn[j+18] = np.array([n6, n7, n_nds+i*6+4])
        tris_conn[j+19] = np.array([n7, n4, n_nds+i*6+4])

        tris_conn[j+20] = np.array([n0, n1, n_nds+i*6+5])
        tris_conn[j+21] = np.array([n1, n2, n_nds+i*6+5])
        tris_conn[j+22] = np.array([n2, n3, n_nds+i*6+5])
        tris_conn[j+23] = np.array([n3, n0, n_nds+i*6+5])

        if disps is not None:
            disps_c[6*i] = np.sum(disps[[n0, n1, n5, n4]], axis=0)/4.
            disps_c[6*i+1] = np.sum(disps[[n1, n2, n6, n5]], axis=0)/4.
            disps_c[6*i+2] = np.sum(disps[[n2, n3, n7, n6]], axis=0)/4.
            disps_c[6*i+3] = np.sum(disps[[n3, n0, n4, n7]], axis=0)/4.
            disps_c[6*i+4] = np.sum(disps[[n4, n5, n6, n7]], axis=0)/4.
            disps_c[6*i+5] = np.sum(disps[[n0, n1, n2, n3]], axis=0)/4.

    return tris_conn, nds_c_crd, nds_c_val, disps_c


# brick20n bricks to tris
def bricks_to_48tris(bricks_conn, nds_crd, nds_val, disps=None):
    """
    Get triangles connectivity, coordinates and new vals at brick face centroids,
    for brick20n

    Args:
        bricks_conn (ndarray):

        nds_crd (ndarray):

        nds_val (ndarray):

    Returns:
        tris_conn, nds_c_crd, nds_c_val (tuple):

    Notes:
        Triangles connectivity array is based on
        stdBricks connectivity.
        Each of 6 brick faces is split into four triangles.
        New nodes are created at the face centroid.

    See also:
        function: bricks_to_24tris
    """
    n_bricks, _ = bricks_conn.shape
    n_nds, _ = nds_crd.shape

    # coordinates and values at brick centroids _c_
    nds_c_crd = np.zeros((n_bricks*6, 3))
    nds_c_val = np.zeros(n_bricks*6)

    disps_c = None
    if disps is not None:
        disps_c = np.zeros((n_bricks*6, 3))

    tris_conn = np.zeros((48*n_bricks, 3), dtype=int)

    for i, brick_conn in enumerate(bricks_conn):
        j = 48*i
        n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, \
            n14, n15, n16, n17, n18, n19 = brick_conn

        # brick centroids
        nds_c_crd[i*6] = np.array([np.sum(nds_crd[[n0, n1, n5, n4], 0])/4.,
                                   np.sum(nds_crd[[n0, n1, n5, n4], 1])/4.,
                                   np.sum(nds_crd[[n0, n1, n5, n5], 2])/4.])
        nds_c_crd[i*6+1] = np.array([np.sum(nds_crd[[n1, n2, n6, n5], 0])/4.,
                                     np.sum(nds_crd[[n1, n2, n6, n5], 1])/4.,
                                     np.sum(nds_crd[[n1, n2, n6, n5], 2])/4.])
        nds_c_crd[i*6+2] = np.array([np.sum(nds_crd[[n2, n3, n7, n6], 0])/4.,
                                     np.sum(nds_crd[[n2, n3, n7, n6], 1])/4.,
                                     np.sum(nds_crd[[n2, n3, n7, n6], 2])/4.])
        nds_c_crd[i*6+3] = np.array([np.sum(nds_crd[[n3, n0, n4, n7], 0])/4.,
                                     np.sum(nds_crd[[n3, n0, n4, n7], 1])/4.,
                                     np.sum(nds_crd[[n3, n0, n4, n7], 2])/4.])
        nds_c_crd[i*6+4] = np.array([np.sum(nds_crd[[n4, n5, n6, n7], 0])/4.,
                                     np.sum(nds_crd[[n4, n5, n6, n7], 1])/4.,
                                     np.sum(nds_crd[[n4, n5, n6, n7], 2])/4.])
        nds_c_crd[i*6+5] = np.array([np.sum(nds_crd[[n0, n1, n2, n3], 0])/4.,
                                     np.sum(nds_crd[[n0, n1, n2, n3], 1])/4.,
                                     np.sum(nds_crd[[n0, n1, n2, n3], 2])/4.])

        nds_c_val[6*i] = np.sum(nds_val[[n0, n1, n5, n4]])/4.
        nds_c_val[6*i+1] = np.sum(nds_val[[n1, n2, n6, n5]])/4.
        nds_c_val[6*i+2] = np.sum(nds_val[[n2, n3, n7, n6]])/4.
        nds_c_val[6*i+3] = np.sum(nds_val[[n3, n0, n4, n7]])/4.
        nds_c_val[6*i+4] = np.sum(nds_val[[n4, n5, n6, n7]])/4.
        nds_c_val[6*i+5] = np.sum(nds_val[[n0, n1, n2, n3]])/4.

        # triangles connectivity
        tris_conn[j] = np.array([n0, n8, n_nds+i*6])
        tris_conn[j+1] = np.array([n8, n1, n_nds+i*6])
        tris_conn[j+2] = np.array([n1, n17, n_nds+i*6])
        tris_conn[j+3] = np.array([n17, n5, n_nds+i*6])
        tris_conn[j+4] = np.array([n5, n12, n_nds+i*6])
        tris_conn[j+5] = np.array([n12, n4, n_nds+i*6])
        tris_conn[j+6] = np.array([n4, n16, n_nds+i*6])
        tris_conn[j+7] = np.array([n16, n0, n_nds+i*6])

        tris_conn[j+8] = np.array([n1, n9, n_nds+i*6+1])
        tris_conn[j+9] = np.array([n9, n2, n_nds+i*6+1])
        tris_conn[j+10] = np.array([n2, n18, n_nds+i*6+1])
        tris_conn[j+11] = np.array([n18, n6, n_nds+i*6+1])
        tris_conn[j+12] = np.array([n6, n13, n_nds+i*6+1])
        tris_conn[j+13] = np.array([n13, n5, n_nds+i*6+1])
        tris_conn[j+14] = np.array([n5, n17, n_nds+i*6+1])
        tris_conn[j+15] = np.array([n17, n1, n_nds+i*6+1])

        tris_conn[j+16] = np.array([n2, n10, n_nds+i*6+2])
        tris_conn[j+17] = np.array([n10, n3, n_nds+i*6+2])
        tris_conn[j+18] = np.array([n3, n19, n_nds+i*6+2])
        tris_conn[j+19] = np.array([n19, n7, n_nds+i*6+2])
        tris_conn[j+20] = np.array([n7, n14, n_nds+i*6+2])
        tris_conn[j+21] = np.array([n14, n6, n_nds+i*6+2])
        tris_conn[j+22] = np.array([n6, n18, n_nds+i*6+2])
        tris_conn[j+23] = np.array([n18, n2, n_nds+i*6+2])

        tris_conn[j+24] = np.array([n3, n11, n_nds+i*6+3])
        tris_conn[j+25] = np.array([n11, n0, n_nds+i*6+3])
        tris_conn[j+26] = np.array([n0, n16, n_nds+i*6+3])
        tris_conn[j+27] = np.array([n16, n4, n_nds+i*6+3])
        tris_conn[j+28] = np.array([n4, n15, n_nds+i*6+3])
        tris_conn[j+29] = np.array([n15, n7, n_nds+i*6+3])
        tris_conn[j+30] = np.array([n7, n19, n_nds+i*6+3])
        tris_conn[j+31] = np.array([n19, n3, n_nds+i*6+3])

        tris_conn[j+32] = np.array([n4, n12, n_nds+i*6+4])
        tris_conn[j+33] = np.array([n12, n5, n_nds+i*6+4])
        tris_conn[j+34] = np.array([n5, n13, n_nds+i*6+4])
        tris_conn[j+35] = np.array([n13, n6, n_nds+i*6+4])
        tris_conn[j+36] = np.array([n6, n14, n_nds+i*6+4])
        tris_conn[j+37] = np.array([n14, n7, n_nds+i*6+4])
        tris_conn[j+38] = np.array([n7, n15, n_nds+i*6+4])
        tris_conn[j+39] = np.array([n15, n4, n_nds+i*6+4])

        tris_conn[j+40] = np.array([n0, n8, n_nds+i*6+5])
        tris_conn[j+41] = np.array([n8, n1, n_nds+i*6+5])
        tris_conn[j+42] = np.array([n1, n9, n_nds+i*6+5])
        tris_conn[j+43] = np.array([n9, n2, n_nds+i*6+5])
        tris_conn[j+44] = np.array([n2, n10, n_nds+i*6+5])
        tris_conn[j+45] = np.array([n10, n3, n_nds+i*6+5])
        tris_conn[j+46] = np.array([n3, n11, n_nds+i*6+5])
        tris_conn[j+47] = np.array([n11, n0, n_nds+i*6+5])

        if disps is not None:
            disps_c[6*i] = np.sum(disps[[n0, n1, n5, n4,
                                         n8, n17, n12, n16]], axis=0)/4.
            disps_c[6*i+1] = np.sum(disps[[n1, n2, n6, n5,
                                           n9, n18, n13, n17]], axis=0)/4.
            disps_c[6*i+2] = np.sum(disps[[n2, n3, n7, n6,
                                           n10, n19, n14, n18]], axis=0)/4.
            disps_c[6*i+3] = np.sum(disps[[n3, n0, n4, n7,
                                           n11, n16, n15, n19]], axis=0)/4.
            disps_c[6*i+4] = np.sum(disps[[n4, n5, n6, n7,
                                           n12, n13, n14, n15]], axis=0)/4.
            disps_c[6*i+5] = np.sum(disps[[n0, n1, n2, n3,
                                           n8, n9, n10, n11]], axis=0)/4.

    return tris_conn, nds_c_crd, nds_c_val, disps_c


def tetra4n_to_4tris(tetras4n_conn):
    """Get triangles connectivity.

    Four-node tetrahedron is subdivided into four triangles

    Args:
        tetra4n_conn (ndarray):

    Returns:
        tris_conn_subdiv (ndarray):
    """
    n_tetras, _ = tetras4n_conn.shape
    # n_nds, _ = nds_crd.shape

    tris_conn_subdiv = np.zeros((4*n_tetras, 3), dtype=int)

    for i, tetra4n_conn in enumerate(tetras4n_conn):
        j = 4*i
        n0, n1, n2, n3 = tetra4n_conn

        # triangles connectivity
        tris_conn_subdiv[j] = np.array([n0, n1, n2])
        tris_conn_subdiv[j+1] = np.array([n0, n3, n1])
        tris_conn_subdiv[j+2] = np.array([n0, n2, n3])
        tris_conn_subdiv[j+3] = np.array([n1, n2, n3])

    return tris_conn_subdiv


def tetra10n_to_16tris(tetras10n_conn):
    """Get triangles connectivity.

    Six-node triangle is subdivided into four triangles

    Args:
        tetra10n_conn (ndarray):

    Returns:
        tris_conn_subdiv (ndarray):
    """
    n_tetras, _ = tetras10n_conn.shape
    # n_nds, _ = nds_crd.shape

    tris_conn_subdiv = np.zeros((16*n_tetras, 3), dtype=int)

    for i, tetra10n_conn in enumerate(tetras10n_conn):
        j = 16*i
        n0, n1, n2, n3, n4, n5, n6, n7, n8, n9 = tetra10n_conn

        # triangles connectivity
        tris_conn_subdiv[j] = np.array([n0, n4, n6])
        tris_conn_subdiv[j+1] = np.array([n1, n5, n4])
        tris_conn_subdiv[j+2] = np.array([n4, n5, n6])
        tris_conn_subdiv[j+3] = np.array([n2, n6, n5])

        tris_conn_subdiv[j+4] = np.array([n0, n7, n4])
        tris_conn_subdiv[j+5] = np.array([n4, n7, n8])
        tris_conn_subdiv[j+6] = np.array([n1, n4, n8])
        tris_conn_subdiv[j+7] = np.array([n3, n8, n7])

        tris_conn_subdiv[j+8] = np.array([n0, n6, n7])
        tris_conn_subdiv[j+9] = np.array([n6, n9, n7])
        tris_conn_subdiv[j+10] = np.array([n2, n9, n6])
        tris_conn_subdiv[j+11] = np.array([n3, n7, n9])

        tris_conn_subdiv[j+12] = np.array([n1, n5, n8])
        tris_conn_subdiv[j+13] = np.array([n5, n9, n8])
        tris_conn_subdiv[j+14] = np.array([n2, n9, n5])
        tris_conn_subdiv[j+15] = np.array([n3, n8, n9])

    return tris_conn_subdiv


def tris6n_to_4tris(tris_conn):
    """Get triangles connectivity.

    Six-node triangle is subdivided into four triangles

    Args:
        tris_conn (ndarray):

    Returns:
        tris_conn_subdiv (ndarray):
    """
    n_tris, _ = tris_conn.shape
    # n_nds, _ = nds_crd.shape

    tris_conn_subdiv = np.zeros((4*n_tris, 3), dtype=int)

    for i, tri_conn in enumerate(tris_conn):
        j = 4*i
        n0, n1, n2, n3, n4, n5 = tri_conn

        # triangles connectivity
        tris_conn_subdiv[j] = np.array([n0, n3, n5])
        tris_conn_subdiv[j+1] = np.array([n3, n1, n4])
        tris_conn_subdiv[j+2] = np.array([n3, n4, n5])
        tris_conn_subdiv[j+3] = np.array([n5, n4, n2])

    return tris_conn_subdiv


# def plot_mesh_2d(nds_crd, eles_conn, lw=0.4, ec='k', ax=False):
def plot_mesh_2d(nds_crd, eles_conn, lw=0.4, ec='k'):
    """
    Plot 2d mesh (quads or triangles) outline.
    """
    nen = np.shape(eles_conn)[1]
    if nen == 3 or nen == 4:
        for ele_conn in eles_conn:
            x = nds_crd[ele_conn, 0]
            y = nds_crd[ele_conn, 1]
            plt.fill(x, y, edgecolor=ec, lw=lw, fill=False)

    elif nen == 6:
        for ele_conn in eles_conn:
            x = nds_crd[[ele_conn[0], ele_conn[3], ele_conn[1], ele_conn[4],
                         ele_conn[2], ele_conn[5]], 0]
            y = nds_crd[[ele_conn[0], ele_conn[3], ele_conn[1], ele_conn[4],
                         ele_conn[2], ele_conn[5]], 1]
            plt.fill(x, y, edgecolor=ec, lw=lw, fill=False)

    elif nen == 9:
        for ele_conn in eles_conn:
            x = nds_crd[[ele_conn[0], ele_conn[4], ele_conn[1], ele_conn[5],
                         ele_conn[2], ele_conn[6], ele_conn[3], ele_conn[7]],
                        0]
            y = nds_crd[[ele_conn[0], ele_conn[4], ele_conn[1], ele_conn[5],
                         ele_conn[2], ele_conn[6], ele_conn[3], ele_conn[7]],
                        1]
            plt.fill(x, y, edgecolor=ec, lw=lw, fill=False)

    elif nen == 8:
        for ele_conn in eles_conn:
            x = nds_crd[[ele_conn[0], ele_conn[4], ele_conn[1], ele_conn[5],
                         ele_conn[2], ele_conn[6], ele_conn[3], ele_conn[7]],
                        0]
            y = nds_crd[[ele_conn[0], ele_conn[4], ele_conn[1], ele_conn[5],
                         ele_conn[2], ele_conn[6], ele_conn[3], ele_conn[7]],
                        1]
            plt.fill(x, y, edgecolor=ec, lw=lw, fill=False)


# def plot_stress_2d(nds_val, mesh_outline=1, cmap='turbo', levels=50,
#                    fig_wi_he=False, fig_lbrt=False, ax=False):
def plot_stress_2d(nds_val, mesh_outline=1, cmap='turbo', levels=50):
    """
    Plot stress distribution of a 2d elements of a 2d model.

    Args:
        nds_val (ndarray): the values of a stress component, which can
            be extracted from sig_out array (see sig_out_per_node
            function)

        mesh_outline (int): 1 - mesh is plotted, 0 - no mesh plotted.

        cmap (str): Matplotlib color map (default is 'turbo')

    Usage:
        See demo_quads_4x4.py example.
    """

    node_tags, ele_tags_all = ops.getNodeTags(), ops.getEleTags()

    ele_tags = stress_2d_ele_tags_only(ele_tags_all)

    n_nodes, n_eles = len(node_tags), len(ele_tags)

    ele_classtag = ops.getEleClassTags(ele_tags[0])[0]

    if (ele_classtag == EleClassTag.tri3n):
        nen = 3

    elif (ele_classtag == EleClassTag.tri6n):
        nen = 6

    elif (ele_classtag == EleClassTag.quad4n):
        nen = 4

    elif (ele_classtag == EleClassTag.quad8n):
        nen = 8

    elif (ele_classtag == EleClassTag.quad9n):
        nen = 9

    # nen = np.shape(ops.eleNodes(ele_tags[0]))[0]

    # idiom coordinates as ordered in node_tags
    # use node_tags.index(tag) for correspondence
    nds_crd = np.zeros((n_nodes, 2))
    for i, node_tag in enumerate(node_tags):
        nds_crd[i] = ops.nodeCoord(node_tag)

    # from utils / sig_out_per_node
    # fixme: if this can be simplified
    # index (starts from 0) to node_tag correspondence
    # (a) data in np.array of integers
    # nodes_tag_count = np.zeros((n_nodes, 2), dtype=int)
    # nodes_tag_count[:, 0] = node_tags
    #
    # correspondence indx and node_tag is in node_tags.index
    # after testing remove the above

    eles_conn = np.zeros((n_eles, nen), dtype=int)

    for i, ele_tag in enumerate(ele_tags):
        ele_node_tags = ops.eleNodes(ele_tag)

        for j, ele_node_tag in enumerate(ele_node_tags):
            eles_conn[i, j] = node_tags.index(ele_node_tag)

    if (ele_classtag == EleClassTag.tri3n):
        tris_conn = eles_conn
        nds_crd_all = nds_crd
        nds_val_all = nds_val

    elif (ele_classtag == EleClassTag.quad4n):
        tris_conn, nds_c_crd, nds_c_val = \
            quads_to_4tris(eles_conn, nds_crd, nds_val)

        nds_crd_all = np.vstack((nds_crd, nds_c_crd))
        nds_val_all = np.hstack((nds_val, nds_c_val))

    elif (ele_classtag == EleClassTag.tri6n):
        tris_conn = tris6n_to_4tris(eles_conn)
        nds_crd_all = nds_crd
        nds_val_all = nds_val

    elif (ele_classtag == EleClassTag.quad8n):
        tris_conn, nds_c_crd, nds_c_val = quads_to_8tris_8n(eles_conn,
                                                            nds_crd, nds_val)

        nds_crd_all = np.vstack((nds_crd, nds_c_crd))
        nds_val_all = np.hstack((nds_val, nds_c_val))

    elif (ele_classtag == EleClassTag.quad9n):
        tris_conn = quads_to_8tris_9n(eles_conn)
        nds_crd_all = nds_crd
        nds_val_all = nds_val



    # 1. plot contour maps
    triangulation = tri.Triangulation(nds_crd_all[:, 0],
                                      nds_crd_all[:, 1],
                                      tris_conn)

    plt.tricontourf(triangulation, nds_val_all, levels=levels, cmap=cmap)

    # 2. plot original mesh (quad) without subdivision into triangles
    if mesh_outline:
        # plot_mesh_2d(nds_crd, eles_conn, ax=ax)
        plot_mesh_2d(nds_crd, eles_conn)

    plt.colorbar()
    plt.axis('equal')
    plt.grid(False)


# see also quads_to_8tris_9n
def quads_to_8tris_8n(quads_conn, nds_crd, nds_val):
    """
    Get triangles connectivity, coordinates and new values at quad centroids.

    Args:
        quads_conn (ndarray):

        nds_crd (ndarray):

        nds_val (ndarray):

    Returns:
        tris_conn, nds_c_crd, nds_c_val (tuple):

    Notes:
        Triangles connectivity array is based on
        quadrilaterals connectivity.
        Each quad is split into eight triangles.
        New nodes are created at the quad centroid.

    See also:
        function: quads_to_8tris_9n, quads_to_4tris
    """
    n_quads, _ = quads_conn.shape
    n_nds, _ = nds_crd.shape

    # coordinates and values at quad centroids _c_
    nds_c_crd = np.zeros((n_quads, 2))
    nds_c_val = np.zeros(n_quads)

    tris_conn = np.zeros((8*n_quads, 3), dtype=int)

    for i, quad_conn in enumerate(quads_conn):
        j = 8*i
        n0, n1, n2, n3, n4, n5, n6, n7 = quad_conn

        # quad centroids
        # nds_c_crd[i] = np.array([np.sum(nds_crd[[n0, n1, n2, n3], 0])/4.,
        #                          np.sum(nds_crd[[n0, n1, n2, n3], 1])/4.])
        # nds_c_val[i] = np.sum(nds_val[[n0, n1, n2, n3]])/4.
        nds_c_crd[i] = quad8n_val_at_center(nds_crd[[n0, n1, n2, n3,
                                                     n4, n5, n6, n7]])
        nds_c_val[i] = quad8n_val_at_center(nds_val[[n0, n1, n2, n3,
                                                     n4, n5, n6, n7]])

        # triangles connectivity
        tris_conn[j] = np.array([n0, n4, n_nds+i])
        tris_conn[j+1] = np.array([n4, n1, n_nds+i])
        tris_conn[j+2] = np.array([n1, n5, n_nds+i])
        tris_conn[j+3] = np.array([n5, n2, n_nds+i])
        tris_conn[j+4] = np.array([n2, n6, n_nds+i])
        tris_conn[j+5] = np.array([n6, n3, n_nds+i])
        tris_conn[j+6] = np.array([n3, n7, n_nds+i])
        tris_conn[j+7] = np.array([n7, n0, n_nds+i])

    return tris_conn, nds_c_crd, nds_c_val


# see also quads_to_8tris_8n
def quads_to_8tris_9n(quads_conn):
    """
    Get triangles connectivity, coordinates and new values at quad centroids.

    Args:
        quads_conn (ndarray):

        nds_crd (ndarray):

        nds_val (ndarray):

    Returns:
        tris_conn, nds_c_crd, nds_c_val (tuple):

    Notes:
        Triangles connectivity array is based on
        quadrilaterals connectivity.
        Each quad is split into eight triangles.
        New nodes are created at the quad centroid.

    See also:
        function: quads_to_8tris_8n, quads_to_4tris
    """
    n_quads, _ = quads_conn.shape
    # n_nds, _ = nds_crd.shape

    tris_conn = np.zeros((8*n_quads, 3), dtype=int)

    for i, quad_conn in enumerate(quads_conn):
        j = 8*i
        n0, n1, n2, n3, n4, n5, n6, n7, n8 = quad_conn

        # center nodeds already present in quad 9n

        # triangles connectivity
        tris_conn[j] = np.array([n0, n4, n8])
        tris_conn[j+1] = np.array([n4, n1, n8])
        tris_conn[j+2] = np.array([n1, n5, n8])
        tris_conn[j+3] = np.array([n5, n2, n8])
        tris_conn[j+4] = np.array([n2, n6, n8])
        tris_conn[j+5] = np.array([n6, n3, n8])
        tris_conn[j+6] = np.array([n3, n7, n8])
        tris_conn[j+7] = np.array([n7, n0, n8])

    # return tris_conn, nds_c_crd, nds_c_val
    return tris_conn


def quad8n_val_at_center(vals):
    """
    Calculate values at the center of 8-node quad element.

    """

    val_c1 = -np.mean(vals[[0, 1, 2, 3]], axis=0) + 2*np.mean(vals[[4, 5, 6, 7]], axis=0)  # noqa: E501
    # val_c1 = -np.sum(vals[[0, 1, 2, 3]], axis=0)/4 + np.sum(vals[[4, 5, 6, 7]], axis=0)/2  # noqa: E501
    # val_c2 = -np.sum(vals[[0, 1, 2, 3]], axis=0)/4 + np.sum(vals[[4, 5, 6, 7]], axis=0)/2  # noqa: E501
    # val_c3 = np.mean(vals[[4, 5, 6, 7]], axis=0)
    # val_c4 = np.mean(vals[[0, 1, 2, 3]], axis=0)
    # val_c5 = np.mean(vals, axis=0)
    # return val_c1, val_c2, val_c3, val_c4, val_c5

    return val_c1


def plot_stress(stress_str, mesh_outline=1, cmap='turbo', levels=50):
    """Plot stress distribution of the model.

    Args:
        stress_str (string): stress component string. Available options are:
            'sxx', 'syy', 'sxy', 'vmis', 's1', 's2', 'alpha'

        mesh_outline (int): 1 - mesh is plotted, 0 - no mesh plotted.

        cmap (str): Matplotlib color map (default is 'turbo')

        levels (int): number and positions of the contour lines / regions.

    Usage:
        ::

            opsv.plot_stress('vmis')
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
            plt.show()

    See also:

    :ref:`opsvis_sig_out_per_node`
    """

    ndim = ops.getNDM()[0]

    if ndim == 2:
        _plot_stress_2d(stress_str, mesh_outline, cmap, levels)
        # if axis_off:
        #     plt.axis('off')

    # not implemented yet
    # elif ndim == 3:
    #     _plot_stress_3d(stress_str, mesh_outline, cmap, levels)

    else:
        print(f'\nWarning! ndim: {ndim} not implemented yet.')

    # plt.show()  # call this from main py file for more control


def _plot_stress_2d(stress_str, mesh_outline, cmap, levels):
    """See documentation for plot_stress command"""

    # node_tags = ops.getNodeTags()
    # ele_tags = ops.getEleTags()
    # n_nodes = len(node_tags)

    # second version - better - possible different types
    # of elements (mix of quad and tri)
    # for ele_tag in ele_tags:
    #     nen = np.shape(ops.eleNodes(ele_tag))[0]

    # avoid calculating and storing all stress components
    # sig_out = sig_out_per_node(stress_str)
    # switcher = {'sxx': 0,
    #             'syy': 1,
    #             'sxy': 2,
    #             'svm': 3,
    #             'vmis': 3,
    #             's1': 4,
    #             's2': 5,
    #             'angle': 6}

    # nds_val = sig_out[:, switcher[stress_str]]

    nds_val = sig_component_per_node(stress_str)
    plot_stress_2d(nds_val, mesh_outline, cmap, levels)
