# default settings

# fmt: format string setting color, marker and linestyle
# check documentation on matplotlib's plot

# element end nodes
fmt_nodes = {'color': 'red', 'linestyle': 'None', 'linewidth': 1.2, 'marker': 's', 'markersize': 6}

# initial model
fmt_model = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '.', 'markersize': 6}
fmt_model_truss = {'color': 'green', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': 'o', 'markersize': 6, 'markerfacecolor': 'white'}
fmt_model_nodes_only = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '.', 'markersize': 6}
fmt_model_loads = {'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}
fmt_model_rigid_offset = {'color': 'black', 'linestyle': 'solid', 'linewidth': 3.2, 'marker': '.', 'markersize': 1}

fmt_model_joint2d = {'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}
fmt_model_secforce = {'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}

# deformed model
fmt_defo = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}
fmt_defo_faces = {'linewidths': 1, 'edgecolors': 'k', 'alpha': 0.5}
fmt_defo_zeroLenght = {'color': 'blue', 'linestyle': 'dashed', 'linewidth': 1.2, 'marker': '*', 'markersize': 6, 'markerfacecolor': 'white'}

# undeformed model
fmt_undefo = {'color': 'green', 'linestyle': (0, (1, 5)), 'linewidth': 1.2, 'marker': '', 'markersize': 1}
fmt_undefo_faces = {'linewidths': 1, 'linestyles': 'dotted', 'edgecolors': 'g', 'facecolors': 'w', 'alpha': 0.5}

# section forces
fmt_secforce1 = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 2.0, 'marker': '',
                 'markersize': 1, 'solid_capstyle': 'round', 'solid_joinstyle': 'round',
                 'dash_capstyle': 'butt', 'dash_joinstyle': 'round'}
fmt_secforce2 = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.0, 'marker': '',
                 'markersize': 1, 'solid_capstyle': 'round', 'solid_joinstyle': 'round',
                 'dash_capstyle': 'butt', 'dash_joinstyle': 'round'}

fmt_secforce_tension = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 2.0, 'marker': '', 'markersize': 1}
fmt_secforce_compression = {'color': 'red', 'linestyle': 'solid', 'linewidth': 2.0, 'marker': '', 'markersize': 1}

fmt_gauss_points = {'color': 'firebrick', 'linestyle': 'None', 'linewidth': 2.0, 'marker': 'X', 'markersize': 5}

# figure left right bottom top offsets
fig_lbrt = (.04, .04, .96, .96)
fig_lbrt_model = (.08, .08, .985, .94)
fig_lbrt_secforces = (.04, .08, .985, .92)
fig_lbrt_defo = (.08, .08, .985, .92)
fig_lbrt_mode = (.08, .08, .985, .92)

# azimuth and elevation in degrees
az_el = (-60., 30.)

# figure width and height in centimeters
fig_wi_he = (16., 10.)


class EleClassTag:
    """ELE_TAG constants defined in SRC/classTags.h"""
    truss = 12
    trussSection = 13
    CorotTruss = 14
    ZeroLength = 19
    ZeroLengthSection = 20
    CoupledZeroLength = 26
    ElasticBeam2d = 3
    ElasticBeam3d = 5
    DispBeamColumn2d = 62
    DispBeamColumn3d = 64
    ForceBeamColumn2d = 73
    ForceBeamColumn3d = 74
    TimoshenkoBeamColumn2d = 63
    TimoshenkoBeamColumn3d = 631
    ElasticTimoshenkoBeam2d = 145
    ElasticTimoshenkoBeam3d = 146
    tri3n = 33
    tri6n = 209
    quad4n = 31
    quad4n3d = 32
    quad9n = 207
    quad8n = 208
    SSPquad = 119
    EnhancedQuad = 59
    brick20n = 49
    brick8n = 56
    SSPbrick = 121
    FourNodeTetrahedron = 179
    TenNodeTetrahedron = 256
    TenNodeTetrahedronSK = 1790
    ASDShellQ4 = 203
    ASDShellT3 = 204
    ShellMITC4 = 53
    ShellMITC9 = 54
    ShellDKGQ = 156
    ShellNLDKGQ = 157
    ShellDKGT = 167
    ShellNLDKGT = 168        
    Joint2D = 82
    Joint3D = 83
    MVLEM = 162
    SFI_MVLEM = 163
    MVLEM_3D = 212
    SFI_MVLEM_3D = 213
    TwoNodeLink = 86


class LoadTag:
    """LOAD_TAG constants defined in SRC/classTags.h"""
    Beam2dUniformLoad = 3
    Beam2dUniformLoad_ndata = 2

    Beam2dPointLoad = 4
    Beam2dPointLoad_ndata = 3

    Beam3dUniformLoad = 5
    Beam3dUniformLoad_ndata = 3

    Beam3dPointLoad = 6
    Beam3dPointLoad_ndata = 4

    BrickSelfWeight = 7
    Beam2dTempLoad = 8
    SurfaceLoader = 9
    SelfWeight = 10
    Beam2dThermalAction = 11

    Beam2dPartialUniformLoad = 12
    Beam2dPartialUniformLoad_ndata = 6

    Beam3dPartialUniformLoad = 121
    Beam3dPartialUniformLoad_ndata = 5

    Beam3dThermalAction = 13
    ShellThermalAction = 14
    NodalThermalAction = 15
    ThermalActionWrapper = 16
    LysmerVelocityLoader = 17
