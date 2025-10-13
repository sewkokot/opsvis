import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt

ops.model('basic', '-ndm', 2, '-ndf', 2)

Ls, Hs = 10, 4
nx, ny = Ls, Hs
E0, Ec, Es = 1000, 20.e6, 200.e6

b, h = .2, .3
A, Iz, As, c = b*h, b*h**3./12., 0.0002, 0.025
y1, z1 = h/2., b/2.
nFibZ, nFib = 1, 20

sec_tag, nip, beam_integ_tag, transf_tag = 1, 5, 1, 1
mat0_tag, mat1_tag, mat2_tag = 1, 2, 3

ops.nDMaterial('ElasticIsotropic', mat0_tag, E0, .3)
ops.uniaxialMaterial('Elastic', mat1_tag, Ec)
ops.uniaxialMaterial('Elastic', mat2_tag, Es)

ops.block2D(nx, ny, 1, 1, 'quad', mat0_tag, 'PlaneStress', 1,
            1, 0., 0.,
            2, Ls, 0.,
            3, Ls, Hs,
            4, 0., Hs)

ops.fixY(0., 1, 1); ops.fixX(0., 1, 0); ops.fixX(Ls, 1, 0)

ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(100, 2., 4.); ops.node(101, 2., 8.)
ops.node(102, 1., 4.); ops.node(103, 3., 4.)
ops.node(110, 8., 4.); ops.node(111, 8., 8.)
ops.node(112, 7., 4.); ops.node(113, 9., 4.)

ops.equalDOF(46, 102, 1, 2); ops.equalDOF(47, 100, 1, 2)
ops.equalDOF(48, 103, 1, 2); ops.equalDOF(52, 112, 1, 2)
ops.equalDOF(53, 110, 1, 2); ops.equalDOF(54, 113, 1, 2)

fib_sec_1 = [['section', 'Fiber', sec_tag, '-GJ', 1.0e6],
             ['patch', 'rect', mat1_tag, nFib, nFibZ, -y1, -z1, y1, z1],
             ['layer', 'straight', mat2_tag, 3, As, y1-c, z1-c, y1-c, c-z1],
             ['layer', 'straight', mat2_tag, 3, As, c-y1, z1-c, c-y1, c-z1]
             ]

matcolor = ['r', 'gold', 'w']
opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
plt.axis('equal')

opsv.fib_sec_list_to_cmds(fib_sec_1)

ops.beamIntegration('Lobatto', beam_integ_tag, sec_tag, nip)

ops.geomTransf('Linear', transf_tag)
ops.element('elasticBeamColumn', 101, 100, 101, A, Ec, Iz, 1)
ops.element('elasticBeamColumn', 111, 110, 111, A, Ec, Iz, 1)

ops.element('forceBeamColumn', 121, 101, 111, transf_tag, beam_integ_tag)

ops.element('elasticBeamColumn', 103, 102, 100, A, Ec, Iz, transf_tag)
ops.element('elasticBeamColumn', 104, 100, 103, A, Ec, Iz, transf_tag)

ops.element('elasticBeamColumn', 113, 112, 110, A, Ec, Iz, transf_tag)
ops.element('elasticBeamColumn', 114, 110, 113, A, Ec, Iz, transf_tag)

Px, Py, Wy = 15., 2., -3.

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(101, Px, -Py, 0.)

ops.eleLoad('-ele', 121, '-type', '-beamUniform', Wy)

ops.algorithm('Linear')
ops.analysis('Static')
ops.analyze(1)

opsv.plot_model()
opsv.plot_load()

opsv.plot_defo(node_supports=0)
plt.title('opsv.plot_defo(node_supports=0)')

sfac = 20.

_, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('sxx', sfac=sfac, ax=ax)
plt.title('Stress sxx')

_, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('syy', sfac=sfac, ax=ax)
plt.title('Stress syy')

_, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('sxy', sfac=sfac, ax=ax)
plt.title('Stress sxy')

_, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('vmis', sfac=sfac, ax=ax)
plt.title('Stress vmis')

opsv.section_force_diagram_2d('M', 1.e-1, number_format='.1f')
plt.title('Bending moment distribution')

plt.show()
