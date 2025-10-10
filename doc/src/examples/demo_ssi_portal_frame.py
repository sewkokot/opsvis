import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt

ops.model('basic', '-ndm', 2, '-ndf', 2)

Ls, Hs = 10, 4
nx, ny = Ls, Hs

ops.nDMaterial('ElasticIsotropic', 1, 1000., .3)

ops.block2D(nx, ny, 1, 1, 'quad', 1, 'PlaneStress', 1,
            1, 0., 0.,
            2, Ls, 0.,
            3, Ls, Hs,
            4, 0., Hs)

ops.fixY(0., 1, 1)
ops.fixX(0., 1, 0)
ops.fixX(Ls, 1, 0)

ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(100, 2., 4.)
ops.node(101, 2., 8.)
ops.node(102, 1., 4.)
ops.node(103, 3., 4.)
ops.node(110, 8., 4.)
ops.node(111, 8., 8.)
ops.node(112, 7., 4.)
ops.node(113, 9., 4.)

ops.equalDOF(46, 102, 1, 2)
ops.equalDOF(47, 100, 1, 2)
ops.equalDOF(48, 103, 1, 2)
ops.equalDOF(52, 112, 1, 2)
ops.equalDOF(53, 110, 1, 2)
ops.equalDOF(54, 113, 1, 2)

b, h = .3, .2
A, Iz = b*h, b*h**3./12.
E = 25.e6

ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 101, 100, 101, A, E, Iz, 1)
ops.element('elasticBeamColumn', 103, 102, 100, A, E, Iz, 1)
ops.element('elasticBeamColumn', 104, 100, 103, A, E, Iz, 1)

ops.element('elasticBeamColumn', 111, 110, 111, A, E, Iz, 1)
ops.element('elasticBeamColumn', 113, 112, 110, A, E, Iz, 1)
ops.element('elasticBeamColumn', 114, 110, 113, A, E, Iz, 1)

ops.element('elasticBeamColumn', 121, 101, 111, A, E, Iz, 1)

Px, Py, Wy = 10., 1., -3.

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

sfac, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('sxx', sfac=sfac, ax=ax)
plt.title('Stress sxx')

sfac, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('syy', sfac=sfac, ax=ax)
plt.title('Stress syy')

sfac, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('sxy', sfac=sfac, ax=ax)
plt.title('Stress sxy')

sfac, ax = opsv.plot_defo(unDefoFlag=0, node_supports=0, sfac=sfac)
opsv.plot_stress('vmis', sfac=sfac, ax=ax)
plt.title('Stress vmis')

plt.show()

exit()
