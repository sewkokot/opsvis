import openseespy.opensees as ops
import opsvis as opsv

import matplotlib.pyplot as plt

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 2)
ops.node(1, 0., 0.)
ops.node(2, 0., 1.)
ops.node(3, 0., 2.)
ops.node(4, 0., 3.)
ops.node(5, 0., 4.)
ops.node(6, 1., 0.)
ops.node(7, 1., 1.)
ops.node(8, 1., 2.)
ops.node(9, 1., 3.)
ops.node(10, 1., 4.)
ops.node(11, 2., 0.)
ops.node(12, 2., 1.)
ops.node(13, 2., 2.)
ops.node(14, 2., 3.)
ops.node(15, 2., 4.)
ops.node(16, 3., 0.)
ops.node(17, 3., 1.)
ops.node(18, 3., 2.)
ops.node(19, 3., 3.)
ops.node(20, 3., 4.)
ops.node(21, 4., 0.)
ops.node(22, 4., 1.)
ops.node(23, 4., 2.)
ops.node(24, 4., 3.)
ops.node(25, 4., 4.)

ops.nDMaterial('ElasticIsotropic', 1, 1000, 0.3)

ops.element('quad', 1, 1, 6, 7, 2, 1, 'PlaneStress', 1)
ops.element('quad', 2, 2, 7, 8, 3, 1, 'PlaneStress', 1)
ops.element('quad', 3, 3, 8, 9, 4, 1, 'PlaneStress', 1)
ops.element('quad', 4, 4, 9, 10, 5, 1, 'PlaneStress', 1)
ops.element('quad', 5, 6, 11, 12, 7, 1, 'PlaneStress', 1)
ops.element('quad', 6, 7, 12, 13, 8, 1, 'PlaneStress', 1)
ops.element('quad', 7, 8, 13, 14, 9, 1, 'PlaneStress', 1)
ops.element('quad', 8, 9, 14, 15, 10, 1, 'PlaneStress', 1)
ops.element('quad', 9, 11, 16, 17, 12, 1, 'PlaneStress', 1)
ops.element('quad', 10, 12, 17, 18, 13, 1, 'PlaneStress', 1)
ops.element('quad', 11, 13, 18, 19, 14, 1, 'PlaneStress', 1)
ops.element('quad', 12, 14, 19, 20, 15, 1, 'PlaneStress', 1)
ops.element('quad', 13, 16, 21, 22, 17, 1, 'PlaneStress', 1)
ops.element('quad', 14, 17, 22, 23, 18, 1, 'PlaneStress', 1)
ops.element('quad', 15, 18, 23, 24, 19, 1, 'PlaneStress', 1)
ops.element('quad', 16, 19, 24, 25, 20, 1, 'PlaneStress', 1)

ops.fix(1, 1, 1)
ops.fix(6, 1, 1)
ops.fix(11, 1, 1)
ops.fix(16, 1, 1)
ops.fix(21, 1, 1)

ops.equalDOF(2, 22, 1, 2)
ops.equalDOF(3, 23, 1, 2)
ops.equalDOF(4, 24, 1, 2)
ops.equalDOF(5, 25, 1, 2)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(15, 0., -1.)

ops.analysis('Static')
ops.analyze(1)

# - plot model
opsv.plot_model()
plt.axis('equal')

# plt.figure()
opsv.plot_load()

# - plot deformation
opsv.plot_defo(unDefoFlag=1)
plt.axis('equal')

jstr = 'sxx'
plt.figure()
opsv.plot_stress(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

jstr = 'syy'
plt.figure()
opsv.plot_stress(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

jstr = 'sxy'
plt.figure()
opsv.plot_stress(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

jstr = 'vmis'
plt.figure()
opsv.plot_stress(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')


jstr = 'exx'
plt.figure()
opsv.plot_strain(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

jstr = 'eyy'
plt.figure()
opsv.plot_strain(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

jstr = 'exy'
plt.figure()
opsv.plot_strain(jstr)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(f'{jstr}')

plt.show()

exit()
