.. opsvis documentation master file, created by
   sphinx-quickstart on Sat Aug 22 13:35:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Opsvis documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   plot_model
   plot_defo
   plot_load
   plot_mode_shape
   section_force_diagram_2d
   section_force_diagram_3d
   plot_stress
   plot_strain
   plot_stress_2d
   plot_extruded_model_rect_section_3d
   anim_defo
   anim_mode
   plot_fiber_section
   fib_sec_list_to_cmds
   sig_out_per_node
   examples

Warning - Incompatible changes!

Starting from Opsvis ver. 1.0.1 (May 2022):

#. the ``plot_supports_and_loads_2d`` function is removed and instead
   ``plot_loads_2d`` is available. The model supports can now be shown in
   the ``plot_model`` function. The additiona function argument
   ``node_supports`` (default is True) can be switched off (e.g.
   ``plot_model(node_supports=False)``).

#. the ``plot_mesh_with_ips_2d()`` function is removed.

#. See and use the updated example ``.py`` files (:doc:`examples`).

#. Now the opsvis plotting functions use the Python dictionary
   (instead of the string) for the maptplotlib line formatting.
   Example: Use ``fmt_model = {'color': 'blue', 'linestyle': 'solid',
   'linewidth': 1.2, 'marker': '.', 'markersize': 6}`` instead of the
   previous ``fmt_defo = 'b-'`` string format. Also the ``lw``
   (linewidth) arguments are removed from the plotting functions as
   they can be defined in the ``fmt_*`` dictionary. This feature has
   been implemented as suggested by `mbbatukan
   <https://github.com/sewkokot/opsvis/issues/11>`_.


Opsvis is an OpenSeesPy postprocessing and visualization module
written by Seweryn Kokot (Opole University of Technology, Poland).

For OpenSeesPy documentation click the following link:
`OpenSeesPy documentation <https://openseespydoc.readthedocs.io>`_

Opsvis can be mainly useful for students when learning the
fundamentals of structural analysis (interpolated deformation of frame
structures (static images or animations), section force distribution
of frame structures, stress distribution in triangle, quadrilateral 2d
elements, orientation of frame members in 3d space, fibers of a cross
section, static and animated eigenvalue mode shapes etc.). This way,
we can lower the bar in teaching and learning OpenSees at earlier
years of civil engineering studies. However the visualization features
for OpenSees can also be helpful for research studies.

Opsvis offers the following plots:

- interpolated deformation of frame structures,
- stresses of triangular and (four, eight and nine-node) quadrilateral
  2d elements (calculation of Huber-Mieses-Hencky equivalent stress,
  principal stresses),
- fibers of cross-sections,
- models with extruded cross sections
- animation of deformation (from time history analysis) and mode shapes.

Installation
============

``pip install opsvis``

Note the name of the PyPi package is without the underscore ``_``.

Usage
=====

To use Opsvis in OpenSeesPy
scripts, your ``.py`` file should start as follows: ::

	import openseespy.opensees as ops
	import opsvis as opsv
	import matplotlib.pyplot as plt

	# ... your OpenSeesPy model and analysis commands ...
	opsv.plot_model()
	sfac = opsv.plot_defo()


Commands
========

The main commands related to various aspects of OpenSees model
visualization are as follows:

#. :doc:`plot_model`
#. :doc:`plot_defo`
#. :doc:`plot_loads_2d`
#. :doc:`plot_mode_shape`
#. :doc:`section_force_diagram_2d`
#. :doc:`section_force_diagram_3d`
#. :doc:`plot_stress_2d`
#. :doc:`plot_extruded_model_rect_section_3d`
#. :doc:`anim_defo`
#. :doc:`anim_mode`
#. :doc:`plot_fiber_section`

Helper functions include:

#. :doc:`fib_sec_list_to_cmds`
#. :doc:`sig_out_per_node`

For examples go to: :doc:`examples`.

Notes:
======

* matplotlib's ``plt.axis('equal')`` does not work for 3d plots
  therefore right angles are not guaranteed to be 90 degrees on the
  plots

* ``plot_fiber_section`` is inspired by Matlab ``plotSection.zip``
  written by D. Vamvatsikos available at
  http://users.ntua.gr/divamva/software.html
