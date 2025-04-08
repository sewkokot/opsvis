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
   plot_reactions
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

Opsvis is a postprocessing and visualization module for OpenSeesPy,
developed by Seweryn Kokot from Opole University of Technology,
Poland.

For OpenSeesPy documentation, visit:
`OpenSeesPy documentation <https://openseespydoc.readthedocs.io>`_

Key Features of Opsvis:

- Orientation of frame members in 3D space,
- Interpolated deformations of frame structures (as static images or
  animations),
- Section force distributions,
- Stress visualization for triangular and quadrilateral 2D elements
  (4-, 8-, and 9-node), including Huber-Mieses-Hencky equivalent stress,
  principal stresses.
- Visualization of cross-section fibers,
- Animations of time-history deformations and mode shapes.

These capabilities help lower the barrier to entry for teaching and
learning OpenSees in civil engineering education. Additionally, Opsvis
provides visualization features that can support advanced research
applications.


Installation
============

``pip install opsvis``


Usage
=====

To use Opsvis in your OpenSeesPy scripts, begin your ``.py`` file with the
following import: ::

	import openseespy.opensees as ops
	import opsvis as opsv
	import matplotlib.pyplot as plt

	# ... your OpenSeesPy model and analysis commands ...
	opsv.plot_model()
	opsv.plot_load()
	opsv.plot_reactions()
	sfac = opsv.plot_defo()


Commands
========

The main commands for visualizing various aspects of an OpenSees
model are as follows:

#. :doc:`plot_model`
#. :doc:`plot_defo`
#. :doc:`plot_load`
#. :doc:`plot_mode_shape`
#. :doc:`section_force_diagram_2d`
#. :doc:`section_force_diagram_3d`
#. :doc:`plot_stress`
#. :doc:`plot_strain`
#. :doc:`plot_stress_2d`
#. :doc:`plot_extruded_model_rect_section_3d`
#. :doc:`anim_defo`
#. :doc:`anim_mode`
#. :doc:`plot_fiber_section`

You can also make use of the following helper functions:

#. :doc:`fib_sec_list_to_cmds`
#. :doc:`sig_out_per_node`

Check out the :doc:`examples` page for usage demonstrations.

Notes:
======

* ``plot_fiber_section`` is inspired by Matlab ``plotSection.zip``
  written by D. Vamvatsikos available at
  http://users.ntua.gr/divamva/software.html
