.. opsvis documentation master file, created by
   sphinx-quickstart on Sat Aug 22 13:35:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Opsvis documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   GettingStarted
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

Check out the :doc:`GettingStarted` page on how to get started in
OpenSeesPy and Opsvis.


Notes:
======

* ``plot_fiber_section`` is inspired by Matlab ``plotSection.zip``
  written by D. Vamvatsikos available at
  http://users.ntua.gr/divamva/software.html
