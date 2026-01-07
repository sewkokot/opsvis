================
 Getting Started
================

Getting started in OpenSeesPy and opsvis is straightforward especially
if you are already using OpenSeesPy

Under Windows do the following steps:

1. download and install Anaconda distribution from `<https://www.anaconda.com/download>`_

  Please check the current Python (e.g. 3.12) version required for
  OpenSeesPy (e. g. 3.7.0.4). More details can be found in `OpenSeesPy documentation <https://openseespydoc.readthedocs.io>`_

2. After installing Anaconda, run its Powershell Prompt and type ::

    pip install openseespy
    pip install opsvis

3. Launch Spyder from Anaconda

4. Download Python example files from this opsvis website

5. Open an example .py file in Spyder editor

6. Note that opsvis commands begins with 'opsv.' as it is specified in
   the 'import opsvis as opsv' at the beginning of each example file.

7. Run the example: (1) throught Menu - Run or (2) green triangle
   button in the toolbar or (3) pressing F1 shortcut key.

8. After running the example the plots should appear in the Plots tab
   of Spyder.

   Alternatively (recommended), the plots can be displayed outside the
   Spyder window. To achieve this, change the following settings.
   Go to Menu - Preferences - Ipython Console - Graphics tab and here
   change the 'Inline' option to 'Automatic'. Click Apply and Ok.

9. Check other opsvis examples and read the opsvis commands (e.g.
   plot_model, plot_defo, plot_load, plot_stress
   documentation.

10. Basic usage

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


#. :doc:`examples`
