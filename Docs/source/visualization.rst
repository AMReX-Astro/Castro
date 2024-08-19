*************
Visualization
*************

There are a large number of tools that can be used to read in Castro
or AMReX data and make plots.  These tools all work from Castro
plotfiles.  Here we give an overview of the variables in plotfiles and
controlling their output, as well as some of the tools that can be used
for visualization.





Visualization Tools
===================

amrvis
------

Our favorite visualization tool is amrvis. We heartily encourage you
to build the amrvis2d and amrvis3d executables, and to try using them
to visualize your data. A very useful feature is View/Dataset, which
allows you to actually view the numbers – this can be handy for
debugging. You can modify how many levels of data you want to see,
whether you want to see the grid boxes or not, what palette you use,
etc.

If you like to have amrvis display a certain variable, at a certain
scale, when you first bring up each plotfile (you can always change it
once the amrvis window is open), you can modify the amrvis.defaults
file in your directory to have amrvis default to these settings every
time you run it. The directories CoreCollapse, HSE_test, Sod and
Sedov have amrvis.defaults files in them. If you are working in a new
run directory, simply copy one of these and modify it.

VisIt
-----

VisIt is also a great visualization tool, and it directly handles our
plotfile format (which it calls Boxlib). For more information check
out ``visit.llnl.gov``.

[Useful tip:] To use the Boxlib3D plugin, select it from File
:math:`\rightarrow` Open file :math:`\rightarrow` Open file as type Boxlib, and
then the key is to read the Header file, ``plt00000/Header``, for example,
rather than telling it to read ``plt00000``.

yt
--

yt is a free and open-source software that provides data analysis and
publication-level visualization tools for astrophysical simulation
results such as those Castro produces. As yt is script-based, it’s not
as easy to use as VisIt, and certainly not as easy as amrvis, but the
images can be worth it! Here we do not flesh out yt, but give an
overview intended to get a person started. Full documentation and
explanations from which this section was adapted can be found at
http://yt-project.org/doc/index.html.

Example notebook
^^^^^^^^^^^^^^^^

Using the plotfiles generated in the example in the :doc:`getting_started` section, here we demonstrate how to use ``yt`` to load and visualize data. This section was generated from a Jupyter notebook which can be found in ``Docs/source/yt_example.ipynb`` in the Castro repo.

.. include:: yt_example.rst
