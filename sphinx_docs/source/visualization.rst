*************
Visualization
*************

There are a large number of tools that can be used to read in Castro or BoxLib data and make plots. Here we give a brief overview of some
of the tools as well as some examples.

Controlling What’s in the PlotFile
==================================

There are a few options that can be set at runtime to control what
variables appear in the plotfile.

-  amr.plot_vars : this controls which of the main
   state variables appear in the plotfile. The default is for all of
   them to be stored. But you can specify a subset by name, e.g.

   ::

             amr.plot_vars = density
           

   to only store that subset.

-  amr.derive_plot_vars : this controls which of the
   derived variables to be stored in the plotfile. Derived variables
   are created only when the plotfile is being created, using the
   infrastructure provided by BoxLib to register variables and the
   associated Fortran routine to do the deriving
   (Derive_nd.F90).

   By default, no derived variables are stored. You can store all
   derived variables that Castro knows about by doing:

   ::

             amr.derive_plot_vars = ALL
           

   or a subset by explicitly listing them, e.g.,

   ::

             amr.derive_plot_vars = entropy pressure
           

amrvis
======

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
=====

VisIt is also a great visualization tool, and it directly handles our
plotfile format (which it calls Boxlib). For more information check
out visit.llnl.gov.

[Useful tip:] To use the Boxlib3D plugin, select it from File
:math:`\rightarrow` Open file :math:`\rightarrow` Open file as type Boxlib, and
then the key is to read the Header file, plt00000/Header, for example,
rather than telling to to read plt00000.

yt
==

yt is a free and open-source software that provides data analysis and
publication-level visualization tools for astrophysical simulation
results such as those CASTRO produces. As yt is script-based, it’s not
as easy to use as VisIt, and certainly not as easy as amrvis, but the
images can be worth it! Here we do not flesh out yt, but give an
overview intended to get a person started. Full documentation and
explanations from which this section was adapted can be found at
http://yt-project.org/doc/index.html.

yt can be installed by the following commands:

$ wget https://raw.githubusercontent.com/yt-project/yt/master/doc/install_script.sh

$ bash install_script.sh

This installs yt in your current directory. To update ytin the
future, simply do

$ conda update yt

assuming you have conda installed.

Castro-Specific Data
--------------------

yt was originally created for simple analysis and visualization of
data from the Enzo code. Since, it has grown to include support for a
variety of codes, including Castro. However, ytwill still sometimes
make assumptions, especially about data field names, that favor Enzo
and cause errors with Castro data. These problems can usually be
avoided by taking care to specify the data fields desired in
visualization. For example, Enzo’s density field is called
“Density,” and is the default for many plotting mechanisms when the
user does not specify the field. However, Castro does not have a field
called “Density”; instead, the density field is called “density.”
If a user does not specify a field while plotting with Castro data,
chances are that yt will try (and fail) to find “Density” and return
an error. As you will see in the examples, however, there is a way to
create your own fields from existing ones. You can use these derived
fields as you would use any other field.

There are also a few imperatives when it comes to reading in your
Castro simulation data and associated information. First and foremost
is that the inputs file for the simulation **must** exist in the
same directory as where the plotfile directory is located, and it
**must** be named “**inputs**.” yt reads information from the
inputs file such as the number of levels in the simulation run, the
number of cells, the domain dimensions, and the simulation time. yt will also optionally parse the probin file for pertinent information
if it is similarly included with the name “**probin**” in the same
directory as the plotfile of interest. When specifying a plotfile as
the data source for plots, you may simply call it by its directory
name, rather than using the Header file as in VisIt. As a final
caveat, the existence of the job_info file within the plotfile
directory is what currently distinguishes Castro data from MAESTRO
data in yt; unless you like surprises, we suggest you ensure your
plotfile has one.

Interacting with yt: Command Line and Scripting
-----------------------------------------------

ytis written completely in python (if you don’t have python, yt will
install it for you) and there are a number of different ways to
interact with it, including a web-based gui. Here we will cover
command-line yt and scripts/the python interactive prompt, but other
methods are outlined on the yt webpage at
http://yt-project.org/doc/interacting/index.html.

The first step in starting up yt is to activate the yt environment:

$ source $YT_DEST/bin/activate

From the command line you can create simple plots, perform simple
volume renderings, print the statistics of a field for your data set,
and do a few other things. Try $ yt to see a list of commands,
and $ yt :math:`<`\ command\ :math:`>` --help
to see the details of a command. The command line is the easiest way
to get quick, preliminary plots – but the simplicity comes at a
price, as yt will make certain assumptions for you. We could plot a
projection of density along the x-axis for the plotfile (yt calls it a
parameter file) plt_def_00020 by doing the following:

$ yt plot -p -a 0 -f density plt_def_00020

Or a temperature-based volume rendering with 14 contours:

$ yt render -f Temp --contours 14 plt_def_00020

Any plots created from the command line will be saved into a
subdirectory called “frames.” The command line is nice for fast
visualization without immersing yourself too much in the data, but
usually you’ll want to specify and control more details about your
plots. This can be done either through scripts or the python
interactive prompt. You can combine the two by running scripts within
the interactive prompt by the command

:math:`>>>` execfile(‘script.py’)

which will leave you in the interactive prompt, allowing you to
explore the data objects you’ve created in your script and debug
errors you may encounter. While in the yt environment, you can access
the interactive prompt by $ *python* or the shortcut

$ pyyt

Once you’re in the yt environment and in a .py script or the
interactive prompt, there are a couple of points to know about the
general layout of yt scripting. Usually there are five sections to a
yt script:

#. Import modules

#. Load parameter files and saved objects

#. Define variables

#. Create and modify data objects, image arrays, plots,
   etc. :math:`\rightarrow` this is the meat of the script

#. Save images and objects

Note that neither saving nor loading objects is necessary, but can be
useful when the creation of these objects is time-consuming, which is
often the case during identification of clumps or contours.

yt Basics
---------

The first thing you will always want to do is to import yt:

:math:`>>>` from yt.mods import \*

Under certain circumstances you will be required to import more, as we
will see in some of the examples, but this covers most of it,
including all of the primary functions and data objects provided by
yt. Next, you’ll need yt to access the plotfile you’re interested in
analyzing. Remember, you must have the “inputs” file in the same
directory:

:math:`>>>` ds = load(‘plt_def_00020’)

When this line is executed, it will print out some key parameters from
the simulation. However, in order to access information about all of
the fluid quantities in the simulation, we must use the “index”
object. (Note that for yt versions earlier than 3.0, this information
was contained in the “hierarchy” object; for these versions, replace
pf.index with pf. h in the following examples. The “hierarchy” object
was removed in yt-3.0 and its associated functionality for accessing data
was moved directly to the datasets themselves.) It contains the geometry
of the grid zones, their parentage relationships, and the fluid states
within each one. It is easily created:

:math:`>>>` ds.index

Upon execution, yt may print out a number of lines saying it’s adding
unknown fields to the list of fields. This is because Castro has
different names for fields than what yt expects. We can see what
fields exist through the commands

:math:`>>>` print ds.index.field_list

:math:`>>>` print ds.index.derived_field_list

There may not be any derived fields for Castro data. We can find out
the number of grids and cells at each level, the simulation time, and
information about the finest resolution cells:

:math:`>>>` ds.index.print_stats()

The dataset itself also stores a number of associated methods; for example,
you can find the value and location of the maximum of a field in the domain:

:math:`>>>` value, location = ds.find_max(‘density’)

(Note that in yt versions before 3.0, this type of method was primarily
associated with the hierarchy object and was accessed with ds.h.find_max.)

The list goes on. A full list of methods and attributes associated
with the index object (and most any yt object or function) can be
accessed by the help function:

:math:`>>>` help(pf.index)

You can also use :math:`>>>` *dir()* on an object or
function to find out which names it defines. Don’t be shy about
searching the yt documentation for help. Note that creating the
index object in its own line is not always needed before calling
functions like find_max; yt will construct it automatically if it
does not already exist.

Data Containers and Selection
-----------------------------

Sometimes, you’ll want to select, analyze, or plot only portions of
your simulation data. To that end, yt includes a way to create data
“containers” that select data based on geometric bounds or fluid
quantity values. There are many, including rays, cylinders, and clumps
(some in the examples, all described in the documentation), but the
easiest to create is a sphere, centered on the location of the maximum
density cell we found above:

:math:`>>>` my_data_container = ds.sphere(location, (5.0e4, ‘km’))

Here, specify that the radius is in units of kilometers using a dimensionful
quantity. When specifying distances in yt, the default is to use the
simulation-native unit named “code_length”, which for Castro is “cm”, and
if you just put in 5.0e4 instead of (5.0e4, ‘km’), you will get a 50,000 cm radius.
The pf.index.print_stats() command lists available units. We can access the data
within the container:

:math:`>>>` print my_data_container[‘density’]

:math:`>>>` print my_data_container.quantities[‘Extrema’]([‘density’, ‘pressure’])

When the creation of objects is time-consuming, it can be convenient
to save objects so they can be used in another session. To save an
object as part of the .yt file affiliated with the index:

:math:`>>>` pf.index.save_object(my_data_container, ‘sphere_to_analyze_later’)

Once it has been saved, it can be easily loaded later:

:math:`>>>` sphere_to_analyze = pf.index.load_object(‘sphere_to_analyze_later’)

Grid Inspection
---------------

yt also allows for detailed grid inspection. The index object
possesses an array of grids, from which we can select and examine
specific ones:

:math:`>>>` print pf.index.grids

:math:`>>>` my_grid = pf.index.grids[4]

Each grid is a data object that carries information about its
location, parentage relationships (grids within which it resides, and
grids that reside within it, at least in part), fluid quantities, and
more. Here are some of the commands:

:math:`>>>` print my_grid.Level

:math:`>>>` print my_grid_ActiveDimensions

:math:`>>>` print my_grid.LeftEdge

:math:`>>>` print my_grid.RightEdge

:math:`>>>` print my_grid.dds

(dds is the size of each cell within the grid).

:math:`>>>` print my_grid.Parent

:math:`>>>` print my_grid.Children[2].LeftEdge

:math:`>>>` print my_grid[‘Density’]

You can examine which cells within the grid have been refined with the
child_mask attribute, a representative array set to zero everywhere
there is finer resolution.To find the fraction of your grid that isn’t
further refined:

:math:`>>>`\ print my_grid.child_mask.sum()/float(my_grid.ActiveDimensions.prod())

Rather than go into detail about the many possibilities for plotting
in yt, we’ll provide some examples.

Example Scripts
---------------

In these examples, we investigate 3-D simulation data of two stars
orbiting in the center of the domain, which is a box of sides
:math:`10^{10}\:cm`.

*# Pressure Contours*

.. raw:: latex

   \setlength{\parskip}{0pt}

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’

pf.index

*# Most Castro fields have no inherent units, so we add them in,
in the form of a raw string*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# with some LaTeX-style formatting.*

pf.field_info[field]._units = r‘\\rm{Ba}’

*# SlicePlot parameters include: parameter file, axis, field, window width (effectively the*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# x and y zoom), and fontsize. We can also create projections with ProjectionPlot().*

p = SlicePlot(pf, ‘z’, field, width=((5.0e9, ‘cm’), (3.0e9, ‘cm’)),

fontsize=13)

*# Zlim is the range of the colorbar. In other words, the range of the data we want to display.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Names for many colormaps can be found at wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps.*

p.set_zlim(field, 2.85e13, 2.95e13)

p.set_cmap(field, ‘jet’)

*# Here we add 5 density contour lines within certain limits on top of the image. We overlay*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# our finest grids with a transparency of 0.2 (lower is more transparent). We add a quiver*

*# plot with arrows every 16 pixels with x_velocity in the x-direction and y_velocity in*

*# the y-direction. We also mark the center with an ‘x’ and label one of our stars.*

p.annotate_contour(‘density’, clim=(1.05e-4, 1.16e-4), ncont=5, label=False)

p.annotate_grids(alpha=0.2, min_level=2)

p.annotate_quiver(‘x_velocity’, ‘y_velocity’, factor=16)

p.annotate_marker([5.0e9, 5.0e9], marker=‘x’)

p.annotate_point([5.95e9, 5.1e9], ‘Star!’)

*# This saves the plot to a file with the given prefix. We can alternatively specify*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# the entire filename.*

p.save(‘contours.press_den\_’)

.. raw:: latex

   \centering

.. figure:: Slice_z_pressure
   :alt: Pressure slice with annotations
   :width: 6in

   Pressure slice with annotations

*#————————*

*# Volume Rendering*

.. raw:: latex

   \setlength{\parskip}{0pt}

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’
dd = pf.all_data()

*# We take the log of the extrema of the pressure field, as well as a couple other interesting*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# value ranges we’d like to visualize.*

h_mi, h_ma = dd.quantities[‘Extrema’](field)[0]

h_mi, h_ma = np.log10(h_mi), np.log10(h_ma)

s_mi, s_ma = np.log10(2.90e13), np.log10(3.10e13)

pf.index

*# We deal in terms of logarithms here because we have such a large range of values.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# It can make things easier, but is not necessary.*

pf.field_info[field].take_log=True

*# This is what we use to visualize volumes. There are a couple of other, more complex*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# ways. We set the range of values we’re interested in and the number of bins in the*

*# function. Make sure to have a lot of bins if your data spans many orders of magnitude!*

*# Our raw data ranges from about :math:`10^{13}` to :math:`10^{22}`.*

tf = ColorTransferFunction((h_mi-1, h_ma+1), nbins=1.0e6)

*# Here we add several layers to our function, either one at a time or in groups. We*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# specify the value-center and width of the layer. We can manipulate the color by*

*# individually setting the colormaps and ranges to spread them over. We can also*

*# change the transparency, which will usually take some time to get perfect.*

tf.sample_colormap(np.log10(2.0e21), 0.006, col_bounds=[h_mi,h_ma],

alpha=[27.0], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.0e19), 0.001, col_bounds=[h_mi,h_ma],

alpha=[5.5], colormap=‘RdBu_r’)

tf.add_layers(6, mi=np.log10(2.95e13), ma=s_ma,

col_bounds=[s_mi,s_ma],

alpha=19*na.ones(6,dtype=‘float64’), colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.95e13), 0.000005, col_bounds=[s_mi,s_ma],

alpha=[13.0], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.90e13), 0.000007, col_bounds=[s_mi,s_ma],

alpha=[11.5], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.85e13), 0.000008, col_bounds=[s_mi,s_ma],

alpha=[9.5], colormap=‘RdBu_r’)

*# By default each color channel is only opaque to itself. If we set grey_opacity=True,*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# this is no longer the case. This is good to use if we want to obscure the inner*

*# portions of our rendering. Here it only makes a minor change, as we must set our*

*# alpha values for the outer layers higher to see a strong effect.*

tf.grey_opacity=True

*# Volume rendering uses a camera object which centers the view at the coordinates we’ve*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# called ‘c.’ ‘L’ is the normal vector (automatically normalized) between the camera*

*# position and ‘c,’ and ‘W’ determines the width of the image—again, like a zoom.*

*# ‘Nvec’ is the number of pixels in the x and y directions, so it determines the actual*

*# size of the image.*

c = [5.0e9, 5.0e9, 5.0e9]

L = [0.15, 1.0, 0.40]

W = (pf.domain_right_edge - pf.domain_left_edge)*0.5

Nvec = 768

*# ‘no_ghost’ is an optimization option that can speed up calculations greatly, but can*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# also create artifacts at grid edges and affect smoothness. For our data, there is no*

*# speed difference, so we opt for a better-looking image.*

cam = pf.camera(c, L, W, (Nvec,Nvec), transfer_function = tf,

fields=[field], pf=pf, no_ghost=False)

*# Obtain an image! However, we’ll want to annotate it with some other things before*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# saving it.*

im = cam.snapshot()

*# Here we draw a box around our stars, and visualize the gridding of the top two levels.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Note that draw_grids returns a new image while draw_box does not. Also, add\_*

*# background_color in front of draw_box is necessary to make the box appear over*

*# blank space (draw_grids calls this internally). For draw_box we specify the left*

*# (lower) and right(upper) bounds as well its color and transparency.*

im.add_background_color(‘black’, inline=True)

cam.draw_box(im, np.array([3.0e9, 4.0e9, 4.0e9]),

np.array([7.0e9, 6.0e9, 6.0e9]), np.array([1.0, 1.0, 1.0, 0.14]))

im = cam.draw_grids(im, alpha=0.12, min_level=2)

im = cam.draw_grids(im, alpha=0.03, min_level=1, max_level=1)

*# ‘im’ is an image array rather than a plot object, so we save it using a different*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# function. There are others, such as ‘write_bitmap.’*

im.write_png(‘pressure_shell_volume.png’)

.. raw:: latex

   \centering

.. figure:: volume
   :alt: Volume rendering
   :width: 3.5in

   Volume rendering

*#————————*

*# Isocontour Rendering*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Here we extract isocontours using some extra modules and plot them using matplotlib.*

from mpl_toolkits.mplot3d import Axes3D

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import matplotlib.pyplot as plt

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’

field_weight = ‘magvel’

contour_value = 2.83e13

domain = pf.all_data()

*# This object identifies isocontours at a given value for a given field. It returns*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# the vertices of the triangles in that isocontour. It requires a data source, which*

*# can be an object—but here we just give it all of our data. Here we find a pressure*

*# isocontour and color it the magnitude of velocity over the same contour.*

surface = pf.surface(domain, field, contour_value)

colors = apply_colormap(np.log10(surface[field_weight]), cmap_name=‘RdBu’)

fig = plt.figure()

ax = fig.gca(projection=‘3d’)

p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

p3dc.set_facecolors(colors[0,:,:]/255.)

ax.add_collection(p3dc)

*# By setting the scaling on the plot to be the same in all directions (using the x scale),*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# we ensure that no warping or stretching of the data occurs.*

ax.auto_scale_xyz(surface.vertices[0,:], surface.vertices[0,:],

surface.vertices[0,:])

ax.set_aspect(1.0)

plt.savefig(‘pres_magvel_isocontours.png’)

.. raw:: latex

   \centering

.. figure:: isocontours
   :alt: Pressure isocontour rendering colored with velocity magnitude
   :width: 4in

   Pressure isocontour rendering colored with velocity magnitude

*#————————*

*#1-D and 2-D Profiles*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Line plots and phase plots can be useful for analyzing data in detail.*

from yt.mods import \*

pf = load(‘plt00020’)

pf.index

*# Just like with the pressure_contours script, we can set the units for fields that*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# have none.*

pf.field_info[‘magvel’]._units = r‘\\rm{cm}/\rm{s}’

pf.field_info[‘kineng’]._units = r‘\\rm{ergs}’

*# We can create new fields from existing ones. ytassumes all units are in cgs, and*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# does not do any unit conversions on its own (but we can make it). Creating new fields*

*# requires us to define a function that acts on our data and returns the new data,*

*# then call add_field while supplying the field name, the function the data comes from,*

*# and the units. Here, we create new fields simply to rename our data to make the plot*

*# look prettier.*

def \_newT(field, data):

return data[‘t’]

add_field(‘X’, function=_newT, units=r‘\\rm{domain} \rm{fraction}’)

def \_newDen(field, data):

return data[‘density’]

add_field(‘Density’, function=_newDen, units=r‘\\rm{g}/\rm{cm}^{3}’)

*# PlotCollections are one of the most commonly used tools in yt, alongside SlicePlots and*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# ProjectionPlots. They are useful when we want to create multiple plots from the same*

*# parameter file, linked by common characteristics such as the colormap, its bounds, and*

*# the image width. It is easy to create 1-D line plots and 2-D phase plots through a*

*# PlotCollection, but we can also create thin projections and so on. When we create a*

*# PlotCollection, it is empty, and only requires the parameter file and the ’center’ that*

*# will be supplied to plots like slices and sphere plots.*

pc = PlotCollection(pf, ‘c’)

*# Now we add a ray—a sample of our data field along a line between two points we define*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# in the function call.*

ray = pc.add_ray([0.0, 5.0e9, 5.0e9],[1.e10, 5.0e9, 5.0e9], ‘magvel’)

*# This is where our derived fields come in handy. Our ray is drawn along the x-axis*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# through the center of the domain, but by default the fraction of the ray we have gone*

*# along is called ‘t.’ We now have the same data in another field we called ‘X,’ whose*

*# name makes more sense, so we’ll reassign the ray’s first field to be that. If we wanted,*

(*# we could also reassign names to ‘magvel’ and ‘kineng.’*

ray.fields = [‘X’, ‘magvel’]

*# Next, we’ll create a phase plot. The function requires a data source, and we can’t*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# just hand it our parameter file, but as a substitute we can quickly create an object*

*# that spans our entire domain (or use the method in the isocontour example). The*

*# specifications of the region (a box) are the center, left bound, and right bound.*

region = pf.region([5.0e9, 5.0e9, 5.0e9], [0.0, 0.0, 0.0],

[1.0e10, 1.0e10, 1.0e10])

*# The phase object accepts a data source, fields, a weight, a number of bins along both*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# axes, and several other things, including its own colormap, logarithm options,*

*# normalization options, and an accumulation option. The first field is binned onto*

*# the x-axis, the second field is binned onto the y-axis, and the third field is*

*# binned with the colormap onto the other two. Subsequent fields go into an underlying*

*# profile and do not appear on the image.*

phase = pc.add_phase_object(region, [‘Density’, ‘magvel’,‘kineng’], weight=None,

x_bins=288, y_bins=288)

pc.save(‘profile’)

.. raw:: latex

   \centering

.. figure:: LineQueryPlot_0_t_magvel
   :alt: Density/velocity magnitude/kinetic energy phase plot
   :width: 4in

   Density/velocity magnitude/kinetic energy phase plot

.. figure:: Profile2D_1_Density_magvel_kineng
   :alt: Density/velocity magnitude/kinetic energy phase plot
   :width: 4in

   Density/velocity magnitude/kinetic energy phase plot

.. raw:: latex

   \quad

*#————————*

*#Off-Axis Projection*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# If we don’t want to take a projection (this can be done for a slice as well) along*

*# one of the coordinate axes, we can take one from any direction using an*

*# OffAxisProjectionPlot. To accomplish the task of setting the view up, the plot*

*# requires some of the same parameters as the camera object: a normal vector, center,*

*# width, and field, and optionally we can set no_ghost (default is False). The normal*

*# vector is automatically normalized as in the case of the camera. The plot also*

*# requires a depth—that is, how much data we want to sample along the line of sight,*

*# centered around the center. In this case ‘c’ is a shortcut for the domain center.*

pf = load(‘plt00020’)

field = ‘density’

L = [0.25, 0.9, 0.40]

plot = OffAxisProjectionPlot(pf, L, field, center=‘c’,

width=(5.0e9, 4.0e9), depth=3.0e9)

*# Here we customize our newly created plot, dictating the font, colormap, and title.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Logarithmic data is used by default for this plot, so we turn it off.*

plot.set_font({‘family’:‘Bitstream Vera Sans’, ‘style’:‘italic’,

‘weight’:‘normal’, ‘size’:14, ‘color’:‘red’})

plot.set_log(field, False)

plot.set_cmap(field, ‘jet’)

plot.annotate_title(‘Off-Axis Density Projection’)

*# The actual size of the image can also be set. Note that the units are in inches.*

.. raw:: latex

   \setlength{\parskip}{0pt}

plot.set_window_size(8.0)

plot.save(‘off_axis_density’)

.. raw:: latex

   \centering

.. figure:: OffAxisProjection_density
   :alt: Off-axis density projection
   :width: 4in

   Off-axis density projection
