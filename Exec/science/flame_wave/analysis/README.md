# Flame Wave Analysis Scripts

### plot_generator.py

Script for generating plots of a sequence of datasets using yt. The variable to plot, whether to use
logscale, the domain and colorbar bounds, streamline options, etc. are all configurable through
command line arguments. For a usage description and a full list of valid parameters, type
`./plot_generator.py -h`. *TODO*: Add output dpi and figure size settings.

### image_animator.py

Uses ffmpeg directly (through an `os.system` call) to generate an animation from a sequence of
images. Allows the user to specify an output file and offers some support for sorting the images.
Advantages: Fast, low memory footprint, configurable framerate. Restrictions: The width and height
of each image must be divisible by 2.

### mpl_image_animator.py

Uses matplotlib to generate an animation from a sequence of images. Allows the user to specify an
output file and offers some support for sorting the images. Advantages: No height and width
restriction, offers a stack feature for vertically stacking multiple images using separate subplots.
Restrictions: Much slower than the ffmpeg version, and much larger memory footprint. *TODO*: Allow
user to specify framerate.

### front_tracker.py
Script for measuring the location of a flame or shock front over a sequence of snapshots. Allows the
user to specify the metrics (will only track 1 / 1000th the enuc maximum by default), domain bounds,
and a few other things. For a usage description and a full list of valid parameters, type
`./front_tracker.py -h`. Should work for any dataset, but has only been tested on flame wave ones.
Restrictions: Currently only tracks along one dimension (the user can tell it how to eliminate the
others - either through slicing or averaging), and only tracks percentages of the maximum of some
field. Outputs to space-delimited data file called front_tracking.dat by default.

### flame_speed.py
Script for reading in the front tracking dataset generated with `front_tracker.py`, plotting it, and
fitting a line to some portion of it.
Usage: `./flame_speed.py [--tmin TMIN] [--tmax TMAX] data_file column [column...]`, where `TMIN` and
`TMAX` specify the times to consider when fitting the line. The script prints out the slope of the
line, the r-squared value, and the fit error. Uses `scipy` and `pandas`.

### multirays.py
Plot 1D vertical slices of axisymmetric datasets. It generates 3 ortho rays - one at each end of
the domain and one along the center. The variable to plot can be supplied as a command line argument
(e.g. `./multirays.py -v Temp`).

### parallel
Run multiple instances of a particular analysis script in parallel. Can only run one instance of
this at a time on a given directory.

### overview.py
Plots temperature, enuc, and z-velocity in a vertical stack. The fields to plot can be set by
modifying the fields list in script.

### time_series.py
Create a stacked plot of abar at a sequence of time points.

### schlieren.py
Make a Schlieren plot (a plot of ln{(∇^2 ρ) / ρ}) of a dataset.
