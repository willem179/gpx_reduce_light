# gpx_reduce_light

* gpx_reduce.py is a version of the [original gpx_reduce][1] with no dependencies.

    The [original version][1] depends on scipy, lxml, numpy and pylab.
    Because the program only does some basic linear algebra, all these dependencies can be easily removed.
    This has two benefits:

    1. easy installation (nothing needed apart from python)
    2. the program is considerably faster than the original.  
        This is because the original uses a numpy array for each trackpoint which
        incurs a large overhead for creating and for all subsequent little computations.

    The one disadvantage of removing all depecencies is that the plot option had to be removed.  
    Usage example:

        > gpx_reduce.py -d 2 -t 30 your_track.gpx

    I made a separate python script with one dependency for plotting tracks:

* gpx_plot.py, a script to plot one or more tracks with "gnuplot".

    Gnuplot has to be installed and the path to the binary has to be changed in the code
    to reflect your installation:

        gnuPlotCmd = 'path/to/gnuplot'

    Usage example that compares a reduced track with the original:

        > gpx_plot.py your_track.gpx your_track_reduced.gpx

[1]: https://github.com/Alezy80/gpx_reduce/