# gpx_reduce_light

* gpx_reduce.py is a modified version of the [original gpx_reduce][1] with no dependencies.

    The original depends on *scipy, lxml, numpy* and *pylab* (which are heavy requirements).
    However, these dependencies can be easily removed, because the program only does basic linear
    algebra.
    This has two benefits:

    1. easy installation (nothing needed apart from python)
    2. the program is considerably faster than the original.
    
    A stand alone windows executable (win32, 3 Mb) can be downloaded from <https://wim.vree.org/sporen/gpx.html>

    Usage example (see [original][1] for a better description):

        > gpx_reduce.py -d 2 -t 30 your_track.gpx

    The disadvantage of removing all dependencies is that the plot option had to be removed.  
    I made a separate python script with one dependency for plotting tracks:

* gpx_plot.py, a script to plot one or more tracks with [gnuplot][3].

    Gnuplot has to be installed and if the path to the gnuplot executable is not */usr/bin/gnuplot*
    you have to specify the path with the command line option *-g /path/to/gnuplot*

    Usage example that compares a reduced track with the original:

        > gpx_plot.py your_track.gpx your_track_reduced.gpx
    
    and with the -g option:

        > gpx_plot.py -g /path/to/gnuplot your_track.gpx your_track_reduced.gpx
    
    A stand alone windows executable (win32, 3 Mb) can be downloaded from <https://wim.vree.org/sporen/gpx.html>
    (does not include [gnuplot][3])

[1]: https://github.com/Alezy80/gpx_reduce/
[3]: https://sourceforge.net/projects/gnuplot/files/gnuplot/
