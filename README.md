# gpx_reduce_light

* gpx_reduce.py, a version of gpx_reduce with no dependencies (only the standard python libraries)
Only the plot option has been removed.
See [original version](https://github.com/Alezy80/gpx_reduce) for more info.

Usage example:
```
> gpx_reduce.py -d 2 -t 30 your_track.gpx
```

* gpx_plot.py, a separate program to plot one or more tracks with "gnuplot".
Gnuplot has to be installed and the path to the binary has to be changed in the code
to reflect your installation:

python```
# the path to the gnuplot binary
gnuPlotCmd = 'path/to/gnuplot'
```

Usage example that compares a reduced track with the original:
```
> gpx_plot.py your_track.gpx your_track_reduced.gpx
```
