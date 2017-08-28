#!/usr/bin/env python
# -*- coding: utf8 -*-

'''
gpx_plot.py calls gnuplot to draw the gpx-tracks given as arguments.
Copyright 2017 willem179
Extracted from the code of "gpx_reduce.py"
Copyright (C) 2011,2012,2013,2015,2016,2017 travelling_salesman on OpenStreetMap

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import datetime
import sys
import time
from math import *
import xml.etree.ElementTree as etree
from optparse import OptionParser

# the path to the gnuplot binary
gnuPlotCmd = 'gnuplot'

parser = OptionParser('usage: %prog [options] input-file.gpx')
(options, args) = parser.parse_args()


if len(args) < 1:
    parser.print_usage()
    exit(2)


# use the WGS-84 ellipsoid
rE = 6356752.314245 # earth's radius
a = 6378137.0
b = 6356752.314245179

timeformat = '%Y-%m-%dT%H:%M:%SZ'

# the linear algebra with lists
norm = lambda p: sqrt (sum (a * a for a in p))

def rotate(x, y, phi):
    return x*cos(phi) - y*sin(phi), x*sin(phi) + y*cos(phi)


def project_to_meters(lat, lon, latm, lonm):
    # azimuthal map projection centered at average track coordinate
    lon -= lonm
    xyz = latlonele_to_xyz(lat, lon, 0.0)
    zy = rotate(xyz[2], xyz[1], radians(90 - latm))
    lat2 = atan2(zy[0], norm([zy[1], xyz[0]]))
    lon2 = atan2(xyz[0], -zy[1])
    x_meters = rE * sin(lon2) * (pi / 2.0 - lat2)
    y_meters = -rE * cos(lon2) * (pi / 2.0 - lat2)
    return x_meters, y_meters


def latlonele_to_xyz(lat, lon, ele):
    s = sin(radians(lat))
    c = cos(radians(lat))
    r = ele + a * b / norm([s*a, c*b])
    lon = radians(lon)
    return r * c * sin(lon), r * c * (-cos(lon)), r * s


def xyz_to_latlonele(x, y, z):
    r = norm([x, y, z])
    if (r == 0):
        return 0.0, 0.0, 0.0
    lat = degrees(atan2(z, norm([x, y])))
    lon = degrees(atan2(x, -y))
    ele = r * (1.0 - a * b / norm([a*z, b*x, b*y]))
    return lat, lon, ele


############################## main function #################################
tracks = []
npoints = []
for fname in args:
    # initialisations
    tracksegs = []
    sumx, sumy, sumz = 0.0, 0.0, 0.0
    ntot = 0    # total number of trackpoints (sum of segments)
    
    # import xml data from files
    print 'opening file', fname
    infile = open(fname)
    tree = etree.parse(infile)
    infile.close()

    gpx = tree.getroot()
    nsurl = gpx.tag.split ('}')[0][1:]  # == 'http://www.topografix.com/GPX/1/1'
    etree.register_namespace('', nsurl) # default namespace -> xmlns:.... in the output
    nsmap = '{' + nsurl + '}'           # prefix for all tags in the tree

    # extract data from xml
    for si, trkseg in enumerate (gpx.findall('.//' + nsmap + 'trkseg')):
        trkpts = trkseg.findall(nsmap + 'trkpt')
        n = len(trkpts)
        
        # extract coordinate values
        lats = [float(trkpt.get('lat')) for trkpt in trkpts]
        lons = [float(trkpt.get('lon')) for trkpt in trkpts]
        eles = [float(trkpt.find(nsmap + 'ele').text) for trkpt in trkpts]
        try:
            times = [datetime.datetime.strptime(trkpt.find(nsmap + 'time'
                                       ).text, timeformat) for trkpt in trkpts]
        except Exception as e:
            print '-- trackpoint without time'
            times = None

        # save original trackseg for plotting
        tracksegs.append ([ (lats[i], lons[i], eles[i]) for i in range(n) ])

        # calculate projected points to work on
        for i in range(n):
            x, y, z = latlonele_to_xyz (lats[i], lons[i], eles[i])
            sumx += x
            sumy += y
            sumz += z

        print 'segment %d, with %d points:' % (si, n)
        ntot += n

    tracks.append (tracksegs)
    npoints.append ((fname.replace ('_','\_'), ntot))   # underscore -> subscripts in gnuplot
    print 'number of points in track %s = %d:' % (fname, ntot)

latm, lonm, elesum = xyz_to_latlonele (sumx, sumy, sumz)

data = []
xmin, xmax = ymin, ymax = float ('inf'), float ('-inf')
for tracksegs in tracks:
    for trkseg in tracksegs:
        for lat, lon, ele in trkseg:
            x, y = project_to_meters (lat, lon, latm, lonm)
            data.append ('%f %f' % (x, y))
            xmin, ymin = min (xmin, x), min (ymin, y)   # determine the x range
            xmax, ymax = max (xmax, x), max (ymax, y)   # and the y range
    data.append ('e')

dx, dy = xmax - xmin, ymax - ymin   # make x and y ranges equal to the largest
if dx > dy:                         # and keep ranges centered
    dr = (dx - dy) / 2
    ymax += dr
    ymin -= dr
else:
    dr = (dy - dx) / 2
    xmax += dr
    xmin -= dr

from subprocess import Popen, PIPE
plot = Popen ([gnuPlotCmd], stdin=PIPE, stdout=PIPE, stderr=PIPE)

range = 'set xrange [%f:%f]\nset yrange [%f:%f]\n' % (xmin, xmax, ymin, ymax)
plot.stdin.write (range)

curves = ','.join ("'-' with linespoints ti '%s: %d punten'" % t for t in npoints)
plot.stdin.write ('plot ' + curves + '\n')

plot.stdin.write ("\n".join (data))
plot.stdin.write ('\n')
plot.stdin.flush ()
raw_input ('druk')
