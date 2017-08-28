#!/usr/bin/env python
# -*- coding: utf8 -*-

'''
gpx_reduce v1.8: removes points from gpx-files to reduce filesize and
tries to keep introduced distortions to the track at a minimum.
Copyright (C) 2011,2012,2013,2015,2016,2017 travelling_salesman on OpenStreetMap

changelog: v1.2: clarity refractoring + speedup for identical points
           v1.3: new track weighting functions, progress display
           v1.4: new track weighting function, restructuring for memory saving
           v1.5: algorithm speedup by roughly a factor of 2 by eliminating some cases.
           v1.6: presets for train etc.
           v1.7: introduced weighting function for elevation errors
           v1.8: speed-dependent distance limit
           v1.9: removal of all dependencies (lxml, scipy, numpy, pylab)

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

import datetime, sys, time, operator
from math import *
import xml.etree.ElementTree as etree
from optparse import OptionParser


parser = OptionParser('usage: %prog [options] input-file.gpx')
parser.add_option('-v', '--verbose', action='store', type='int',
    dest='verbose', default=1, help='verbose=[0,1]')
parser.add_option('-d', '--dist', action='store', type='float', dest='max_dist',
    default=0.5, help='Maximum distance of line from original points in meters')
parser.add_option('-o', '--out', action='store', type='string',
    dest='ofname', default=None, help='Output file name')
parser.add_option('-m', '--maxsep', action='store', type='float', dest='max_sep',
    default=200.0, help='Absolute maximum separation of points. No points will be deleted where the resulting distance would become greater than maxsep. Standard JOSM settings will not display points spaced more than 200m. -1 means no limit.')
parser.add_option('-n', '--maxsep0', action='store', type='float', dest='max_sep0',
    default=4.0, help='Maximum separation of points at zero speed.')
parser.add_option('-t', '--maxseptime', action='store', type='float', dest='max_sep_time',
    default=3.0, help='Maximum time separation of points, which will be added to maxsep0. Maximum allowed point separation will be min(m, n + t*v) where v is the speed in m/s.')
parser.add_option('-e', '--ele-weight', action='store', type='float',
    dest='ele_weight', default=0.0,
    help='Weighting of elevation errors vs. horizontal errors. Default is 0.')
parser.add_option('-b', '--bend', action='store', type='float', dest='bend',
    default=1.0, help='Penalty value for large bending angles at each trackpoint. Larger values (1 or more) make the track smoother.')
parser.add_option('-w', '--weighting', action='store', type='string',
    dest='weighting', default='exp',
    help='''Weighting function to be minimized for track reduction:
pnum (number of points),
sqrdistsum (number of points plus sum of squared distances to leftout points),
sqrdistmax (number of points plus sum of squared distances to each maximally separated leftout point per new line segment),
sqrlength (number of points plus sum of squared new line segment lengths normalized by maxsep),
mix (number of points plus sum of squared distances to each maximally separated leftout point per new line segment weighted with corresponding segment length),
exp (number of points plus sum of squared distances to leftout points with exponential weighting of 1/2, 1/4, 1/8... from furthest to closest point). exp=standard''')
parser.add_option ('-r', dest='remove', default='extensions,hdop', help='remove tags T1,T2,..,Tn from every trackpoint', metavar='T1 T2 Tn')

options, args = parser.parse_args()

if len(args) < 1:
    parser.print_help ()
    sys.exit (2)

# use the WGS-84 ellipsoid
rE = 6356752.314245 # earth's radius
a = 6378137.0
b = 6356752.314245179

timeformat = '%Y-%m-%dT%H:%M:%SZ'

# the linear algebra with lists
dot = lambda p1, p2: reduce (operator.add, [a * b for a, b in zip (p1, p2)])
vmin  = lambda p1, p2: [a - b for a, b in zip (p1, p2)]
smul = lambda p, c: [a * c for a in p]
sdiv = lambda p, c: [a / c for a in p]
norm = lambda p: sqrt (sum (a * a for a in p))

def distance (p1, pm, p2, ele_weight=1.0):
    # returns distance of pm from line between p1 and p2

    h1, hm, h2 = norm (p1), norm (pm), norm (p2)
    if ele_weight != 1.0 and min(h1, hm, h2) > 0.0:
        hmean = (h1 + hm + h2) / 3.0
        p1 = smul (p1, ele_weight + (1.0 - ele_weight) * hmean / h1)
        pm = smul (pm, ele_weight + (1.0 - ele_weight) * hmean / hm)
        p2 = smul (p2, ele_weight + (1.0 - ele_weight) * hmean / h2)
    line = vmin (p2, p1)
    linel = norm (line)
    vm = vmin (pm, p1)
    if linel == 0.0:
        return norm (vm)
    linem = sdiv (line, linel)
    
    position = dot (vm, linem) / linel
    if position < 0.0:
        return norm (vm)
    elif position > 1.0:
        return norm (vmin (pm, p2))
    else:
        return norm (vmin (vm, smul (line, position)))

def latlonele_to_xyz (lat, lon, ele):
    s = sin (radians (lat))
    c = cos (radians (lat))
    r = ele + a * b / norm ([s*a, c*b])
    lon = radians (lon)
    return r * c * sin (lon), r * c * (-cos (lon)), r * s


def reduced_track_indices(coordinate_list, timesteps=None):
    # returns a list of indices of trackpoints that constitute the reduced track
    # takes a list of kartesian coordinate tuples
    m = len(coordinate_list)
    if (m == 0): return []
    if timesteps != None and len(timesteps) != len(coordinate_list):
        timesteps = None
    
    # number of dimensions
    d = len(coordinate_list[0])
    
    # remove identical entries (can speed up algorithm considerably)
    original_indices = [0]
    points = [{'p': coordinate_list[0], 'weight':1}]
    if timesteps != None: points[0]['t'] = timesteps[0]
    for i in range(1, m):
        if False in [coordinate_list[i-1][j] == coordinate_list[i][j] for j in range(d)]:
            original_indices.append(i)
            points.append({'p': coordinate_list[i], 'weight':1})
            if timesteps != None: points[-1]['t'] = timesteps[i]
        else:
            points[-1]['weight'] += 1
    n = len(points)
    
    # progress printing initialisations
    progress_printed = False
    progress = None
    tprint = time.time()
    
    # execute Dijkstra-like algorithm on points
    points[0]['cost'] = 1.0
    points[0]['prev'] = -1
    
    for i2 in range(1, n):
        penalties = {}
        imin = None
        costmin = float('inf')
        for i1 in reversed(range(i2)):
            p1 = points[i1]['p']
            p2 = points[i2]['p']
            seglength = norm (vmin (p2, p1))
            
            # estimate speed between p1 and p2
            if timesteps != None:
                dt = (points[i2]['t'] - points[i1]['t']).total_seconds()
                v = seglength / max(0.1, dt)
            else:
                v = seglength / float(i2 - i1) # assume 1s time spacing
            
            max_sep = options.max_sep0 + v * options.max_sep_time
            if options.max_dist >= 0:
                max_sep = min(max_sep, options.max_sep)
            
            if (seglength >= max_sep and i1 != i2 - 1):
                # point separation is too far
                # but always accept direct predecessor i1 = i2 - 1
                if (seglength >= max_sep + options.max_dist):
                    # no chance to find a valid earlier predecessor point
                    break
                else:
                    continue
            
            if points[i1]['cost'] + 1.0 > costmin:
                # the possible predecessor i1 is already too bad.
                continue
            
            i1_i2_segment_valid = True
            lower_i1_possible = True
            distance_squaremax = 0.0
            distance_squaresum = 0.0
            distances_squared = []
            # iterate all medium points between i1 and i2
            for im in range(i1+1, i2):
                pm = points[im]['p']
                d = distance (p1, pm, p2, options.ele_weight)
                if (d <= options.max_dist):
                    d_sq = (d / options.max_dist) ** 2
                    distance_squaremax = max(distance_squaremax, d_sq)
                    distance_squaresum += points[im]['weight'] * d_sq
                    distances_squared.append(d_sq)
                else:
                    i1_i2_segment_valid = False
                
                    # check if connection to any further point i1 is impossible
                    d1 = dot (vmin (p1, p2), vmin (p1, p2))
                    d2 = dot (vmin (pm, p2), vmin (pm, p2))
                    dd = options.max_dist ** 2
                    d1d2 = dot (vmin (p1, p2), vmin (pm, p2))
                    # formula from cosines of point separation angle and cone-opening angles around points
                    if (d1 > dd and d2 > dd and (d1d2 + dd)**2 < (d2 - dd) * (d1 - dd)):
                        lower_i1_possible = False
                        break
            
            if (lower_i1_possible == False):
                break
            
            if i1_i2_segment_valid:
                if options.weighting == 'sqrdistmax':
                    penalties[i1] = distance_squaremax
                elif options.weighting == 'sqrdistsum':
                    penalties[i1] = distance_squaresum
                elif options.weighting == 'sqrlength':
                    penalties[i1] = (seglength / max_sep) ** 2
                elif options.weighting == 'mix':
                    penalties[i1] = (distance_squaremax * (1.0 + seglength / max_sep))
                elif options.weighting == 'exp':
                    penalties[i1] = 0.5 * sum([0.5**i * d for i, d in
                        enumerate(sorted(distances_squared, reverse=True))])
                else:
                    penalties[i1] = 0.0
                
                # add a penalty for kinks
                if options.bend > 0.:
                    if points[i1]['prev'] != -1:
                        p0 = points [points [i1]['prev']] ['p']
                        v0 = vmin (p1, p0)
                        v1 = vmin (p2, p1)
                        if norm (v0) > 0. and norm (v1) > 0.:
                            v0 = sdiv (v0, norm (v0))
                            v1 = sdiv (v1, norm (v1))
                            kink = (1.0 - dot (v0, v1)) / 2.0
                            penalties[i1] += options.bend * kink
        
        # find best predecessor
        imin = None
        costmin = float('inf')
        for prev, penalty in penalties.iteritems():
            # cost function is sum of points used (1.0) plus penalties
            cost = points[prev]['cost'] + 1.0 + penalty
            if cost < costmin:
                imin = prev
                costmin = cost
        points[i2]['cost'] = costmin
        points[i2]['prev'] = imin
        
        # print progess
        if options.verbose == 1 and (100 * i2) / n > progress and time.time() >= tprint + 1:
            tprint = time.time()
            progress = (100 * i2) / n
            print '\r', progress, '%',
            sys.stdout.flush()
            progress_printed = True
    
    if progress_printed:
        print '\r',
    
    # trace route backwards to collect final points
    final_pnums = []
    i = n-1
    while i >= 0:
        final_pnums = [i] + final_pnums
        i = points[i]['prev']
    
    return [original_indices[i] for i in final_pnums]



############################## main function #################################
for fname in args:
    # initialisations
    ntot = 0    # total number of trackpoints (sum of segments)
    newtot = 0 # total number of trackpoints after removal
    
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

        # calculate projected points to work on
        coords = []
        for i in range(n):
            x, y, z = latlonele_to_xyz (lats[i], lons[i], eles[i])
            coords.append((x, y, z))
        
        # execute the reduction algorithm
        final_pnums = reduced_track_indices(coords, times)

        n_new = len (final_pnums)
        print 'segment %d, with %d - %d = %d points' % (si, n, n - n_new, n_new)
        ntot += n
        newtot += n_new

        # delete certain points from original data
        delete_pnums = [i for i in range(n) if i not in final_pnums]
        for i in reversed(delete_pnums):
            trkseg.remove (trkpts[i])   # remove from the xml-tree
            del trkpts [i]              # also remove from the list

        # remove certain sub-elements from remaining points
        options.remove = options.remove.replace ('+','extensions,hdop,')
        taglist = options.remove.split (',')
        if taglist: print 'remove %s subelements from points' % ', '.join (taglist)
        for trkpnt in trkpts:
            for tag in taglist:
                e = trkpnt.find (nsmap + tag)
                if e != None: trkpnt.remove (e)

    print 'total number of points: %d - %d = %d' % (ntot, ntot - newtot, newtot)

    # export data to file
    if options.ofname != None:
        ofname = options.ofname
    elif fname.endswith('.gpx'):
        ofname = fname[:-4] + '_reduced.gpx'
    else:
        ofname = fname + '_reduced.gpx'
    outfile = open(ofname, 'w')
    tree.write (outfile, encoding="utf-8", xml_declaration=True, default_namespace=None, method="xml")
    outfile.close()
    print 'modified copy written to', ofname
