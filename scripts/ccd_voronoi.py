#!/usr/bin/env python3
# Creates Periodic Voronoi tesellations from a given set of points
# This uses the Python module Freud, which in turn uses Voro++. Freud is multithreaded.
# Usage: ccd_voronoi.py <box> <pts_in_file_path> <cells_out_file_path>

import sys
import numpy as np
import freud

boxlen = float(sys.argv[1]) # Access box length from first command-line argument
pts_in = sys.argv[2] # Input file containing XY coordinates of the points
cells_out = sys.argv[3] # Output file containing Voronoi cell vertices

points = np.loadtxt(pts_in) # Loading the input points (2D)

# Periodic Voronoi construction using Freud
box = freud.box.Box.square(boxlen)
voro = freud.locality.Voronoi()
voro.compute((box, points))
cells = voro.polytopes

# Output Voronoi cell vertices. Cells are separated by a blank line.
f = open(cells_out, 'w')
for cell in cells:
    np.savetxt(f, cell)
    f.write("-\n") # Blank line to separate cells
f.close()

# Output neighborlist
for i, j in voro.nlist[:]:
    print(i+1, j+1) # +1 to change 0-based array to 1-based (because cells are numbered 1 to m) 
