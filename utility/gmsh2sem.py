#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Read a gmsh (http://www.geuz.org/gmsh/) mesh file (.msh extension
# typical but we do not check), assumed to supply an all-quad mesh,
# and from this, produce a semtex session file skeleton, containing
# NODES and ELEMENTS sections. Wrap this with other standard semtex
# feml file sections in order to generate a valid semtex session file.
#
# Print all this to standard output.
#
# Usage: gmsh2sem.py gmshfile.msh > session
# Usage: python gmsh2sem.py gmshfile.msh > session
#
# NB: use Mesh.SubdivisionAlgorithm=1; in gmsh in order to generate quads.
# Standard semtex utility meshpr must be found in your shell's PATH variable.
#
# $Id: gmsh2sem.py,v 8.1 2015/04/20 11:14:18 hmb Exp $
##############################################################################

import sys, os
import pdb
import numpy as np
import logging

class Node(object):
    """Node object, holds tag and position. Initialized from string."""
    def __init__(self, line):
        tag, x, y, z = line.split()
        self.tag = int(tag)
        self.X = np.array((float(x), float(y), float(z)))

    def __str__(self):
        #sc = 0.001 # dirty hack: scale nodes
        sc = 1.
        return "\t%i\t%g\t%g\t%g\n" % (self.tag, self.X[0]*sc, self.X[1]*sc, self.X[2]*sc)

class Element(object):
    def __init__(self, tag, nodes):
        self.tag = tag
        self.nodes = map(int, nodes)

    def __str__(self):
        return '\t% 4i <Q> % 4i % 4i % 4i % 4i </Q>\n' % (self.tag, \
            self.nodes[0], self.nodes[1], self.nodes[2], self.nodes[3])

class Version(object):
    def __init__(self, version_string):
        if version_string == '2.1':
            self.set_for_2_1()
        else:
            logging.warn('Unknown version %s, defaulting to 2.1' % version_string)
            self.set_for_2_1()
            
    def set_for_2_1(self):
        self.physidx_col = 3
        self.PhysNamesTag = 1
        

def find_element_and_edge(Elements, node0, node1):
    """find edge containing given nodes. Return element/edge number"""
    for e in Elements:
        if node0 in e.nodes and node1 in e.nodes:
            elmt = Elements.index(e)
            i0 = e.nodes.index(node0)
            i1 = e.nodes.index(node1)
            lo = min(i0, i1)
            hi = max(i0, i1)
            if hi == lo + 1:
                return elmt, lo
            else:
                return elmt, hi

    raise ValueError("nodes %i and %i are not on a common element." % (node0, node1))




# -- main

# -- "parse" command line
writeTokens = True
writeBC = True
if '-t' in sys.argv: writeTokens = False
if '-b' in sys.argv: writeBC = False

fields = ["u", "v", "w", "p"]
if '-f' in sys.argv:
    fieldstr = sys.argv[sys.argv.index('-f') + 1]
    fields = []
    for i in range(len(fieldstr)):
        fields.append(fieldstr[i])

try:
    infilename = sys.argv[-1]
except:
    print "Usage: gmsh2sem [options] infile.msh"; sys.exit(1)

ifile = open (infilename, 'r')
t1    = open (".t1", 'w')

genericTokens = """# Skeleton semtex session file generated from gmsh .msh file

<TOKENS>
        N_P    = 5
        N_STEP = 20
        D_T    = 0.01
        KINVIS = 1.
</TOKENS>

"""
genericGroup = """
<GROUPS NUMBER=1>
        1       w       wall
</GROUPS>
"""
genericBC = """
<BCS NUMBER=1>
        1       w       3
                <D> u = 0. </D>
                <D> v = 0. </D>
                <H> p      </H>
</BCS>
"""

genericHis = """
<HISTORY NUMBER=1>
    1 0 0 0
</HISTORY>
"""

# -- We will assume that the file is a valid gmsh mesh file with quad meshes.
#    First print a generic semtex header.

if writeTokens: t1.write(genericTokens)

Phys = {}
# -- Extract NODES and ELEMENTS from given gmsh file.

while 1:
    line = ifile.readline()

    if '$MeshFormat' in line:
        line = ifile.readline()
        version = Version(line.split()[0])

    # -- PhysicalNames are used to assign a name to lines/surfaces/volumes
    #    useful for assigning BCs
    if '$PhysicalNames' in line:
        line = ifile.readline()
        Nphys = int (line)
        for i in range(Nphys):
            line = ifile.readline()
            words = line.split()
            # -- PhysicalNames can contain 1/2/3-D items. We only
            #    care for 1-D here
            if int(words[0]) == version.PhysNamesTag:
                tag = words[2][1:-1]   # cut leading/trailing "
                key = int(words[1])
                Phys[key] = tag

    if '$Nodes' in line:
        # -- Store node objects in list, write to t1 later on.
        Nodes = []
        line = ifile.readline()
        Nnode = int (line)
        for i in range (Nnode):
            line = ifile.readline()
            Nodes.append(Node(line))
    elif '$Elements' in line:
        # -- Since we don't know the number of valid (i.e. quad) elements yet,
        #    we store them in a list first and write to t1 later on.
        Elements = []
        rawsurf = []
        line = ifile.readline()
        Nel = int (line)
        Nquads = 0
        for i in range (Nel):
            line = ifile.readline()
            words = line.split()
            if int(words[1]) == 3 and (len(words) == 10 or len(words) == 9):
                # -- Process quad.
                #    Quad element lines have 10 numbers (9 in gmsh2.6).
                #    The first one is a tag, and the last 4 are node numbers.
                Nquads += 1
                element = Element(Nquads, (words[-4], words[-3], words[-2], words[-1]))
                Elements.append(element)
            elif int(words[1]) == 1:
                # -- process physical lines, i.e. BC
                #    For now, we merely cache all lines,
                #    as we don't know of all quads yet
                rawsurf.append(line)

        break

ifile.close()
# -- check for inverted elements (node order not CCW)

if 1:
    for el in Elements:
        el_nodes = np.array( [Nodes[n-1].X for n in el.nodes] )
        center = el_nodes.sum(axis=0)/4.
        vec = el_nodes - center
        # -- compare, swap
        if np.cross(vec[0], vec[1])[2] < 0:
            sys.stderr.write("fixing inverted element %i\n" % el.tag)
            el.nodes = el.nodes[::-1]

# -- FIELDS section

t1.write("<FIELDS NUMBER=%i>\n\t" % len(fields))
for f in fields:
    t1.write('%s ' % f)
t1.write("\n</FIELDS>\n\n")


# -- If gmsh file contains PhysicalNames, we can generate
#    appropriate GROUPS and BCS section from it. Otherwise,
#    we write generic sections.

if len(Phys):

    # -- GROUPS section

    t1.write("<GROUPS NUMBER=%i>\n" % len(Phys))
    i = 0
    for p in Phys:
        i += 1
        t1.write("\t%i\t%s\t%s\n" % (i, Phys[p][0], Phys[p]))
        #    1       w       wall
    t1.write("</GROUPS>\n\n")

    # -- BCS section

    if writeBC:
        i = 0
        t1.write("<BCS NUMBER=%i>\n" % len(Phys))
        for p in Phys:
            i += 1
            t1.write('\t%i\t%s\t%i\n' % (i, Phys[p][0], len(fields)))
            if Phys[p][0] == 'a':
                # -- axis group
                for f in fields:
                    if f == 'p':
                        t1.write('\t\t<A> p      </A>\n')
                    else:
                        t1.write('\t\t<A> %s = 0. </A>\n' % f)
            else:
                # -- other
                for f in fields:
                    if f == 'p':
                        t1.write('\t\t<H> p      </H>\n')
                    else:
                        t1.write('\t\t<D> %s = 0. </D>\n' % f)
        t1.write("</BCS>\n\n")
else:
    t1.write(genericGroup)
    if writeBC: t1.write(genericBC)

t1.write(genericHis+"\n")
# -- NODES section

t1.write("<NODES NUMBER=%i>\n" % len(Nodes))
for n in Nodes:
    t1.write(str(n))
t1.write('</NODES>\n\n')

# -- ELEMENTS section

t1.write("<ELEMENTS NUMBER=%i>\n" % len(Elements))
for e in Elements:
    t1.write(str(e))
t1.write('</ELEMENTS>\n\n')

t1.close()

if len(Phys):
    # -- SURFACES section from read BC
    #    Now that we know all quads, we can match BC to surfaces
    #    FIXME: surfaces list is unsorted
    surflines = ""
    i = 0
    for line in rawsurf:
        words = line.split()
        physidx = int(words[version.physidx_col])
        try:
            elmt, edge = find_element_and_edge(Elements, int(words[-2]), int(words[-1]))
            i += 1
            surflines += "\t% 4i % 4i %i <B> %s </B>\n" % (i, elmt+1, edge+1, Phys[physidx][0])
        except ValueError:
            sys.stderr.write("ignoring unmatched nodes %s" % line)

    t2 = open (".t2", 'w')
    t2.write("<SURFACES NUMBER=%i>\n" % i)
    t2.write(surflines)
    t2.write("</SURFACES>\n")
    t2.close()
else:

    # -- Generate a SURFACES section that contains all the unfilled
    #    surface valencies (those element edges that do not mate to
    #    another element). You can then edit this (by hand, or other
    #    means).

    os.system("meshpr -s .t1 > .t2")

t1 = open (".t1", 'r')
t2 = open (".t2", 'r')

# -- Print up all of t1 and t2

for line in t1.readlines(): print line,
for line in t2.readlines(): print line,

t1.close();       t2.close()
os.remove('.t1'); os.remove('.t2')
