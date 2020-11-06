#!/usr/bin/env python

from fipy.meshes import Gmsh3D 

cfgfile = 'parameters'
from parameters import *

def macroMeshSphere(rad, dx):
	print rad, dx
        MESHstr = """
	rad = %(rad)g;
	cl1 = %(dx)g;
	Point(1) = {0, 0, 0, cl1};
	Point(2) = {rad, 0, 0, cl1};	
	Point(3) = {0, rad, 0, cl1};
	Point(4) = {0, 0, rad, cl1};
	Point(5) = {-rad, 0, 0, cl1};
	Point(6) = {0, -rad, 0, cl1};
	Point(7) = {0, 0, -rad, cl1};
	Circle(1) = {2, 1, 7};
	Circle(2) = {7, 1, 5};
	Circle(3) = {5, 1, 4};
	Circle(4) = {4, 1, 2};
	Circle(5) = {2, 1, 3};
	Circle(6) = {3, 1, 5};
	Circle(7) = {5, 1, 6};
	Circle(8) = {6, 1, 2};
	Circle(9) = {7, 1, 3};
	Circle(10) = {3, 1, 4};
	Circle(11) = {4, 1, 6};
	Circle(12) = {6, 1, 7};
	Line Loop(14) = {5, 10, 4};
	Ruled Surface(14) = {14};
	Line Loop(16) = {9, -5, 1};	
	Ruled Surface(16) = {16};
	Line Loop(18) = {12, -1, -8};
	Ruled Surface(18) = {18};
	Line Loop(20) = {8, -4, 11};
	Ruled Surface(20) = {20};
	Line Loop(22) = {-10, 6, 3};
	Ruled Surface(22) = {22};
	Line Loop(24) = {-11, -3, 7};
	Ruled Surface(24) = {24};
	Line Loop(26) = {-2, -12, -7};
	Ruled Surface(26) = {26};
	Line Loop(28) = {-6, -9, 2};
	Ruled Surface(28) = {28};
	Surface Loop(30) = {28, 22, 14, 16, 18, 26, 24, 20};
	Volume(30) = {30};"""
        mesh = Gmsh3D(MESHstr % locals())
        meshstyle = 'sph'
        return mesh, meshstyle

def macroMeshRect(dx,L):
    MESHstr = """
	dx = %(dx)g;
	L = %(L)g;
	Point(7) = {0, 0, 0, dx};
	Point(16) = {L, 0, 0, dx};	
	Point(24) = {0, L, 0, dx};
	Point(28) = {L, L, 0, dx};
	Point(104) = {0, 0, L, dx};
	Point(105) = {0, L, L, dx};
	Point(109) = {L, L, L, dx};
	Point(113) = {L, 0, L, dx};
	Line(105) = {16, 7};
	Line(107) = {7, 24};
	Line(108) = {24, 28};
	Line(109) = {28, 16};
	Line(131) = {104, 105};
	Line(132) = {105, 109};
	Line(133) = {109, 113};
	Line(134) = {113, 104};
	Line(136) = {7, 104};
	Line(137) = {24, 105};
	Line(141) = {28, 109};
	Line(145) = {16, 113};
	Line Loop(106) = {107, 108, 109, 105};
	Plane Surface(106) = {106};
	Line Loop(138) = {107, 137, -131, -136};	
	Ruled Surface(138) = {138};
	Line Loop(142) = {108, 141, -132, -137};
	Ruled Surface(142) = {142};
	Line Loop(146) = {109, 145, -133, -141};
	Ruled Surface(146) = {146};
	Line Loop(150) = {105, 136, -134, -145};
	Ruled Surface(150) = {150};
	Line Loop(151) = {131, 132, 133, 134};
	Plane Surface(151) = {151};
	Surface Loop(130) = {106, 151, 138, 142, 146, 150};
	Volume(130) = {130};
	Physical Volume(1) = {130};"""
    mesh = Gmsh3D(MESHstr % locals())
    meshstyle = 'rect'
    return mesh, meshstyle

