SparseGrids - a set of MATLAB/Octave scripts and functions
for implementing spare grid finite element methods.

Authors:  Stephen Russell and Niall Madden, NUI Galway. 

Date:     September 2016.

Download: https://github.com/niallmadden/SparseGrids

This is a set of 5 MATLAB/Octave scripts and function files for solving a 
linear, two-dimensional partial differential equation on the unit square,
using a standard Galerkin FEM with bilinear elements, and two sparse grid methods.

This code is used to generate results in Russell, S., and Madden, N. 
An analysis and implementation of sparse grid finite element methods. 
http://arxiv.org/abs/1511.07193

The files are: Test_FEM.m, FEM_System_Matrix.m, FEM_RHS.m, TwoScale_Projector.m, MultiScale_Projector.m

MIT License

Copyright (c) 2016 Niall Madden

Permission is hereby granted, free of charge, to any person obtaining a 
copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE