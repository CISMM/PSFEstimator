Last update: 2012-10-19
Author: Cory Quammen <cquammen@cs.unc.edu>

This document gives instructions on how to build PSF Estimator from source.


PRELIMINARIES

You'll need to obtain CMake, Qt, VTK, and ITK.

1. Download CMake 2.8.5 or higher from
http://www.cmake.org/cmake/resources/software.html

2. Download Qt libraries and applications by going to
http://qt-project.org/downloads and choosing Qt SDK appropriate for
the platform on which you are developing.

3. Download the VTK source code. The best way to do this is to access it from 
the git repository with

git clone git://vtk.org/VTK
git checkout 8ab89ff701ee6dd7d256398867bd0931c26da1ac

Configure VTK with CMake and build it. You'll need to change these options
from the default:

VTK_Group_Imaging=ON
Module_vtkGUISupportQt=ON
QT_QMAKE_EXECUTABLE=<path to the qmake executable that came with the Qt SDK you downloaded/built

4. Download the ITK source code. The best way to do this is to access it from
the git repository with

git clone git://itk.org/ITK
git checkout 9dbde859b0062f0b888129baf264537e6028bca8

5. Configure the PSF Estimator project with CMake. Set the following variables:

VTK_DIR=<path to the VTK build directory you specified in step 3>
ITK_DIR=<path to the ITK build directory you specified in step 4>


