Last update: 2009-07-05
Author: Cory Quammen <cquammen@cs.unc.edu>

This document gives instructions on how to create a new Qt/VTK/ITK-based
application from this generic example.


PRELIMINARIES

You'll need to obtain CMake, Qt, VTK, and ITK.

1. Download CMake 2.6 or higher from
http://www.cmake.org/cmake/resources/software.html

2. Download Qt libraries and applications by going to
http://www.qtsoftware.com/downloads/downloads and choosing the LGPL/Free tab.
Download the Qt SDK appropriate for the platform on which you are developing.

3. Download the VTK source code. The best way to do this is to access it from 
the CVS repository. Instructions for doing so are at 
http://www.vtk.org/VTK/resources/software.html#cvs. Configure VTK with CMake
and build it. You'll need to enable GUI support (option VTK_USE_GUI_SUPPORT),
Qt (option VTK_USE_QVTK), and Qt version 4 (set option DESIRED_QT_VERSION to 4).

4. Download the ITK source code at
http://www.itk.org/ITK/resources/software.html. Version 3.14 is known to work 
with this example. Choose either InsightToolkit*.tar.gz or InsightToolkit*.zip 
depending on what extraction applications you have available on your development
system. Configure ITK with CMake and build it. No non-default options are 
required for this application example to work.


CREATING A NEW PROJECT

Copy this project root directory to a new location in your filesystem where you 
will be able to import it into CVS. Make sure you delete the CVS directories in
this project. A one-liner for this task in linux is

  find -type d -wholename '*CVS' | xargs /bin/rm -rf

Next, you'll need to change your project name. Find and replace these strings 
(exluding quotes) throughout the source tree (within files and file names) with
a string appropriate to your project:

  "VisualPSFOptimizer"   - used in variables and paths
  "Visual PSF Optimizer" - used for presentation to humans
  "vpsfo"                        - prefix for library names
  

CHANGE ICONS

There are two icon files, one for Mac OS X and one for Windows. These are
located in 
  apps/VisualPSFOptimizer/installer/MacOSX/VisualPSFOptimizer.icns
  apps/VisualPSFOptimizer/installer/Win32/Win32Icon.ico

Modify these icons using appropriate tools. On Mac OS X, Icon Composer works
well. On Windows, you can download one of numerous icon editors.


BUILD YOUR NEW APPLICATION

Create a directory for the build and configure the new project with CMake.
Congratulations, you are now ready to specialize your new Qt/VTK/ITK-based
application!
