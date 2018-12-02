# celltrack
celltrack is software that finds continuous cells in 2D fields and tracks them
in time. The primary use case is rain cell tracking.

Version 0.6 can be considered as a first usable version with the following features:
 - NetCDF4 support (with compression) for input and output
 - a 2D clustering algorithm to find continuous cells
 - a linking procedure to establish connections in time between cells
 - iterative advection correction to account for displacement of cells with the mean wind field
 - support for splitting/merging tracks.
 - missing value support 
 - output of summary statistics for cells and tracks
 - support for periodic boundary conditions
 - option to filter out cells which are smaller than a user defined threshold
 
See doc/celltrack_doc.pdf for installing/usage infos and more detailed explanations of how it works.

You will need the cdi library to compile and run celltrack: https://code.zmaw.de/projects/cdi/wiki

This software comes with absolutely no warranty! It is still under development!
