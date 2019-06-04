# celltrack
celltrack is software that finds continuous cells in 2D fields and tracks them
in time. The primary use case is rain cell tracking.

For use in scientific publications, please cite as:
Lochbihler, K., Lenderink, G., and Siebesma, A. P. ( 2017), The spatial extent of rainfall events and its relation to precipitation scaling, Geophys. Res. Lett., 44, 8629– 8636, doi:10.1002/2017GL074857.

Version 0.8 can be considered as a usable version with the following features:
 - NetCDF4 support (with compression) for input and output
 - a 2D clustering algorithm to find continuous cells
 - a linking procedure to establish connections in time between cells
 - iterative advection correction to account for displacement of cells with the mean wind field
 - support for splitting/merging tracks.
 - missing value support 
 - detection of subcells (sub-division of cells): based on an algorithm which is inspired by the watershed principle
 - output of summary statistics for cells, subcells and tracks
 - support for periodic boundary conditions
 - option to filter out cells which are smaller than a user defined threshold
 
See doc/celltrack_doc.pdf for installing/usage infos and more detailed explanations of how it works.

You will need the cdi library to compile and run celltrack: https://code.zmaw.de/projects/cdi/wiki

This software comes with absolutely no warranty! It is still under development!
