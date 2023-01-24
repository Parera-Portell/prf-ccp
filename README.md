# rf-ccp
Common Conversion Point Stacking for receiver function (RF) analysis. Input data must be in SAC format.

Example parameter file:

    lat0,lon0 (initial latitude and longitude)
    lat1,lon1 (final latitude and longitude)
    t0 (initial time of data)
    dx,dz (lateral and vertical resolution of profile, or cell size, in km)
    zmin,zmax (initial and final depth of profile, in km)
    hw (half width of lateral sampling, in km)
    outfile (output file path)
    model (path to earth model text file, in format Z,Vp,Vs)
    p (name of ray parameter variable in SAC header)
    zv (exponential term for depth-amplitude scaling. As amplitude decreases with depth, it is sometimes useful to apply a scale factor. The amplitude (A) is recalculated as A=A*(z**zv), where z is depth. Leave as 0 if you do not want depth scaling at all)
    v (exponential term for phase weighting. 0 for linear stacking, i.e. no phase weighting)
    
Required libraries

    sacio
    fftw3
    
This program migrates Ps waves along the first fresnel zone and outputs two text files:
1) X (distance in km), Z (depth in km), A (amplitude), LAT, LON, Z (depth in deg)
2) X (distance in km), Z (depth in km), N (number of stacked RFs)

To run the program:

    ccp [path to parameter file] [path to RF list]
    
RFs in the list must include the full path.
