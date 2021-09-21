# rf-ccp
Common Conversion Point Stacking for receiver functions analysis. Input data must be in SAC format. Requires the sacio and fftw3 libraries.

Exampe parameter file:

    lat0,lon0 (initial latitude and longitude)
    lat1,lon1 (final latitude and longitude)
    t0,t1 (desired initial and final time of data)
    dx,dz (lateral and vertical sampling of profile)
    zmin,zmax (initial and final depth of profile)
    hw (half width of lateral sampling)
    outfile (output file with profile data)
    model (earth model, in format Z,Vp,Vs)
    p (name of ray parameter variable in SAC header)
    v (exponential term for phase weighting. 0 for linear stacking, i.e. no phase weighting)
    
This program migrates RFs laterally along the first fresnel zone by default.
