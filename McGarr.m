function [Mo]=McGarr(V,G,type_flag)
  % Simple function that will output the McGarr [1976; 2014] upperbound, given injection volumes (m^3).
  %
  % V  - Cumulative volume injected axis (units of m^3).
  % G  - Apparent Shear Modulus (GPa)
  %
  % References:
  % 
  % McGarr, A. (1976). Seismic moments and volume changes. Journal of Geophysical Research, 81(8), 1487-1494, doi: 10.1029/JB081i008p01487.
  % McGarr, A. (2014). Maximum magnitude earthquakes induced by fluid injection. Journal of Geophysical Research: solid earth, 119(2), 1008-1019, doi: 10.1002/2013JB010597.
  %
  
  % Compute McGarr upperbound (Eqn 13).
  Mo=G*V*1e9;

  % If this is a cumulative Mo plot, double it (Eqn 7).
  if(strcmpi(type_flag,'cumulative'))
      Mo=2*Mo;
  end
  
return;