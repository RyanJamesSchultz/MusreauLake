function [E]=compute_hydraulic_energy(OP,Tv)
  % Simple routine to compute the cumulaitve hydraulic energy
  
  % Initialize.
  E=zeros(size(Tv));
  
  % Loop over the operational data structure.
  for i=1:length(OP)
      
      % Get injeciton information for this one well.
      ei=OP(i).V.*OP(i).P*1e6; % Units of J.
      et=interp1([0;OP(i).T],[0;cumsum(ei)],Tv,'linear',0);
      et=[et(1), diff(et)];
      et(et<0)=0;
      
      % Add this well into the whole dataset.
      E=E+et;
  end
  
return;