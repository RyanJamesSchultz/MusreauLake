function [Vr]=getVdata(OP,Tv)
  % Simple routine to get injected/extracted volumes from operational data structures.
  
  % Preallocate space for the output vector.
  Vr=zeros(size(Tv));
  
  % Loop over the operational data structure.
  for i=1:length(OP)
      
      % Get injeciton information for this one well.
      vt=interp1(OP(i).T,cumsum(OP(i).V),Tv,'linear',0);
      vt=[vt(1), diff(vt)];
      vt(vt<0)=0;
      
      % Add this well into the whole dataset.
      Vr=Vr+vt;
  end
  
return