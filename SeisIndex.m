function [SI]=SeisIndex(OP,OPtype,Tv,dT,SAF,Meq,Teq,b,Mc)
  % A function to perform an assessment of the seismogenic index.
  
  % Predefine the output structure.
  SI=struct('Veq',[],'SIGMA',[],'sigma',[],'sigma_err',[],'R2',[],'Vfit',[],'Nfit',[],'Vs',[],'I',[]);
  
  % Time shift the volume data.
  Tv=Tv-dT;
  
  % loop over all of the associated wells.
  for i=1:length(SAF)
      
      % Get the input values for computing the seismogenic index.
      Vc=cumsum(getVdata(OP(i),Tv,OPtype));
      T=Teq(SAF(i).Ieq);
      M=Meq(SAF(i).Ieq);
      Vstart=interp1(Tv, Vc, min(Tv), 'linear');
      Vend=0;
      sigma0=log10(0.1)-log10(max(Vc))+b*Mc;
      sigma0(isinf(sigma0))=NaN;
      sigma0=NaN;
      
      % Skip cases with no earthquakes associated to the well.
      if(isempty(T))
          SI(i).Veq=[]; SI(i).SIGMA=[]; SI(i).sigma=sigma0; SI(i).sigma_err=[]; SI(i).R2=[]; SI(i).Vfit=[]; SI(i).Nfit=[]; SI(i).Vs=[];
          SI(i).I=false;
          continue;
      end
      
      % Get the seismogenic index results.
      [Veq,SIGMA,sigma,sigma_err,R2,Vfit,Nfit,Vs]=SeismogenicIndex(Tv,Vc,T,M,b,Mc,Vstart,Vend);
      
      % Stuff results into the output structure.
      SI(i).Veq=Veq;
      SI(i).SIGMA=SIGMA;
      SI(i).sigma=sigma;
      SI(i).sigma_err=sigma_err;
      SI(i).R2=R2;
      SI(i).Vfit=Vfit;
      SI(i).Nfit=Nfit;
      SI(i).Vs=Vs;
      SI(i).I=true;
      
      % Flag if no results came back.
      if(isempty(sigma))
          SI(i).sigma=sigma0;
          SI(i).I=false;
      end
      
  end
  
return