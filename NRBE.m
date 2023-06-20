function [Mrb1, Mrb2, Mul]=NRBE(M,Mc,b)
  % Function to determine the next record breaking event and upper-limits
  % [Mendecki, 2016; Cao et al., 2020; Verdon & Bommer, 2021]. 
  % This routine simply applies the relationships for record-breaking 
  % statistics theory [Tata, 1969; Cooke, 1979], as applied to earthquake sequences.
  % 
  % M  - Earthquake magnitudes (same units as Mc, ordered temporally).
  % Mc - Magnitude cut-off (same units as M).
  % b  - Gutenberg-Richter b-value.
  % 
  % References:
  % 
  % Cao, N.T., Eisner, L., & Jechumtálová, Z. (2020). Next record breaking magnitude for injection induced seismicity. First Break, 38(2), 53-57, doi: 10.3997/1365-2397.fb2020010.
  % Cooke, P. (1979). Statistical inference for bounds of random variables. Biometrika, 66(2), 367-374, doi: 10.1093/biomet/66.2.367.
  % Tata, M.N. (1969). On outstanding values in a sequence of random variables. Probability Theory and Related Fields, 12 (1), 9–20,doi:10.1007/BF00538520.
  % Mendecki, A.J. (2016), Mine Seismology Reference Book: Seismic Hazard, Institute of Mine Seismology, ISBN 978-0-9942943-0-2.
  % Verdon, J.P., & Bommer, J.J. (2021). Green, yellow, red, or out of the blue? An assessment of Traffic Light Schemes to mitigate the impact of hydraulic fracturing-induced seismicity. Journal of Seismology, 25(1), 301-326, doi: 10.1007/s10950-020-09966-9.
  %
  
  % Check dimensions of input vector.
  if(~isrow(M))
      M=M';
  end
  
  % Keep only earthquakes above given threshold.
  M=M(M>=Mc);
  
  % Check to see there are still earthquakes to use.
  if(isempty(M))
      fprintf('No earthquakes above given magnitude threshold.\n');
      Mul=0; Mrb1=0; Mrb2=0;
      return;
  end
  
  % Initialize.
  Mul=zeros(size(M)); Mul(1)=M(1);
  Mrb1=Mul; Mrb2=Mul;
  
  % Convert b-value into potency Beta-value (near Eqn 11 in Mendecki, 2016).
  beta=b;
  
  % Loop over every unique entry in the catalogue.
  [~,I,~]=unique(cummax(M));
  for i=2:length(I)
      
      % Call the computing subroutine for each pointwise addition to the catalogue.
      [mrb1, mrb2, mul]=nrbe_sub(M(1:I(i)),beta);
      Mul(I(i):end)=mul;
      Mrb1(I(i):end)=mrb1;
      Mrb2(I(i):end)=mrb2;
      
  end
  
return;



% Subroutine.
function [Mrb1, Mrb2, Mul]=nrbe_sub(M,beta)
  % Subrountine to compute the upper-limit and next record-breaking event 
  % forecasts, based on order statistics [Cooke, 1976; Mendecki, 2016].
  % Call this routine n-times to get the time series progression.
  %
  
  % Compute the vectors needed for the upper-limit magnitude (Cooke, 1976; near Eqn 3).
  M_maxo=sort(unique(cummax(M)),'descend');
  n=length(M_maxo);
  i=0:(n-1); i=((1-i/n).^n-(1-(i+1)/n).^n);
  Sul=sum(M_maxo(1:end).*i);
  
  % Compute the upper-limit magnitude (Eqn 38).
  if(isempty(M_maxo))
      Mul=max(M);
  else
      Mul=2*M_maxo(1)-Sul;
  end

  % Compute the vectors needed for the next record breaking event.
  dM_max=sort(diff(unique(cummax(M))),'descend');
  n=length(dM_max);
  i=0:(n-1); i=((1-i/n).^n-(1-(i+1)/n).^n);
  Srb=sum(dM_max(1:end).*i);
  
  % Compute the next record-breaking event upperbound (Eqn 42).
  if(isempty(dM_max))
      Mrb2=max(M);
  else
      Mrb2=(2*dM_max(1)-Srb)+max(M);
  end
  
  % Compute the next record-breaking event (Eqn 40).
  Mm=10^max(M); Mr=10^(Mrb2);
  if(beta==1)
      Mrb1=log(Mr/Mm)/(Mm^(-1)-Mr^(-1));
  else
      Mrb1=beta*(Mr^(1-beta)-Mm^(1-beta))/((1-beta)*(Mm^(-beta)-Mr^(-beta)));
  end
  Mrb1=log10(Mrb1);
  
return;


