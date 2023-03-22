function [Veq, Moc, M, R2, Gamma, Gamma_err, Vfit,Mfit, Vs, Tm]=IS_Galis_fit(Tv, V, Tm, M,  Mc, Vstart,Vend)
  % Function to determine the best fit to a set of induced seismicity data 
  % to the Galis relationship [Galis et al., 2017].  This routine simply 
  % fits the data to a linear relationship between max Mo vs V^(3/2).
  %
  % This routine converts Mw magnitudes to Mo via a scaling relationship 
  % [Hanks & Kanamori, 1979].
  %
  % Tv   - Cumulative volume injected time axis (same units as Tm).
  % V    - Cumulative volume injected im m^3.
  % Tm   - Earthquake time axis.
  % M    - Earthquake magnitudes (same units as Mc).
  % b    - Seismic b-value.
  % Mc   - Magnitude cut-off (same units as M).
  %
  % References:
  % 
  % Hanks, T. C., and H. Kanamori (1979), A moment magnitude scale, Journal of Geophysical Research, 84, 2348?2350, doi:10.1029/JB084iB05p02348.
  % Galis, M., Ampuero, J. P., Mai, P. M., & Cappa, F. (2017). Induced seismicity provides insight into why earthquake ruptures stop. Science advances, 3(12), eaap7528, doi: 10.1126/sciadv.aap7528.
  %

  % Define vector of indecies to output to user discarded
  %ID=1:length(Tv);
  
  % Check dimensions and sizes of input vectors.
  if( (length(Tv)~=length(V))||(length(Tm)~=length(M))   )
      fprintf('improper input lengths.\n');
      Veq=0; Moc=0; R2=0; ShearM=0; ShearMerr=0; Vfit=0; Mfit=0;
      return;
  end
  if(~isrow(Tv))
      Tv=Tv';
  end
  if(~isrow(V))
      V=V';
  end
  if(~isrow(Tm))
      Tm=Tm';
  end
  if(~isrow(M))
      M=M';
  end
  
  % Sort earthquakes by time.
  [Tm, I]=sort(Tm);
  M=M(I);
  
  % Keep only earthquakes above given threshold.
  Tm=Tm(M>=Mc); M=M(M>=Mc);
  
  % Check to see there are still earthquakes to use.
  if(isempty(M))
      fprintf('No earthquakes above given magnitude threshold.\n');
      Veq=0; Moc=0; R2=0; Gamma=0; Gamma_err=0; Vfit=0; Mfit=0;
      return;
  end
  
  % Determine total injected volume at the times of the earthquakes.
  Veq=interp1( Tv, V, Tm, 'linear');
  
  % Keep only earthquakes happening before user defined shut-in volume.
  if(Vend~=0)
      Tm=Tm(Veq<=Vend);  M=M(Veq<=Vend);  Veq=Veq(Veq<=Vend);
  end
  
  % Keep only earthquakes happening during injection interval or user defined start volume.
  Vs=Veq(1);
  if(Vstart==0)
      Veq=Veq-Veq(1);
      Tm=Tm(Veq>0);   M=M(Veq>0);  Veq=Veq(Veq>0);
  else
      Tm=Tm(Veq>=Vstart);  M=M(Veq>=Vstart);  Veq=Veq(Veq>=Vstart);
      Veq=Veq-Vstart;
  end
  
  % Keep only eartqhaukes happening during stimulation periods.
  %[~,I,~]=unique(Veq);
  %Tm=Tm(I);  M=M(I);  Veq=Veq(I);
  
  % Keep only earthquakes that have interpolated injection volumes.
  Tm=Tm(~isnan(Veq));  M=M(~isnan(Veq));  Veq=Veq(~isnan(Veq));
  
  %check to see there are still earthquakes to use.
  if(isempty(M))
      fprintf('No earthquakes during injection interval.\n');
      Veq=0; Moc=0;
      return;
  end
  
  % Convert Mw to seismic moment [Hanks & Kanamori, 1976].
  Moc=10.^(1.5*M+9.1); % Seismic moment in N m.
  Moc=cumsum(Moc);
  
  % Perform linear fit to data.
  X=Veq.^(3/2);
  Y=Moc;
  a=(X*X')\(X*Y');
  p=[a,0];
  
  % Determine goodness-of-fit.
  Yfit=polyval(p,X);
  SSres=sum( (Y-Yfit).^2 );
  SStot=sum( (Y-mean(Y)).^2  );
  R2=1-(SSres/SStot);
  
  % Get the output fit parameters.
  Gamma=p(1); % In units of N/m^7/2.
  Gamma_err=sqrt(inv(X*X')); % In units of N/m^7/2.
  
  % Linear data fit.
  Vfit=linspace(Veq(1),Veq(end),length(Veq));
  Mfit=polyval(p,Vfit.^(3/2));
  
return;

