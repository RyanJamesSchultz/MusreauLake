function [lat,lon,dep,T,M,ID]=parseRYN(filename,lat_L,lon_L,dep_L,T_L,M_L,ID_list,DDflag)
  % Simple routine to parse csv files for Muesreau (by Ryan, final version).
  
  % Load in the EQ data.
  % yyyy, mm, dd, hh, mm, ss, lat, lon, dep, mag
  command=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk -F[,]  ''{if($6==-9.0) $6="NaN"} 1'' | ', ... % Parse csv files, and put in NaN for null values.
           'awk -F[,TZ]  ''{print $1, $2, $3, $4, $5, $6, $7, $8, $9}'' | ' ... % Grab just the necessary fields.
           'awk  ''{split($2,d,"-"); split($3,t,":"); print $1, d[1],d[2],d[3], t[1],t[2],t[3], $4, $5, $6, $7, $8}'' > temp.EventParse']; % Split the time/date field into its respective parts.
  system(command);
  data=load('temp.EventParse');
  system('rm -f temp.EventParse');
  ID=data(:,1);
  T=datenum(data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,7));
  lat=data(:,8);
  lon=data(:,9);
  dep=data(:,10);
  M=data(:,11);
  dd=data(:,12);
  
  % Filter spatially (lateral).
  if(lat_L~=-1)
      I=inpolygon(lon,lat,lon_L,lat_L);
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end
  
  % Filter spatially (depth).
  if(dep_L~=-1)
      I=(dep>=min(dep_L))&(dep<=max(dep_L));
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end
      
  % Filter temporally.
  if(T_L~=-1)
      I=(T>=min(T_L))&(T<=max(T_L));
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end
  
  % Filter by magnitudes.
  if(M_L~=-1)
      I=(M>=min(M_L))&(M<=max(M_L));
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end

  % Filter by given IDs.
  if(ID_list~=-1)
      I=ismember(ID,ID_list);
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end

  % Filter to just the events that have been double-difference relocated.
  if(DDflag~=-1)
      I=(dd==1);
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
      ID=ID(I);
      dd=dd(I);
  end
  
  % Sort chronologically.
  [~,I]=sort(T);
  T=T(I);
  M=M(I);
  lat=lat(I);
  lon=lon(I);
  dep=dep(I);
  ID=ID(I);
  dd=dd(I);
  
return


