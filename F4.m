% Script to make Figure 4 - the Mohr-Coloumb and FSP-dP plots.
clear;

% Get the stress/friction values.
[Sh,Sv,SH,Pp]=getS(4.084,'Shen2019'); % overpressured.
Pp=mean([Pp 1000*4.084*9.81/1e3]); % hydrostatic porepressure.
mu=getMU(Sh,Sv,SH,Pp);
azi=38;

% Get the strike-dip angles of the faults.
%[strike,dip]=parseMT('/Users/rjs10/Desktop/Musreau/data/MTs/Yu2021_MT.csv','Yu');
%[strike,dip]=parseMT('/Users/rjs10/Desktop/Musreau/data/MTs/Li2022_MT.csv','Li');
[strike,dip,IDmt]=parseMT('CompositeMT.csv','both');
%I=[12]; strike=strike(I); dip=dip(I); IDmt=IDmt(I);

% Get the overpressure values for the known faults.
[nf,sf]=computeNS(strike,dip,Sh,Sv,SH,azi,'strike-slip'); 
Pf=-(sf-mu*nf+mu*Pp)/mu;

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
%latB=[  54.482   54.502   54.502   54.482]; % Big fault.
%lonB=[-118.508 -118.508 -118.501 -118.501]; % Big fault.
%latB=[  54.4916   54.480   54.48   54.5119]; % Base of L-shape.
%lonB=[-118.5520 -118.552 -118.46 -118.4680]; % Base of L-shape.
%latB=[  54.4934   54.5122   54.57   54.56]; % N-S stem of L-shape.
%lonB=[-118.5490 -118.4680 -118.54 -118.65]; % N-S stem of L-shape.
lonDB=[-118.58 -118.45]-0.005; latDB=[54.48 54.56]-0.005;
depB=-1;
tB=-1;
mB=[-4 8];
idB=-1;
ddB=1;
Mc=0.2;

% Load in catalogue, injection, and station data.
file='TableS1.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('INJ.mat','S');
load('STN.mat','D');

% Make the magnitude size axis.
Rs=getMscale(Meq)/5;
%Rs=ones(size(Rs));

% Get the UTM coordinates for the earthquakes.
[Xeq,Yeq]=getXY(lat_eq,lon_eq);

% Get the overpressure required for slip, for each earthquake.
[strike_eq,dip_eq]=faultOrientation(Xeq,Yeq,-dep_eq*1000);
[nf_eq,sf_eq]=computeNS(strike_eq,90*ones(size(dip_eq)),Sh,Sv,SH,azi,'strike-slip'); 
Pf_eq=-(sf_eq-mu*nf_eq+mu*Pp)/mu;

% Mohr-Coloumb plot.
figure(42); clf;
% Circles.
plot([Pp Sh Sv SH],[0 0 0 0],'xk'); hold on;
[x1,y1]=makeCircle(Sh,Sv); plot(x1,y1,'-k');
[x2,y2]=makeCircle(Sv,SH); plot(x2,y2,'-k');
[x3,y3]=makeCircle(Sh,SH); plot(x3,y3,'-k');
% Friction line.
n=Pp:105; plot(n,mu*(n-Pp),'-k');
% Overpressure area colors.
[Po,n,s]=makeGrid([min(y3) max(y3)],[min(x3) max(x3)],Pp,mu,[x2,x1,fliplr(x3)],[y2,y1,fliplr(y3)]); contourf(n,s,Po,256,'LineColor','none');
% Stress locations of known faults.
plot(nf,sf,'ok');
% Labels.
xlabel('Normal Stress (MPa)'); ylabel('Shear Stress (MPa)');
cMap = interp1([0;0.05;0.10;0.25;0.30;0.5;1.0],[0.4 0 0; 1 0 0; 1 0.5 0; 1 0.9 0; 1 1 0; 0 0.95 0; 0 0.1 0],linspace(0,1,256));
colormap(cMap); caxis([0 max(Po(:))]);
h=colorbar(); ylabel(h, 'Overpressure for Slip (MPa)');
axis equal;

% The map-view of the clusters.
figure(4); clf;
plot(-118.617274,54.543391,'dc'); hold on;
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b');
    plot(S(i).Wlon,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','g');
end
[~,I]=sort(Pf_eq,'descend');
scatter(lon_eq(I), lat_eq(I), Rs(I), Pf_eq(I),'filled');
%scatter(Xeq(I),Yeq(I),Rs(I),'ok','MarkerFaceColor','m');
xlabel('Longitude'); ylabel('Latitude');
cMap = interp1([0;0.05;0.10;0.25;0.30;0.5;1.0],[0.4 0 0; 1 0 0; 1 0.5 0; 1 0.9 0; 1 1 0; 0 0.95 0; 0 0.1 0],linspace(0,1,256));
colormap(cMap); caxis([0 max(Po(:))]);
ylim([min(latDB) max(latDB)]);
xlim([min(lonDB) max(lonDB)]);


% Histogram of event overpressures.
figure(43); clf;
histogram(Pf_eq,round(sqrt(length(Pf_eq))));
xlim([0 max(Po(:))])
xlabel('Overpressure required for slip (MPa)'); ylabel('Counts');

I=(~isnan(Pf_eq))&(Pf_eq<10);
sum(Meq(I).*Pf_eq(I))/sum(Meq(I))





%%%% SUBROUNTINES.

% Get stress state information from prior study models.
function [Sh,Sv,SH,Pp]=getS(z,study_flag)
  if(strcmpi('Shen2019',study_flag))
      % Shen, Schmitt, & Haug (2019). Quantitative constraints to the complete state of stress from the combined borehole and focal mechanism inversions: Fox Creek, Alberta. Tectonophysics, 764, 110-123, doi: 10.1016/j.tecto.2019.04.023.
      Sh=32.1*z-41.8;
      Sv=24.5*z-00.0;
      SH=mean([14.3*z+40  14.3*z+80]);
      Pp=29.1*z-39.6;
  else
      % Shen, Schmitt, Wang, & Hauck (2021). States of in situ stress in the Duvernay East Shale Basin and Willesden Green of Alberta, Canada: Variable in situ stress states effect fault stability. Journal of Geophysical Research: Solid Earth, 126(6), e2020JB021221, doi: 10.1029/2020JB021221.
      Sh=22.2*z-12.8;
      Sv=24.5*z-00.0;
      SH=mean([14.3*z+40  14.3*z+80]);
      Pp=24.8*z-23.8;
  end
  return
end

% Get critical friction value, given the stress state.
function [mu]=getMU(Sh,Sv,SH,Pp)
  nc=mean([Sh,SH]);
  r=SH-nc;
  %c=Pp;
  a=(8*Pp*nc-4*nc^2+4*r^2-4*Pp^2);
  c=(4*r^2);
  mu=sqrt(-c)/sqrt(a);
  return
end

% Load in strike/dip info for MT faults.
function [strike,dip,ids]=parseMT(filename,file_flag)
  if(strcmpi('Yu',file_flag))
      command=['cat ', filename, ' | awk -F[,]  ''/Mag/{next} {print $13, $14}'' > temp.MTParse '];
      %command=['cat ', filename, ' | awk -F[,]  ''/Mag/{next} {print $16, $17}'' > temp.MTParse '];
  elseif(strcmpi('Li',file_flag))
      command=['cat ', filename, ' | awk -F[,]  ''/Mag/{next} {print $6, $7}'' > temp.MTParse '];
  else
      command=['cat ', filename, ' | awk -F[,]  ''/Mag/{next} {print $4, $5, $1}'' > temp.MTParse '];
  end

  system(command);
  data=load('temp.MTParse');
  system('rm -f temp.MTParse');
  strike=data(:,1);
  dip=data(:,2);
  ids=data(:,3);
end

% Get the boundaries of the Mohr-Coloumb circle.
function [x,y]=makeCircle(xl,xr)
  xc=mean([xl xr]);
  radius=abs(xc-xl);
  theta=linspace(0,pi(),200);
  x=radius*cos(theta)+xc;
  y=radius*sin(theta);
  return
end

% Get the overpressure required for slip inside the Mohr-Coloumb circle.
function [Po,n,s]=makeGrid(Sb,Nb,Pp,f,np,sp)
  n=0.99*min(Nb):0.05:1.01*max(Nb);
  s=0.99*min(Sb):0.05:1.01*max(Sb);
  [N,S]=meshgrid(n,s);
  Po=-(S-f*N+f*Pp)/f;
  I=inpolygon(N,S,np,sp);
  Po(~I)=NaN;
  Po(Po<0)=0;
  return
end

% Get the shear and normal stresses, given a fault dip/strike and the stress field.
function [Sn,T]=computeNS(str,dip,Sh,Sv,SH,azi,regime_flag)
  % https://dnicolasespinoza.github.io/node38.html

  if(strcmpi('reverse',regime_flag))
      S=diag([SH,Sh,Sv]);
      a=azi;
      b=0;
      c=0;
  elseif(strcmpi('strike-slip',regime_flag))
      S=diag([SH,Sv,Sh]);
      a=azi;
      b=0;
      c=90;
  elseif(strcmpi('normal',regime_flag))
      S=diag([Sv,SH,Sh]);
      a=azi-90;
      b=90;
      c=0;
  end
  R=[cosd(a)*cosd(b),                         sind(a)*cosd(b),                         -sind(b);
     cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c), sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c), cosd(b)*sind(c);
     cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c), sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c), cosd(b)*cosd(c)];
  Sg=R'*S*R;
  
  Sn=zeros(size(str));
  T=Sn;
  
  for i=1:length(str)

      nn=[-sind(str(i))*sind(dip(i)); cosd(str(i))*sind(dip(i)); -cosd(dip(i))];
      ns=[ cosd(str(i));              sind(str(i));               0           ];
      nd=[-sind(str(i))*cosd(dip(i)); cosd(str(i))*cosd(dip(i));  sind(dip(i))];
  
      t=Sg*nn;
      Sn(i)=dot(t,nn);

      Td=dot(t,nd);
      Ts=dot(t,ns);
    
      T(i)=sqrt(Ts^2+Td^2);
      rake=atan2d(Td,Ts);
      if(rake<0)
          rake=180+rake;
      end
  end
  
  return
end

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

% Get the x,y coords from lat,lon.
function [x,y]=getXY(lat,lon)
  z1=utmzone([mean(lat), mean(lon)]);
  [ellipsoid,estr] = utmgeoid(z1);
  utmstruct = defaultm('utm');
  utmstruct.zone = z1;
  utmstruct.geoid = ellipsoid;
  utmstruct = defaultm(utmstruct);
  [x,y]=projfwd(utmstruct,lat,lon);
end

% Get the fault orientation from x/y locations.
function [STR,DIP]=faultOrientation(X,Y,Z)
  
  STR=zeros(size(X)); DIP=STR;
  n_e=[0,0,1];
  
  for i=1:length(X)
      x=X(i); y=Y(i); z=Z(i);
      r=sqrt((x-X).^2+(y-Y).^2+(z-Z).^2);
      I=(r<=300);
      Xd=X(I); Yd=Y(I); Zd=Z(I);
      DM=[Xd, Yd, ones(size(Zd))];
      B=DM\Zd;
      %Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));
      P1 = [0, 0, B(1)*0+B(2)*0+B(3)];
      P2 = [1, 0, B(1)*1+B(2)*0+B(3)];
      P3 = [0, 1, B(1)*0+B(2)*1+B(3)];
      normal=cross(P1-P2,P1-P3); normal=normal/norm(normal);
      STR(i)=atan2d(normal(1),normal(2))+90;
      DIP(i)=acosd(dot(normal,n_e));
  end
  
end


