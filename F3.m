% Script to make Figure 3 - the locations.
clear;

% Boundary polygons and stuff.
%latB=[  55   55.00   54.00   54]; % General area.
%lonB=[-119 -117.75 -117.75 -119]; % General area.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
depB=-1;
tB=-1;
mB=[-4 8];
idB=-1;
ddB=1;
lonDB=[-119 -118.2]; latDB=[54.35 54.75];
Mc=0.2;

% Load in catalogue, injection, and station data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('/Users/rjs10/Desktop/Musreau/data/Injection/INJ.mat','S');
load('/Users/rjs10/Desktop/Musreau/data/Stations/Station_XMLs/STN.mat','D');

% Make the magnitude size axis.
Rs=getMscale(Meq)/2;
I=(Meq>=Mc);

% Plot location info.
figure(3); clf;
% The map-view of the clusters.
ax1=subplot(3,3,[1 2 4 5]);
plot(-118.617274,54.543391,'dc'); hold on;
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b');
    plot(S(i).Wlon,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','g');
end
scatter(lon_eq(~I),lat_eq(~I),Rs(~I),'ok','MarkerFaceColor','m');
scatter(lon_eq(I), lat_eq(I), Rs(I), 'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude');

% The N-S depth-profile of the clusters.
ax2=subplot(3,3,[3 6]);
for i=1:length(S)
    plot(0,S(i).Slat,'ok','MarkerFaceColor','b'); hold on;
    plot(S(i).Wtvd,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(0,D(i).lat,'<k','MarkerFaceColor','g');
end
scatter(dep_eq(~I),lat_eq(~I),Rs(~I),'ok','MarkerFaceColor','m');
scatter(dep_eq(I), lat_eq(I), Rs(I), 'ok','MarkerFaceColor','r');
xlabel('Depth (km)'); ylabel('Latitude');
xlim([0 15]);

% The E-W depth-profile of the clusters.
ax3=subplot(3,3,[7 8]);
for i=1:length(S)
    plot(S(i).Slon,0,'ok','MarkerFaceColor','b'); hold on;
    plot(S(i).Wlon,S(i).Wtvd,'-k');
end
for i=1:length(D)
    plot(D(i).lon,0,'^k','MarkerFaceColor','g');
end
scatter(lon_eq(~I),dep_eq(~I),Rs(~I),'ok','MarkerFaceColor','m');
scatter(lon_eq(I), dep_eq(I), Rs(I), 'ok','MarkerFaceColor','r');
ylabel('Depth (km)'); xlabel('Longitude');
ylim([0 15]);
set(gca, 'YDir','reverse')

% The depth histogram.
ax4=subplot(3,3,[9]);
histogram(dep_eq,round(sqrt(length(dep_eq))));
xlabel('Depth (km)'); ylabel('Counts');
xlim([0 15]);
camroll(-90);

% Link up axes.
linkaxes([ax1,ax2],'y');
linkaxes([ax1,ax3],'x');
set(ax1, 'Ylim',[min(latDB) max(latDB)]);
set(ax1, 'Xlim',[min(lonDB) max(lonDB)]);
%xlim([min(lonDB) max(lonDB)]); ylim([min(latDB) max(latDB)]);




%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end