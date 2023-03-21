% Script to make (part of) Figure 1 - the station locations.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.9;
GREY=[0.85,0.85,0.85];

% Boundary polygons and stuff.
%latB=[  55   55.00   54.00   54]; % General area.
%lonB=[-119 -117.75 -117.75 -119]; % General area.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
%latB=[  54.65   54.65   54.44    54.44]; % Main + W.
%lonB=[-118.90 -118.40 -118.40  -118.90]; % Main + W.
%latB=[  54.42   54.42   54.26    54.26]; % S cluster.
%lonB=[-118.65 -118.40 -118.40  -118.65]; % S cluster.
%latB=[  54.61   54.61   54.40    54.40]; % W cluster.
%lonB=[-118.90 -118.60 -118.60  -118.90]; % W cluster.
%latB=[  54.73   54.73   54.62    54.62]; % NW1 cluster.
%lonB=[-118.96 -118.81 -118.81  -118.96]; % NW1 clsuter.
%latB=[  54.85   54.85   54.70    54.70]; % NW2 cluster.
%lonB=[-118.84 -118.68 -118.68  -118.84]; % NW2 clsuter.
depB=-1;
%y=2022; tB=[datenum(y,1,1,0,0,0) datenum(y+1,1,1,0,0,0)];
tB=-1;
mB=[-4 8];
idB=-1;
ddB=-1;

% Load in catalogue data.
file='Table1.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);

% Load in injection & station data.
load('INJ.mat','S');
load('STN.mat','D');



% Spatial plots.
figure(1); clf;
plot(lon_eq,lat_eq,'ok','MarkerFaceColor',GREY); hold on;
plot(lon_eq(Meq>=Mc),lat_eq(Meq>=Mc),'ok','MarkerFaceColor','g');
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b');
    plot(S(i).Wlon,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','r');
end
xlabel('Longitude'); ylabel('Latitude');
xlim([-119 -118.2]); ylim([54.35 54.75]);


