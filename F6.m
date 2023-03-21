% Script to make Figure 6 - the fault 'wake-up' plot.
clear;

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
%latB=[  54.482   54.502   54.502   54.482]; % Big fault.
%lonB=[-118.508 -118.508 -118.501 -118.501]; % Big fault.
lonDB=[-118.58 -118.45]-0.005; latDB=[54.48 54.56]-0.005;
depB=-1;
%y=2022.00; tB=[datenum(y,01,01,00,00,00) datenum(y+0.25,01,01,00,00,00)];
tB=-1;
mB=[-4 8];
idB=-1;
ddB=1;

% Load in catalogue, injection, and station data.
file='TableS1.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('INJ.mat','S');
load('STN.mat','D');

% Flip the time ordering to show newer events on top.
I=fliplr(1:length(Teq));
lat_eq=lat_eq(I); lon_eq=lon_eq(I); dep_eq=dep_eq(I); Teq=Teq(I); Meq=Meq(I); IDeq=IDeq(I);

% Make time axis.
Ts=datenum(2016,01,01,00,00,00); % Injection start date (1 July 2017).
Teq=2016+(Teq-Ts)/365;

% Define the centre point.
i=1:4;
Clat=mean([S(i).Slat]);
Clon=mean([S(i).Slon]);
Ctime=datenum(2016,1,1,0,0,0);
Tsp=2016+(min(S(1).T)-Ctime)/365;

% Get the distances and times.
Tx=2016+(Teq-Ctime)/365;
Re=Geoid_Distance(lat_eq,lon_eq,Clat,Clon,'elliptical')*6371*pi()/180;

% Compute the diffusion curves.
Df=0.04;
tc=Tsp:1/365:max(Teq);
rc=sqrt(4*pi()*Df*(tc-min(tc))*3600*24*365)/1000;

% Make the magnitude size axis.
Rs=getMscale(Meq);

% Plot.
figure(6); clf;
% The spatiotemporal events.
subplot(121);
plot(-118.617274,54.543391,'dc'); hold on;
plot(Clon,Clat,'xk');
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b');
    plot(S(i).Wlon,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','g');
end
scatter(lon_eq, lat_eq, Rs, Teq,'filled');
xlabel('Longitude'); ylabel('Latitude');
colormap(jet); h=colorbar(); ylabel(h, 'Year');
ylim([min(latDB) max(latDB)]);
xlim([min(lonDB) max(lonDB)]);
% The distance-time plot.
subplot(122);
scatter(Teq,Re,Rs,Teq,'filled'); hold on;
plot(tc,rc,':k');
plot(Tsp*[1 1], ylim,'--k');
xlabel('Year'); ylabel('Distance (km)');
xlim([2016 2023]);



% scatter3(lon_eq,lat_eq,dep_eq,Rs,Teq, 'filled'); hold on;
% for i=1:length(S)
%     plot3(S(i).Slon,S(i).Slat,0,'ok','MarkerFaceColor','b'); hold on;
%     plot3(S(i).Wlon,S(i).Wlat,S(i).Wtvd,'-k');
% end
% ylabel('Latitude'); xlabel('Longitude'); zlabel('Depth (km)');
% colormap(jet); h=colorbar(); ylabel(h, 'Year');
% set(gca,'Zdir','reverse');
% ylim([min(latDB) max(latDB)]);
% xlim([min(lonDB) max(lonDB)]);




%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end
