% Script to make Figure S7 - the plot of episodic slip on the (other) big fault.
clear;

% Boundary polygons and stuff.
latB=[  54.485   54.505   54.505   54.485]; % Other big fault.
lonB=[-118.496 -118.496 -118.485 -118.485]; % Other big fault.
depB=-1;
tB=-1;
mB=[-4 8];
idB=-1;
ddB=1;

% Load in catalogue, injection, and station data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('/Users/rjs10/Desktop/Musreau/data/Injection/INJ.mat','S');
load('/Users/rjs10/Desktop/Musreau/data/Stations/Station_XMLs/STN.mat','D');

% Flip the time ordering to show newer events on top.
I=fliplr(1:length(Teq));
lat_eq=lat_eq(I); lon_eq=lon_eq(I); dep_eq=dep_eq(I); Teq=Teq(I); Meq=Meq(I); IDeq=IDeq(I);

% Make time axis.
Ts=datenum(2016,01,01,00,00,00); % Injection start date (1 July 2017).
Teq=2016+(Teq-Ts)/365;

% Define the centre point.
i=1;
Clat=mean([S(i).Slat]);
Clon=mean([S(i).Slon]);
Ctime=datenum(2016,1,1,0,0,0);
Tsp=2016+(min(S(1).T)-Ctime)/365;

% Get the distances and magnitude size axis.
Re=Geoid_Distance(lat_eq,lon_eq,Clat,Clon,'elliptical')*6371*pi()/180;
Rs=getMscale(Meq)/4;

% Plot.
figure(57); clf;
% The distance-time plot for the big fault.
subplot(1,10,1:4);
scatter(Teq,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Year'); ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim([2019.0 2022.5]); ylim([0.7 2.5]);
% Episode 1.
subplot(1,10,5);
Xlim=[2019.100 2019.970];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);
% Episode 2.
subplot(1,10,6);
Xlim=[2019.970 2020.020];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);
% Episode 3.
subplot(1,10,7);
Xlim=[2020.100 2020.700];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);
% Episode 4.
subplot(1,10,8);
Xlim=[2020.800 2021.300];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);
% Episode 5.
subplot(1,10,9);
Xlim=[2021.300 2021.700];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);
% Episode 6.
subplot(1,10,10);
Xlim=[2021.700 2022.400];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.7 2.5]);


% The spatiotemporal events.
%subplot(121);
%plot(-118.617274,54.543391,'dc'); hold on;
%plot(Clon,Clat,'xk');
%for i=1:length(S)
%    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b');
%    plot(S(i).Wlon,S(i).Wlat,'-k');
%end
%for i=1:length(D)
%    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','g');
%end
%scatter(lon_eq, lat_eq, Rs, Teq,'filled');
%xlabel('Longitude'); ylabel('Latitude');
%colormap(jet); h=colorbar(); ylabel(h, 'Year');
%ylim([min(latDB) max(latDB)]);
%xlim([min(lonDB) max(lonDB)]);




%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end
