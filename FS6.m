% Script to make Figure S6 - the plot of episodic slip on the big fault.
clear;

% Boundary polygons and stuff.
latB=[  54.482   54.502   54.502   54.482]; % Big fault.
lonB=[-118.508 -118.508 -118.501 -118.501]; % Big fault.
depB=-1;
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
i=1;
Clat=mean([S(i).Slat]);
Clon=mean([S(i).Slon]);
Ctime=datenum(2016,1,1,0,0,0);
Tsp=2016+(min(S(1).T)-Ctime)/365;

% Get the distances and magnitude size axis.
Re=Geoid_Distance(lat_eq,lon_eq,Clat,Clon,'elliptical')*6371*pi()/180;
Rs=getMscale(Meq)/4;

% Plot.
figure(56); clf;
% The distance-time plot for the big fault.
subplot(1,10,1:4);
scatter(Teq,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Year'); ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim([2019.5 2022.5]); ylim([0.6 2.8]);
% Episode 1.
subplot(1,10,5);
Xlim=[2019.997 2020.004];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);
% Episode 2.
subplot(1,10,6);
Xlim=[2020.15 2020.20];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);
% Episode 3.
subplot(1,10,7);
Xlim=[2020.480 2020.486];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);
% Episode 4.
subplot(1,10,8);
Xlim=[2020.855 2020.880];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);
% Episode 5.
subplot(1,10,9);
Xlim=[2021.425 2021.460];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);
% Episode 6.
subplot(1,10,10);
Xlim=[2021.995 2022.040];
scatter((Teq-min(Xlim))*365,Re,Rs,Teq,'filled','MarkerEdgeColor','k'); hold on;
xlabel('Time (days)'); %ylabel('Distance from the 05-02 Well (km)');
colormap(jet); %h=colorbar(); ylabel(h, 'Year');
xlim((Xlim-min(Xlim))*365); ylim([0.6 2.8]);

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
