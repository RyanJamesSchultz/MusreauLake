% Script to make Figure S2 - the other injection rate information.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.2;
GREY=[0.85,0.85,0.85];

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
depB=-1;
tB=-1;
mB=[-4 8];
idB=-1;

% Load in catalogue & injection data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB);
load('/Users/rjs10/Desktop/Musreau/data/Injection/INJ.mat','S');

% Make time axis.
Ts=datenum(2016,01,01,00,00,00); % Injection start date (1 July 2017).
Tx=2016+(Teq-Ts)/365;
for i=1:length(S)
    S(i).T=2016+(S(i).T-Ts)/365;
end
[~,Tbins]=histcounts(Tx,round(1.5*sqrt(length(Tx))));
Tc=2016:1/365:2023; Tc(Tc>max(S(1).T))=[];
Vt=getVdata(S,Tc);
Et=compute_hydraulic_energy(S,Tc);

% Make the max M time series.
Mx=cummax(Meq);

% Temporal rate plots.
figure(52); clf;
% M-T plot.
ax1=subplot(511);
plot(Tx,Meq,'ok','MarkerFaceColor',GREY); hold on;
plot(Tx(Meq>=Mc),Meq(Meq>=Mc),'ok','MarkerFaceColor','r');
plot(Tx,Mx,'-r');
xlabel('Year'); ylabel('Magnitude (M_L)');
xlim([2016 2023]);
% EQ count plot.
ax2=subplot(512);
histogram(Tx, Tbins, 'FaceColor', GREY); hold on;
histogram(Tx(Meq>=Mc), Tbins, 'FaceColor', 'r');
xlabel('Year'); ylabel('Count');
xlim([2016 2023]);
% Volume Rates.
ax3=subplot(513);
for i=1:length(S)
    plot(S(i).T,S(i).V,'DisplayName',[S(i).UWI]); hold on;
end
plot(Tc,Vt,'-k','DisplayName','Total');
xlabel('Year'); ylabel('Injeciton Rate (m^3/day)');
xlim([2016 2023]);
legend('Location','northwest');
% Presssures.
ax4=subplot(514);
for i=1:length(S)
    plot(S(i).T,S(i).P,'DisplayName',[S(i).UWI]); hold on;
end
xlabel('Year'); ylabel('Pressure (MPa)');
xlim([2016 2023]);
legend('Location','northwest');
% Hydraulic energy plot.
ax5=subplot(515);
for i=1:length(S)
    plot(S(i).T,S(i).P.*S(i).V*1e6,'DisplayName',[S(i).UWI]); hold on;
end
plot(Tc,Et,'-k','DisplayName','Total');
xlabel('Year'); ylabel('Energy (J)');
xlim([2016 2023]);
legend('Location','northwest');
% Link 'em up.
linkaxes([ax1,ax2],'x');
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax4],'x');
linkaxes([ax1,ax5],'x');



