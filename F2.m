% Script to make Figure 2 - the time series plots & GR-MFD.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.2;
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
ddB=1;

% Load in catalogue data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
%file='/Users/rjs10/Desktop/Musreau/data/Catalogue/Kakwa_Event_Summary.csv';
%[lat_eq,lon_eq,dep_eq,Teq,Meq]=parseNMX(file,latB,lonB,depB,tB,mB);
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/AGScatalogue.csv';
[lat_ag,lon_ag,dep_ag,Tag,Mag]=parseAGS(file,latB,lonB,depB,tB,mB);

% Load in injection data.
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

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(Meq, Mc,dM);
[~,~,~,~,~,Mgr_ag,Ngr_ag,ngr_ag]=Bval(Mag, Mc,dM);
po=[-b,a];
Mgr_fit=[Mc, max(Meq)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Plot F2.
figure(2); clf;
% The GR-MFD.
subplot(3,3,[3 6 9]);
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
plot(Mc*[1 1],ylim,'--k');
semilogy(Mgr_ag, Ngr_ag, 'o', 'Color', 'r'); hold on;
bar(Mgr_ag,ngr_ag, 'FaceColor', 'r');
plot(1.85*[1 1],ylim,'--r');
xlabel('Magnitude (M_L)'); ylabel('Count');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);
% The M-T plot.
ax1=subplot(3,3,[1 2]);
plot(Tx,Meq,'ok','MarkerFaceColor',GREY); hold on;
plot(Tx(Meq>=Mc),Meq(Meq>=Mc),'ok','MarkerFaceColor','g');
xlabel('Year'); ylabel('Magnitude (M_L)');
xlim([2016 2023]);
% The n-T & v-T plots.
ax2=subplot(3,3,[4 5]);
histogram(Tx, Tbins, 'FaceColor', GREY); hold on;
histogram(Tx(Meq>=Mc), Tbins, 'FaceColor', 'g');
plot(S(1).T,S(1).V/5,'DisplayName',[S(1).UWI]);
plot(S(2).T,S(2).V/5,'DisplayName',[S(2).UWI]);
plot(S(3).T,S(3).V/5,'DisplayName',[S(3).UWI]);
plot(S(4).T,S(4).V/5,'DisplayName',[S(4).UWI]);
plot(Tc,Vt/5,'-k','DisplayName','Total');
set(gca, 'YScale', 'log')
xlabel('Year'); ylabel('Count');
%ylabel('Injection Rate (m^3/day)');
xlim([2016 2023]);
ylim([0.7 1.1*(max(Vt)/5)]);
% The N-T & V-T plots.
ax3=subplot(3,3,[7 8]);
plot(Tc, cumsum([histcounts(Tx,Tc),0]), '-r','DisplayName','All EQs'); hold on;
plot(Tc, cumsum([histcounts(Tx(Meq>=Mc),Tc),0]), '-g','DisplayName','EQs above Mc');
plot(S(1).T,cumsum(S(1).V)/300,'DisplayName',[S(1).UWI]);
plot(S(2).T,cumsum(S(2).V)/300,'DisplayName',[S(2).UWI]);
plot(S(3).T,cumsum(S(3).V)/300,'DisplayName',[S(3).UWI]);
plot(S(4).T,cumsum(S(4).V)/300,'DisplayName',[S(4).UWI]);
plot(Tc,cumsum(Vt)/300,'-k','DisplayName','Total');
set(gca, 'YScale', 'log')
xlabel('Year'); ylabel('Count');
%ylabel('Cumulative Injection (m^3)');
xlim([2016 2023]);
legend('Location','northwest');
linkaxes([ax1,ax2],'x');
linkaxes([ax1,ax3],'x');
