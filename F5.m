% script to make Figure 5 - the params vs volume figures.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.2;
GREY=[0.85,0.85,0.85];
OJ=[234 182 54]/255;
v1=2936020;
v2=3495180;
v3=3767050;
v4=4622300;
v5=6113950;

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
depB=-1;
tB=-1;
mB=-1;
idB=-1;
ddB=-1;

% Load in catalogue & injection data.
file='TableS1.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('INJ.mat','S');

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(Meq, Mc,dM);

% Make time axis.
Ts=datenum(2016,01,01,00,00,00); % Injection start date (1 July 2017).
Teq=2016+(Teq-Ts)/365;
for i=1:length(S)
    S(i).T=2016+(S(i).T-Ts)/365;
end
Tv=2016:1/365:2023; Tv(Tv>max(S(1).T))=[];
Vc=cumsum(getVdata(S,Tv));

% Get the SI & Mo results.
Vstart=0;
Vend=0;
[Vs,SIGMA,sigma,sigma_err,R2s,Vfit_s,Nfit,~]=SeismogenicIndex(Tv,Vc,Teq,Meq,b,Mc,Vstart,Vend);
[Vm,Mom,R2m,ShearMc,ShearMerrc,Vfit_m,Mfitm,~]=IS_moment_fit(Tv,Vc,Teq,Meq,Mc,Vstart,Vend);
[Vx,Mox,Mx,R2x,ShearMx,ShearMerrx,Vfit_x,Mfitx,~]=IS_maxM_fit(Tv,Vc,Teq,Meq,Mc,Vstart,Vend);
[Vg,Mog,Mg,R2g,gamma,  gamma_err, Vfit_g,Mfitg,~,Tv2]=IS_Galis_fit(Tv,Vc,Teq,Meq,Mc,Vstart,Vend);
Eh=cumsum(compute_hydraulic_energy(S,Tv));
[Mbr1,Mbr2,Mul]=NRBE(Meq,Mc,b);

%%% Plot.
figure(5); clf;
% Plot SI stuff (linear scale).
axL1=subplot(321);
plot(Vs,SIGMA,'-o'); hold on;
plot(v1*[1 1],ylim(),'--k');
plot(v2*[1 1],ylim(),'--k');
ylabel('Seismogenic Index'); xlabel('Cumulative Injected Volume (m^3)');
xlim([0 max(Vs)]);
grid on;
% Plot SI stuff (log scale).
axR1=subplot(322);
semilogx(Vs,SIGMA,'-o'); hold on;
semilogx(v1*[1 1],ylim(),'--k');
semilogx(v2*[1 1],ylim(),'--k');
ylabel('Seismogenic Index'); xlabel('Cumulative Injected Volume (m^3)');
xlim([min(Vs) max(Vs)]);
grid on;
% Plot Mo stuff (linear scale).
axL2=subplot(323);
plot(Vm,Mom,'-o'); hold on;
plot(Vfit_m,Mfitm,'-k');
plot(Vfit_g,Mfitg,':k');
plot(v1*[1 1],ylim(),'--k');
plot(v2*[1 1],ylim(),'--k');
%plot(v3*[1 1],ylim(),'--k');
%plot(v4*[1 1],ylim(),'--k');
%plot(v5*[1 1],ylim(),'--k');
ylabel('Cumulative Seismic Moment Release (Nm)'); xlabel('Cumulative Injected Volume (m^3)');
xlim([0 max(Vs)]);
grid on;
% Plot Mo stuff (log scale).
axR2=subplot(324);
loglog(Vm,Mom,'-o'); hold on;
loglog(Vfit_m,Mfitm,'-k');
loglog(Vfit_g,Mfitg,':k');
loglog(Vfit_m,McGarr(Vfit_m,30,'cumulative'),'--k');
loglog(v1*[1 1],ylim(),'--k');
loglog(v2*[1 1],ylim(),'--k');
%loglog(v3*[1 1],ylim(),'--k');
%loglog(v4*[1 1],ylim(),'--k');
%loglog(v5*[1 1],ylim(),'--k');
ylabel('Cumulative Seismic Moment Release (Nm)'); xlabel('Cumulative Injected Volume (m^3)');
xlim([min(Vs) max(Vs)]);
grid on;
% Plot SE stuff (linear scale).
axL3=subplot(325);
plot(Vm,(0.46/(20*1000))*Mom./interp1(Vc+(1:length(Vc))*min(Vc(Vc>0))/1e6,Eh,Vm,'linear',0),'-o'); hold on;
plot(v1*[1 1],ylim(),'--k');
plot(v2*[1 1],ylim(),'--k');
ylabel('Seismic Efficiency (-)');
xlabel('Cumulative Injected Volume (m^3)');
xlim([0 max(Vs)]);
grid on;
% Plot SE stuff (log scale).
axR3=subplot(326);
loglog(Vm,(0.46/(20*1000))*Mom./interp1(Vc+(1:length(Vc))*min(Vc(Vc>0))/1e6,Eh,Vm,'linear',0),'-o'); hold on;
loglog(v1*[1 1],ylim(),'--k');
loglog(v2*[1 1],ylim(),'--k');
ylabel('Seismic Efficiency (-)');
xlabel('Cumulative Injected Volume (m^3)');
xlim([min(Vs) max(Vs)]);
grid on;
% Link up all of the axes.
linkaxes([axL1,axL2],'x');
linkaxes([axL1,axL3],'x');
linkaxes([axR1,axR2],'x');
linkaxes([axR1,axR3],'x');


% Print out some values that I'll need.
v1
Vm(find(Vs<=v1,1,'last')+1)
datestr(((Tv2(find(Vs<=v1,1,'last')+1)-2016)*365)+Ts)
v2
Vm(find(Vs<=v2,1,'last')+1)
datestr(((Tv2(find(Vs<=v2,1,'last')+1)-2016)*365)+Ts)
%v3
%Vm(find(Vs<=v3,1,'last')+1)
%datestr(((Tv2(find(Vs<=v3,1,'last')+1)-2016)*365)+Ts)
%v4
%Vm(find(Vs<=v4,1,'last')+1)
%datestr(((Tv2(find(Vs<=v4,1,'last')+1)-2016)*365)+Ts)
v5
Vm(find(Vs<=v5,1,'last')+1)
datestr(((Tv2(find(Vs<=v5,1,'last')+1)-2016)*365)+Ts)

