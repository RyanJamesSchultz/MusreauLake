% script to make Figure 8 - the NRBE predictions.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.9;
GREY=[0.85,0.85,0.85];
ShearMc=0.2311;
gamma=1.5981e+05;
sigma=-3.0986;
b=0.8816;

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
depB=-1;
tB=-1;
mB=-1;
idB=-1;
ddB=-1;

% Load in catalogue & injection data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);
load('/Users/rjs10/Desktop/Musreau/data/Injection/INJ.mat','S');

% Make time axis.
Ts=datenum(2016,01,01,00,00,00); % Injection start date (1 July 2017).
Teq=2016+(Teq-Ts)/365;
for i=1:length(S)
    S(i).T=2016+(S(i).T-Ts)/365;
end
Tv=2016:1/365:2023; Tv(Tv>max(S(1).T))=[];
Vc=cumsum(getVdata(S,Tv));

% Make the Mmax predictions.
[Mbr1,Mbr2,Mul]=NRBE(Meq,Mc,b);
[temp,it]=unique(Teq(Meq>=Mc));
Mmax_nrbe=interp1(temp,Mbr2(it),Tv,'previous');
Mmax_mcgr=Mmax(Vc,[ShearMc],'McGarr');
Mmax_elst=Mmax(Vc,[sigma b],'van der Elst');
Mmax_gals=Mmax(Vc,[gamma],'Galis');

% Get the real maximum magnitudes.
[temp,it]=unique(Teq);
Mmax_real=interp1(temp,cummax(Meq(it)),Tv,'previous');

% Make a mixed guess.
Mmax_mean=mean([Mmax_nrbe;Mmax_mcgr;Mmax_elst;Mmax_gals]);

% Make the Mnrbe predictions.
[temp,it]=unique(Teq(Meq>=Mc));
Mmax_nrbe_n=interp1(temp,Mbr1(it),Tv,'previous');
Mmax_mcgr_n=getNRBE(Mmax_mcgr,Mmax_real,b);
Mmax_elst_n=getNRBE(Mmax_elst,Mmax_real,b);
Mmax_gals_n=getNRBE(Mmax_gals,Mmax_real,b);

% Make a mixed guess, based on the average.
Mmax_mean_n=mean([Mmax_nrbe_n;Mmax_mcgr_n;Mmax_elst_n;Mmax_gals_n]);

%%% Plot NRBE stuff.
figure(8); clf;
plot(Teq,Meq,'ok','MarkerFaceColor','g','DisplayName','Data'); hold on;
plot(Tv,Mmax_real,'-g', 'DisplayName','Real M_{NRBE}');
plot(Tv,Mmax_nrbe_n,'-b', 'DisplayName','Naive M_{NRBE}');
plot(Tv,Mmax_mcgr_n,':k', 'DisplayName','McGarr M_{NRBE}');
plot(Tv,Mmax_elst_n,'-k', 'DisplayName','van der Elst M_{NRBE}');
plot(Tv,Mmax_gals_n,'--k', 'DisplayName','Galis M_{NRBE}');
plot(Tv,Mmax_mean_n,'-r', 'DisplayName','Mean M_{NRBE}');
xlabel('Year'); ylabel('Magnitude (M_L)');
xlim([2016 2023]); ylim([Mc-0.5 4.5]);
legend('Location','Northwest');




%%%% SUBROUNTINES.

% Get the mean value of the NRBE, based on a given Mmax.
function [Mmax_next]=getNRBE(Mmax_pred,Mmax_real,b)
  Mmax_pred=10.^Mmax_pred;
  Mmax_real=10.^Mmax_real;
  Mmax_next=b*(Mmax_pred.^(1-b)-Mmax_real.^(1-b))./((1-b)*(Mmax_real.^(-b)-Mmax_pred.^(-b)));
  Mmax_next=log10(Mmax_next);
  return
end
