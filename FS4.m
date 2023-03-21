% Script to make Figure S4 - the depth split GR-MFD.
clear;

% Predefine some stuff.
dM=0.1;
Mc=0.15;
GREY=[0.85,0.85,0.85];
Dsplit=4.5;

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
tB=-1;
mB=[-4 8];
idB=-1;
ddB=1;

% Load in catalogue data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,idB,ddB);

% GR-MFD stats.
[b_D,b_err_D,a_D,R2_D,~,Mgr_D,Ngr_D,ngr_D]=Bval(Meq(dep_eq>Dsplit), Mc,dM);
[b_S,b_err_S,a_S,R2_S,~,Mgr_S,Ngr_S,ngr_S]=Bval(Meq(dep_eq<Dsplit), Mc,dM);
%[~,~,~,~,~,Mgr_ag,Ngr_ag,ngr_ag]=Bval(Mag, Mc,dM);
po_D=[-b_D,a_D]; Mgr_fit_Da=[Mc, max(Meq)]; Ngr_fit_Da=10.^polyval(po_D,Mgr_fit_Da);
po_S=[-b_S,a_S]; Mgr_fit_D=[Mc, max(Meq)]; Ngr_fit_D=10.^polyval(po_S,Mgr_fit_D);

% Plot F2.
figure(54); clf;
% The GR-MFD.
%subplot(3,3,[3 6 9]);
semilogy(Mgr_D, Ngr_D, 'o', 'Color', 'k'); hold on;
bar(Mgr_D,ngr_D, 'FaceColor', GREY);
semilogy(Mgr_fit_Da, Ngr_fit_Da, '-', 'Color', 'k');
plot(Mc*[1 1],ylim,'--k');
semilogy(Mgr_S, Ngr_S, 'o', 'Color', 'r'); hold on;
bar(Mgr_S,ngr_S, 'FaceColor', 'r');
semilogy(Mgr_fit_D, Ngr_fit_D, '-', 'Color', 'r');
xlabel('Magnitude (M_L)'); ylabel('Count');
xlim([min(Mgr_D)-dM/2 max(Mgr_D)+dM/2]); ylim([0.7 1.3*max(Ngr_D)]);