% Script to make Figure S3 - the CC clustering.
clear;

% Define some values.
fileCC='/Users/rjs10/Desktop/Musreau/data/eqcorrscan/CCclust/CCmat_20.npy';
fileID='/Users/rjs10/Desktop/Musreau/data/eqcorrscan/CCclust/IDs_20.csv';
clust_type='weighted';
f=0.0;
GREY=[0.85,0.85,0.85];
lonDB=[-119 -118.2]; latDB=[54.35 54.75];

% Load in the cross-correlation similarity matrix, saved as a NumPy array.
CCm=readNPY(fileCC);
CCm=abs(1-(CCm));
Ne=length(CCm);

% Load in the event IDs.
IDs=load(fileID);

% Boundary polygons and stuff.
latB=[  54.65   54.65   54.44    54.44]; % Main area.
lonB=[-118.62 -118.40 -118.40  -118.62]; % Main area.
depB=-1;
tB=-1;
mB=[-4 8];
ddB=1;

% Load in catalogue, injection, and station data.
file='/Users/rjs10/Desktop/Musreau/data/Catalogue/final/Musreau_Cat.csv';
%file='/Users/rjs10/Desktop/Musreau/data/Catalogue/Kakwa_Event_Details.csv';
[lat_eq,lon_eq,dep_eq,Teq,Meq,IDeq]=parseRYN(file,latB,lonB,depB,tB,mB,IDs,ddB);
load('/Users/rjs10/Desktop/Musreau/data/Injection/INJ.mat','S');
load('/Users/rjs10/Desktop/Musreau/data/Stations/Station_XMLs/STN.mat','D');

% Determine cut-off threshold.
d=1-squareform(1-CCm);
med=max([median(d),0]);
MAD=median(abs(d-med));
cut=f*MAD+med;

% Cluster the data.
Tree=linkage(1-d,clust_type);
c=cluster(Tree,'Cutoff',1-cut,'Criterion','distance');
%leafOrder = optimalleaforder(Tree,pdist(1-CCm));

% Plot clustering info.
figure(53); clf;
% The dendrogram.
ax2=subplot(1,6,[3]);
[H,T,Ic]=dendrogram(Tree,Ne,'ColorThreshold',1-cut,'Labels',cellstr(int2str((1:Ne)'-1)),'Orientation','right'); hold on;
plot((1-cut)*[1 1],ylim(),':k');
set(gca, 'YDir','reverse');
xlim([0 1]);
xlabel('Dissimilarity (1-CC)');
ylabel('Event ID');
% The similarity matrix.
ax1=subplot(1,6,[1 2]);
pcolor(abs(CCm(Ic,Ic)));
set(gca, 'YDir','reverse');
colormap((jet()));
h=colorbar();
ylabel(h, 'Cross-Correlation (CC)')
linkaxes([ax1,ax2],'y');
% The map-view of the clusters.
ax3=subplot(1,6,[4 5]);
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ok','MarkerFaceColor','b'); hold on;
    plot(S(i).Wlon,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(D(i).lon,D(i).lat,'^k','MarkerFaceColor','r');
end
for i=1:max(c)
    %ids=IDs(Ic);
    I=ismember(IDeq,IDs(c==i));
    scatter(lon_eq(I),lat_eq(I),'o','filled');
    %plot(lon_eq(Meq>=Mc),lat_eq(Meq>=Mc),'ok','MarkerFaceColor','g');
end
xlabel('Longitude'); ylabel('Latitude');
xlim([min(lonDB) max(lonDB)]); ylim([min(latDB) max(latDB)]);
% The N-S depth-profile of the clusters.
ax4=subplot(1,6,[6]);
for i=1:length(S)
    plot(0,S(i).Slat,'ok','MarkerFaceColor','b'); hold on;
    plot(S(i).Wtvd,S(i).Wlat,'-k');
end
for i=1:length(D)
    plot(0,D(i).lat,'<k','MarkerFaceColor','r');
end
for i=1:max(c)
    %ids=IDs(Ic);
    I=ismember(IDeq,IDs(c==i));
    scatter(dep_eq(I),lat_eq(I),'o','filled');
    %plot(lon_eq(Meq>=Mc),lat_eq(Meq>=Mc),'ok','MarkerFaceColor','g');
end
xlabel('Depth (km)'); ylabel('Latitude');
%xlim([-119 -118.2]); ylim([54.35 54.75]);
linkaxes([ax3,ax4],'y');

% Plot the distribution of CC values.
figure(54); clf;
histogram(d); hold on;
plot(med*[1 1],ylim(),'--k');
plot(cut*[1 1],ylim(),':k');
xlabel('Cross-Correlation (CC)'); ylabel('Counts');




