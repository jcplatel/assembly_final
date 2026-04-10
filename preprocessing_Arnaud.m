function [Tr1b,Raster,MAct,Race,RasterRace,TRace,Acttmp2] = preprocessing_Arnaud(F,opts,Speed,path,sampling_rate,namefull)

%clear
%close all
%% Load current data

%load('WinRest')
%load('Tr1b')
%load('Cells')
%load('MovT')
%load('Speed')


Tr1b=F; % to check
%% Detect small calcium transients
[NCell,Nz] = size(Tr1b);

% Savitzky-Golay filter
Tr1b = sgolayfilt(Tr1b',3,5)';

% figure
% for i = 1:NCell
%     plot(MovT,Tr1b(i,:)+i-1)
%     hold on
% end
% WinRest=find(speed<=1);
% WinActive=find(speed>1);
% Detect Calcium Transients using a sliding window
% TrRest = Tr1b(:,WinRest);
Raster = zeros(NCell,Nz);
WinSize = 40;
parfor i=1:NCell
    Acttmp = zeros(1,Nz);
    Sigtmp = zeros(1,Nz);
    Trtmp = Tr1b(i,:);
    %Remove points with high baseline
    ThBurst = median(Trtmp) + iqr(Trtmp)/2;
    for j = WinSize+1:Nz-WinSize
        if Speed(1,j)<1
            Wintmp = j-WinSize:j+WinSize;
            Mediantmp = median(Trtmp(Wintmp));
            %Not active in 10 last frames and not within burst activity
            if sum(Acttmp(j-10:j-1)) == 0 && Mediantmp < ThBurst 
                Acttmp(j) = Trtmp(j) - Mediantmp > 3*iqr(Trtmp(Wintmp));
                Sigtmp(j) = (Trtmp(j) - Mediantmp) / iqr(Trtmp(Wintmp));
            end
        end
    end
    Acttmp2{i} = find(Acttmp);
    Sigtmp2{i} = Sigtmp(Acttmp2{i});
end
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;
    % plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over two consecutive frames
MAct = zeros(1,Nz-1);
for i=1:Nz-1
    MAct(i) = sum(max(Raster(:,i:i+1),[],2));
end

% Select synchronies (RACE)
Th = 5;
[~,TRace] = findpeaks(MAct,'MinPeakHeight',Th,'MinPeakDistance',4);
NRace = length(TRace);

% Create RasterPlots
Race = zeros(NCell,NRace);
RasterRace = zeros(NCell,Nz);
for i = 1:NRace
    Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);
    RasterRace(Race(:,i)==1,TRace(i)) = 1;
end

% Display race
% for i = 1:length(TRace)
%     line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
% end
% break
% savefig('traces-onset-race');
%% Save
% save('Acttmp2.mat','Acttmp2')
% save('Race.mat','Race')
% save('TRace.mat','TRace')

%% Clustering
[NCell,NRace] = size(Race);
[IDX2,sCl,M,S] = kmeansopt(Race,100,'var');
% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);

[~,x2] = sort(IDX2);
MSort = M(x2,x2);

%Race clusters
R = cell(0);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end
%Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
%Remove cells with less than 2 spikes in a given cluster
CellCl(max(CellScore,[],2)<2) = 0;
[X1,x1] = sort(CellCl);

figure(visible="on")
subplot(1,2,1)
imagesc(MSort)
colormap jet
axis image
xlabel('RACE #')
ylabel('RACE #')

subplot(1,2,2)
imagesc(Race(x1,x2),[-1 1.2])
axis image
xlabel('RACE #')
ylabel('Cell #')

% savefig('RACE_fig1');
% %% Save Clusters
% save('Clusters.mat','IDX2')

%% Remove cluster non-statistically significant

sClrnd = zeros(1,20);
for i = 1:20
    sClrnd(i) = kmeansoptrnd(Race,10,NCl);
end
NClOK = sum(sCl>max(sClrnd));
sClOK = sCl(1:NClOK)';

% save('NClustersOK.mat','NClOK')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);