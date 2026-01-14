function Raster = transient_detection_std (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);

[NCell, Nz] = size(Tr1b);
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);
ampli = cell(1,NCell);
minithreshold=0.1; 

for i=1:NCell    
    
    th(i)=threshold_peak*std(Tr1b(i,WinRest));
    % th(i)=max ([2*iqr(Tr1b(i,WinRest)) ,2.576*std(Tr1b(i,WinRest))]) ;
    % [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i),'MinPeakDistance',MinPeakDistance);
    [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i) ,'MinPeakDistance',MinPeakDistance);%2.576= 99%, 3.291=99.9
    valeurs_identiques = intersect (locs,WinActive);
    % ampALL{i}=amplitude;
    [locs_sans_ide , idx ]=setdiff(locs(:), valeurs_identiques);
    ampli_sans_ide=amplitude(idx);
    % ampALL{i}=ampli_sans_ide;
    % Acttmp2{i}=locs_sans_ide(ampli_sans_ide>0.05 & ampli_sans_ide<1);
    Acttmp2{i}=locs_sans_ide;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
    %Acttmp2{i}=locs;
    % ampli{i}=amplitude;

end


for i = 1:NCell
    if Acttmp2{i}>0
        Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
    end
end