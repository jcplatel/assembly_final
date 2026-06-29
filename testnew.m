lengthTreadmill=192; %Length of Bruker 's treadmill in cm
PosT=position/lengthTreadmill; %Normalise the position
PosT=PosT'; %position transpose
PosT=round(PosT,2,TieBreaker="plusinf");
PosRun=PosT(find(speed>2));%position during run
if isempty(PosRun)
    nb_lap = 0;
    return
end
tmp=diff(PosRun);
[peaks,locs]=findpeaks(-tmp,'MinPeakHeight',0.9); %find the beginning of each lap
locs(2:length(locs)+1) = locs;
locs(1) = 0 ;
nb_lap=length(locs);