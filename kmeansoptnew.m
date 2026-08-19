function [IDXs,sCl,M,S,NClini,sumd] = kmeansoptnew(E,N,type,NClini, Mcons) % on clusterize les SCE sur la base des cellules qui y participe

%E p parameters (cells) by N Events
%N number of trials per cluster number

Ne = size(E,2);

%% Covariance matrix
if strcmp(type,'var')
    if isempty (Mcons)
        M = CovarMjc(E);
    else 
        M = Mcons;
    end
end

%IDX0=zeros(N,Ne);
NCl=NClini;
%cM = parallel.pool.Constant(M);
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);

% [IDX,~,sumd]  = kmeans(M,NCl,'Options',options,"MaxIter",1000,'OnlinePhase','on','Replicates',N); % Kmeans on distance of covariance metric
[IDX,~,sumd] = kmeans(M,NCl,'Options',options,"MaxIter",1000,'Replicates',N); % Kmeans on distance of covariance metric
%%%%
% nb_rep = 1000;
% inerties_totales = zeros(nb_rep, 1);
% 
% for i = 1:nb_rep
%     [~, ~, sumd] = kmeans(M, NCl,"MaxIter",1000, 'Replicates', 1);
%     inerties_totales(i) = sum(sumd); % On somme l'inertie de tous les clusters
% end
% 
% % Tracer l'histogramme pour visualiser la stabilité
% histogram(inerties_totales,1000)


%%%%
IDX = IDX';
% S = median(silh(M,IDX));%original
%correction calcul silhouette
D_square = 1 - M;
D_square(logical(eye(size(D_square)))) = 0; % Diagonale à 0
D_square = (D_square + D_square') / 2;      % Symétrie parfaite
D_vec = squareform(D_square); 
s = silhouette([], IDX, D_vec);
S = mean (s);
% S = mean(silh(M,IDX));%original

% s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    % sCl(i) = median(s(IDX==i));%original
    sCl(i) = mean(s(IDX==i));
end

%sort RACE/silhouette of best cluster
[sCl,xCl] = sort(sCl,'descend');
IDXs = zeros(1,Ne);
for i = 1:NCl
    IDXs = IDXs + (IDX == xCl(i))*i;
end

