function sCl = kmeansoptrndnew(E,N,NCl)   %Race,10,NCl

[Np,Ne] = size(E);

%% Randomization
ERnd = zeros(Np,Ne);

for i = 1:Ne
    ERnd(:,i) = E(randperm(Np),i);    %random perm of sce for each cell
end

% %% Covariance matrix
% M = zeros(Ne,Ne);
% for i = 1:Ne
%     for j = 1:Ne
%         M(i,j) = covnorm(ERnd(:,i),ERnd(:,j),0);
%     end
% end
% M(isnan(M)) = 0;
M_shuff = CovarMjc(ERnd); 


%% k-means
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);


%IDX = kmeans(M,NCl); %modified on 2023-10-22
IDX = kmeans(M_shuff,NCl,'Options',options,"MaxIter",1000,'Replicates',N); % Kmeans on distance of covariance metric
IDX = IDX';
%S = median(silh(M,IDX));%original

D_square = 1 - M_shuff;
D_square(logical(eye(size(D_square)))) = 0; % Diagonale à 0
D_square = (D_square + D_square') / 2;      % Symétrie parfaite
D_vec = squareform(D_square); 
s = silhouette([], IDX, D_vec);

% s = silh(M_shuff,IDX);                %calculate silhouette of this best clustering with a value by SCE
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = mean(s(IDX==i));    %original     %calculate silhouette for each cluster
end
%sCl = max(sCl);         