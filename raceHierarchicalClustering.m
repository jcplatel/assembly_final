function results = raceHierarchicalClustering(RACE, colorLabel, varargin)
% RACEHIERARCHICALCLUSTERING  Hierarchical clustering of synchronous
% calcium events (SCE), replacing a kmeans + successive statistical
% thresholds pipeline.
%
%   RACE       : [nSCE x nCells] binary matrix. Row i = one synchronous
%                calcium event (SCE); RACE(i,j) = 1 if cell j is active
%                during SCE i. This is your existing RACE matrix.
%   colorLabel : [nCells x 1] logical/categorical vector, e.g. true for
%                "colored"/GFP+ cells, false otherwise. Pass [] to skip
%                the color-bias step.
%
% Name-value options:
%   'NSurrogates'   (default 1000)  number of swap-randomized surrogate
%                                   RACE matrices used to build the null
%                                   distribution of pairwise distances.
%   'Quantile'      (default 0.05)  quantile of the null distribution
%                                   used as the significance cutoff.
%   'MinSCE'        (default 3)     minimum number of SCE per cluster to
%                                   keep it as a valid assembly.
%   'MinParticipation' (default 0.05) alpha for the binomial test used
%                                   to call a cell "core" of an assembly.
%   'Linkage'       (default 'average') linkage method for LINKAGE.
%   'Distance'      (default 'cosine')  distance metric for PDIST.
%
% OUTPUT results (struct):
%   .Z            linkage matrix (as returned by LINKAGE) on the SCE x SCE
%                 distance matrix -> the dendrogram clusters the SCE, i.e.
%                 the result is a clustering of the events, not directly
%                 of the cells.
%   .threshold    statistically-defined cutoff distance
%   .sceClusterID cluster id (integer) assigned to each SCE (row of RACE)
%   .assemblies   struct array, one entry per cluster with >= MinSCE SCEs:
%                     .sceIdx      indices of member SCE (rows of RACE)
%                     .coreCells   indices of cells significantly
%                                  over-represented in this cluster
%                     .pCells      binomial p-value per cell
%   .colorBias    (if colorLabel given) table of color-bias test results
%                 at every node of the dendrogram (see COLORBIASHIERARCHY)
%
% Reference: Feldt Muldoon S, Soltesz I, Cossart R. Spatially clustered
% neuronal assemblies comprise the microstructure of synchrony in
% chronically epileptic networks. PNAS 2013;110(9):3567-3572.
% Re-implemented in Python by Goirand-Lopez et al., Neurobiol Dis 2023
% (https://github.com/lgoirandlopez/Goirand-Lopez_NBD2023).

p = inputParser;
addParameter(p, 'NSurrogates', 1000);
addParameter(p, 'Quantile', 0.05);
addParameter(p, 'MinSCE', 3);
addParameter(p, 'MinParticipation', 0.05);
addParameter(p, 'Linkage', 'average');
addParameter(p, 'Distance', 'cosine');
parse(p, varargin{:});
opt = p.Results;

[nSCE, nCells] = size(RACE);

%% 1) Hierarchical clustering of the SCE (rows of RACE)
D = pdist(double(RACE), opt.Distance);
D(isnan(D)) = 1;                      % empty SCE -> maximal distance
Z = linkage(D, opt.Linkage);

%% 2) Statistical threshold from swap-randomized surrogate RACE matrices
% Null model: randomize cell identity within events while preserving
% both the number of active cells per SCE and the number of SCE each
% cell participates in (degree-preserving swap / "curveball" algorithm),
% as in Feldt Muldoon et al. 2013.
nullQuantiles = zeros(opt.NSurrogates, 1);
for s = 1:opt.NSurrogates
    Rsurr = swapRandomizeBinaryMatrix(RACE);
    Dsurr = pdist(double(Rsurr), opt.Distance);
    Dsurr(isnan(Dsurr)) = 1;
    nullQuantiles(s) = quantile(Dsurr, opt.Quantile);
end
threshold = quantile(nullQuantiles, opt.Quantile);

%% 3) Cut the dendrogram at the statistically-defined distance
sceClusterID = cluster(Z, 'cutoff', threshold, 'criterion', 'distance');

%% 4) Keep clusters with >= MinSCE events, identify core cells per cluster
uniqueClust = unique(sceClusterID);
assemblies = struct('sceIdx', {}, 'coreCells', {}, 'pCells', {});
overallRate = mean(RACE, 1);   % baseline participation rate of each cell

for c = 1:numel(uniqueClust)
    idx = find(sceClusterID == uniqueClust(c));
    if numel(idx) < opt.MinSCE
        continue
    end
    nInCluster = numel(idx);
    kActive = sum(RACE(idx, :), 1);            % per-cell count in cluster
    pCells = 1 - binocdf(kActive - 1, nInCluster, overallRate);
    coreCells = find(pCells < opt.MinParticipation);

    a = numel(assemblies) + 1;
    assemblies(a).sceIdx = idx;
    assemblies(a).coreCells = coreCells;
    assemblies(a).pCells = pCells;
end

%% 5) Color-bias test at every node of the hierarchy (optional)
colorBias = [];
if nargin >= 2 && ~isempty(colorLabel)
    colorBias = colorBiasHierarchy(Z, RACE, colorLabel);
end

results.Z = Z;
results.threshold = threshold;
results.sceClusterID = sceClusterID;
results.assemblies = assemblies;
results.colorBias = colorBias;
end


function Rsurr = swapRandomizeBinaryMatrix(R, nSwaps)
% Degree-preserving randomization of a binary matrix (rows and columns
% sums preserved) via random 2x2 swaps ("curveball"-like algorithm).
if nargin < 2
    nSwaps = 5 * nnz(R);
end
Rsurr = R;
[nRow, nCol] = size(Rsurr);
count = 0;
maxTry = nSwaps * 20;
tries = 0;
while count < nSwaps && tries < maxTry
    tries = tries + 1;
    r = randperm(nRow, 2);
    c = randperm(nCol, 2);
    block = Rsurr(r, c);
    % swap only if it changes marginals-preserving pattern (1 0 / 0 1)
    if isequal(block, [1 0; 0 1])
        Rsurr(r, c) = [0 1; 1 0];
        count = count + 1;
    elseif isequal(block, [0 1; 1 0])
        Rsurr(r, c) = [1 0; 0 1];
        count = count + 1;
    end
end
end


function T = colorBiasHierarchy(Z, RACE, colorLabel)
% For every internal node of the dendrogram (= every possible cut/merge
% of the SCE hierarchy), test whether the set of cells recruited by that
% SCE cluster is biased for colorLabel, using a Fisher exact test on the
% 2x2 table [colored cells recruited, non-colored cells recruited;
% colored cells not recruited, non-colored cells not recruited].
% p-values are corrected for multiple comparisons across all nodes
% (Benjamini-Hochberg FDR).

nSCE = size(RACE, 1);
nCells = size(RACE, 2);
colorLabel = logical(colorLabel(:));
nColorTotal = sum(colorLabel);
nNonColorTotal = nCells - nColorTotal;

sceMembers = cell(nSCE, 1);
for i = 1:nSCE
    sceMembers{i} = i;
end

nNodes = size(Z, 1);
nodeID = zeros(nNodes, 1);
mergeDist = zeros(nNodes, 1);
nSCEinCluster = zeros(nNodes, 1);
pVal = ones(nNodes, 1);

for step = 1:nNodes
    a = Z(step, 1);
    b = Z(step, 2);
    members = [sceMembers{a}, sceMembers{b}];
    node = nSCE + step;
    sceMembers{node} = members;

    % cells recruited at least once by this SCE cluster
    recruited = any(RACE(members, :) > 0, 1);
    recruitedColor = sum(recruited(colorLabel));
    recruitedNonColor = sum(recruited(~colorLabel));
    notRecruitedColor = nColorTotal - recruitedColor;
    notRecruitedNonColor = nNonColorTotal - recruitedNonColor;

    tbl = [recruitedColor, recruitedNonColor; ...
           notRecruitedColor, notRecruitedNonColor];

    pVal(step) = fisherExact2x2(tbl);
    nodeID(step) = node;
    mergeDist(step) = Z(step, 3);
    nSCEinCluster(step) = numel(members);
end

pAdj = benjaminiHochberg(pVal);

T = table(nodeID, mergeDist, nSCEinCluster, pVal, pAdj, ...
    'VariableNames', {'NodeID', 'MergeDistance', 'nSCE', 'pValue', 'pAdj'});
T = sortrows(T, 'pAdj');
end


function p = fisherExact2x2(tbl)
% Two-sided Fisher exact test for a 2x2 contingency table.
% Uses MATLAB's fishertest if the Statistics and Machine Learning Toolbox
% is available; otherwise falls back to a direct hypergeometric
% computation (equivalent result).
if exist('fishertest', 'file') == 2
    [~, p] = fishertest(tbl);
else
    a = tbl(1,1); b = tbl(1,2); c = tbl(2,1); d = tbl(2,2);
    n = a+b+c+d; K = a+c; Nrow = a+b;
    x = max(0, Nrow-(n-K)):min(Nrow, K);
    probs = hygepdf(x, n, K, Nrow);
    obsProb = hygepdf(a, n, K, Nrow);
    p = sum(probs(probs <= obsProb * (1 + 1e-10)));
    p = min(1, p);
end
end


function pAdj = benjaminiHochberg(p)
% Benjamini-Hochberg FDR correction.
n = numel(p);
[psorted, idx] = sort(p);
padj_sorted = psorted .* n ./ (1:n)';
padj_sorted = min(1, cummin(padj_sorted(end:-1:1)));
padj_sorted = padj_sorted(end:-1:1);
pAdj = zeros(n, 1);
pAdj(idx) = padj_sorted;
end
