function [curvesOut, summaryOut] = organelleShuffleTest(morphologyStats, targetMask, cellIDIn, Cidx, calibration, nSimulations, code, opts)
%ORGANELLESHUFFLETEST  Are organelles closer to (or further from) a target
% structure, e.g. the ER, than expected by chance? Per-cell object-shuffle
% Monte Carlo test.
%
%   [curvesOut, summaryOut] = organelleShuffleTest(morphologyStats, ...
%       targetMask, cellIDIn, Cidx, calibration, nSimulations, code, opts)
%
% Object-based analogue of DiAna's shuffle test (Gilles et al. 2017,
% Methods 115:55-64; SDI statistic from Andrey et al. 2010, PLoS Comput
% Biol 6:e1000853), adapted for organelles against a network target:
%
%   1. Descriptor: for every organelle, the minimum edge-to-target
%      distance, read from a Euclidean distance transform of the target
%      mask (0 if the organelle overlaps the target). DiAna uses
%      centre-to-nearest-centre, which is meaningless against a network
%      and poor for elongated objects.
%   2. Null model: within each cell, every organelle is re-placed at a
%      random position (shape preserved, no mutual overlap, optionally a
%      random 90-degree rotation/flip -- see shuffleObjectsInWindow) inside
%      the cell mask minus an optional exclusion mask (e.g. nucleus); the
%      target stays fixed. The descriptor is recomputed NSIMULATIONS times.
%      The same distance-transform estimator is used for observed and
%      shuffled data, so any estimator bias cancels in the comparison --
%      the spline ray-cast of organelleD2ErCompute is deliberately NOT
%      used here (too slow for hundreds of shuffles, and a different
%      estimator for observed vs null would bias the test).
%   3. Statistics: the empirical CDF of the per-object distances, the mean
%      and pointwise [alpha/2, 1-alpha/2] envelope of the shuffled CDFs,
%      and a global Monte Carlo rank test on the signed maximum CDF
%      deviation from the reference curve (the mean of all nSim+1 curves,
%      observed included, which keeps the observed and simulated curves
%      exchangeable under the null so the rank test is exact).
%
%      NOTE: the pointwise envelope is for display only. The observed
%      curve leaving the envelope somewhere is NOT a 5% test (multiple
%      comparisons across r); use sdi / pCloser / pFurther for inference.
%
% INPUTS
%   morphologyStats - {nC_in x nZ x nT} cell array of organelle stats
%                     tables. Each needs .organellePixelIdxList (linear
%                     indices into an [nY nX] plane) and .cellID.
%   targetMask      - [nY x nX x tC x tZ x tT] logical target, e.g.
%                     erSkeleton | erCisternae. Indexed per plane with the
%                     same min(iCh,tC) channel convention as
%                     organelleD2ErCompute's erSkeleton.
%   cellIDIn        - [nY x nX x cC x cZ x cT] label image of cells (same
%                     convention as analyzerRipleyL: channel min(cC, Cidx)).
%   Cidx            - [1 x nCh] channel indices into morphologyStats.
%   calibration     - scalar microns per pixel.
%   nSimulations    - number of shuffles per cell (e.g. 99 or 199; the
%                     smallest attainable one-sided p is 1/(nSim+1)).
%   code            - filename string stamped into each output row.
%   opts            - (optional) struct:
%     .excludeMask          - [nY x nX x eC x eZ x eT] logical region
%                             objects may NOT be shuffled into (e.g. the
%                             nucleus), or [] (default). Same indexing as
%                             cellIDIn.
%     .restrictTargetToCell - true (default): only target pixels inside
%                             the same cell count, so ER belonging to a
%                             neighbouring cell cannot be "nearest".
%     .dihedral             - random 90-degree rotation/flip per shuffled
%                             object (default false = translation only,
%                             as in DiAna).
%     .touchDistance        - microns; an object counts as "touching" the
%                             target if its distance <= this (default 0,
%                             i.e. overlapping).
%     .alpha                - envelope width (default 0.05 -> 2.5-97.5%).
%     .nBins                - r-grid resolution for the CDF curves
%                             (default 200).
%     .maxTries             - see shuffleObjectsInWindow (default 1000).
%     .rngSeed              - seed for reproducibility (default [] = use
%                             the current global stream). The global RNG
%                             state is restored afterwards when set.
%     .progressFcn          - optional @(fracDone, message) -> cancelled
%                             callback, called once per cell. Returning
%                             true stops early, keeping results so far.
%
% OUTPUTS
%   curvesOut  - {nCh x nZ x nT} tables, one block of nBins rows per cell:
%                filename, channel, section, frame, cellID, r (microns),
%                cdfObs, cdfNullMean, cdfNullLo, cdfNullHi.
%   summaryOut - {nCh x nZ x nT} tables, one row per cell:
%                filename, channel, section, frame, cellID, nObjects,
%                nSimulations,
%                sdi          - fraction of shuffles whose signed max CDF
%                               deviation exceeds the observed one (mid-rank
%                               for ties). <=0.025: closer to the target
%                               than chance; >=0.975: further (DiAna's
%                               convention).
%                pCloser      - one-sided Monte Carlo p, alternative
%                               "closer than chance": (1+#{Tsim>=Tobs})/(nSim+1).
%                pFurther     - one-sided p, alternative "further than
%                               chance": (1+#{Tsim<=Tobs})/(nSim+1).
%                maxDeviation - observed signed max CDF deviation (>0:
%                               closer, <0: further).
%                meanDistObs/Null, medianDistObs/Null (microns; Null =
%                               averaged over shuffles),
%                fracTouchingObs/Null - fraction of objects within
%                               touchDistance of the target,
%                fallbackFraction - fraction of shuffled placements that
%                               could not be placed and stayed put (large
%                               values mean the window is too crowded for
%                               a meaningful null),
%                windowAreaMicron2.
%                Cells with no target pixels or no objects get a row of
%                NaN statistics rather than being dropped.

if nargin < 8 || isempty(opts)
    opts = struct();
end
excludeMask   = getOpt(opts, 'excludeMask', []);
restrictTgt   = getOpt(opts, 'restrictTargetToCell', true);
touchDistance = getOpt(opts, 'touchDistance', 0);
alpha         = getOpt(opts, 'alpha', 0.05);
nBins         = getOpt(opts, 'nBins', 200);
progressFcn   = getOpt(opts, 'progressFcn', []);
shuffleOpts   = struct('dihedral', getOpt(opts, 'dihedral', false), ...
                       'maxTries', getOpt(opts, 'maxTries', 1000));

if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    prevState  = rng(opts.rngSeed);
    restoreRng = onCleanup(@() rng(prevState)); %#ok<NASGU>
end

[nC_in, nZ, nT]      = size(morphologyStats);
[nY, nX, tC, tZ, tT] = size(targetMask);
[~, ~, cC, cZ, cT]   = size(cellIDIn);
if ~isempty(excludeMask)
    [~, ~, eC, eZ, eT] = size(excludeMask);
end
nCh = numel(Cidx);

curvesOut  = cell(nCh, nZ, nT);
summaryOut = cell(nCh, nZ, nT);

% total cell count, for progress reporting only
nCellsTotal = 0;
for iT = 1:nT, for iZ = 1:nZ, for iCh = 1:nCh
    s = morphologyStats{min(Cidx(iCh), nC_in), iZ, iT};
    if ~isempty(s) && ismember('cellID', s.Properties.VariableNames)
        nCellsTotal = nCellsTotal + numel(unique(s.cellID(s.cellID > 0)));
    end
end, end, end
nCellsDone = 0;
cancelled  = false;

for iT = 1:nT
    for iZ = 1:nZ
        for iCh = 1:nCh
            iC = min(Cidx(iCh), nC_in);
            sC = Cidx(iCh);

            statsIn = morphologyStats{iC, iZ, iT};
            if isempty(statsIn) || ~ismember('cellID', statsIn.Properties.VariableNames)
                continue
            end

            cellIDPlane = cellIDIn(:,:, min(cC,sC), min(cZ,iZ), min(cT,iT));
            targetPlane = logical(targetMask(:,:, min(tC,iCh), min(tZ,iZ), min(tT,iT)));
            if isempty(excludeMask)
                exclPlane = false(nY, nX);
            else
                exclPlane = logical(excludeMask(:,:, min(eC,sC), min(eZ,iZ), min(eT,iT)));
            end
            if restrictTgt
                distPlane = [];
            else
                distPlane = double(bwdist(targetPlane)) * calibration;
            end

            groups = unique(statsIn.cellID);
            groups = groups(groups > 0);
            curveBlocks = cell(numel(groups), 1);
            sumBlocks   = cell(numel(groups), 1);
            for iG = 1:numel(groups)
                g = groups(iG);
                meta = struct('code', code, 'channel', sC, 'section', iZ, ...
                    'frame', iT, 'cellID', g);
                [curveBlocks{iG}, sumBlocks{iG}] = cellShuffleTest( ...
                    statsIn.organellePixelIdxList(statsIn.cellID == g), ...
                    cellIDPlane == g, exclPlane, targetPlane, distPlane, ...
                    calibration, nSimulations, touchDistance, alpha, nBins, ...
                    shuffleOpts, meta);

                nCellsDone = nCellsDone + 1;
                if ~isempty(progressFcn)
                    cancelled = progressFcn(nCellsDone / max(nCellsTotal,1), ...
                        sprintf('Shuffle test: cell %d of %d', nCellsDone, nCellsTotal));
                    if cancelled
                        break
                    end
                end
            end

            keep = ~cellfun(@isempty, sumBlocks);
            if any(keep)
                curvesOut{iCh, iZ, iT}  = cat(1, curveBlocks{keep});
                summaryOut{iCh, iZ, iT} = cat(1, sumBlocks{keep});
            end
            if cancelled, return, end
        end
    end
end

end % organelleShuffleTest


% =========================================================================
function [Tcurve, Tsum] = cellShuffleTest(pixLists, cellMask, exclPlane, ...
    targetPlane, distPlane, calibration, nSim, touchDistance, alpha, nBins, ...
    shuffleOpts, meta)
%CELLSHUFFLETEST  Single-cell worker -- see organelleShuffleTest for docs.

[nY, nX] = size(cellMask);
pixLists = pixLists(~cellfun(@isempty, pixLists));
nObj     = numel(pixLists);

% --- crop to the cell + its objects (block-scramble-style speedup: all
% per-shuffle work scales with the array size, not the mask area) --------
[rM, cM] = find(cellMask);
allPix   = cat(1, pixLists{:});
[rO, cO] = ind2sub([nY nX], allPix(:));
r0 = max(1,  min([rM; rO]) - 1);   r1 = min(nY, max([rM; rO]) + 1);
c0 = max(1,  min([cM; cO]) - 1);   c1 = min(nX, max([cM; cO]) + 1);
if isempty(rM)
    [r0, r1, c0, c1] = deal(1, 0, 1, 0);
end
nYc = r1 - r0 + 1;
nXc = c1 - c0 + 1;

windowMask = cellMask(r0:r1, c0:c1) & ~exclPlane(r0:r1, c0:c1);
windowArea = nnz(windowMask) * calibration^2;

if isempty(distPlane)
    tgt = targetPlane(r0:r1, c0:c1) & cellMask(r0:r1, c0:c1);
    hasTarget = any(tgt(:));
    D = double(bwdist(tgt)) * calibration;
else
    D = distPlane(r0:r1, c0:c1);
    hasTarget = any(targetPlane(:));
end

if nObj == 0 || ~hasTarget || ~any(windowMask(:))
    Tcurve = table.empty;
    Tsum   = summaryRow(meta, nObj, nSim, nan(1,11), windowArea);
    return
end

% object pixel lists -> crop-local linear indices
cropIdx = cell(nObj, 1);
for k = 1:nObj
    [r, c] = ind2sub([nY nX], pixLists{k}(:));
    cropIdx{k} = (r - r0 + 1) + (c - c0) * nYc;
end

dObs = cellfun(@(x) min(D(x)), cropIdx);

dSim = nan(nSim, nObj);
nFallback = 0;
for s = 1:nSim
    [placed, fb] = shuffleObjectsInWindow(cropIdx, windowMask, shuffleOpts);
    dSim(s, :) = cellfun(@(x) min(D(x)), placed);
    nFallback  = nFallback + nnz(fb);
end

% --- CDF curves + global rank test --------------------------------------
rMax = max([dObs(:); dSim(:)]);
if rMax <= 0
    rMax = calibration;                 % everything overlapping: degenerate grid
end
rGrid = linspace(0, rMax, nBins);

cdfObs = mean(dObs(:) <= rGrid, 1);                 % [1 x nBins]
cdfSim = zeros(nSim, nBins);
for s = 1:nSim
    cdfSim(s, :) = mean(dSim(s, :)' <= rGrid, 1);
end

allCdf = [cdfObs; cdfSim];
ref    = mean(allCdf, 1);
dev    = allCdf - ref;
[~, iMax] = max(abs(dev), [], 2);
Tall   = dev(sub2ind(size(dev), (1:nSim+1)', iMax));
Tobs   = Tall(1);
Tsim   = Tall(2:end);

sdi      = (sum(Tsim > Tobs) + 0.5 * sum(Tsim == Tobs)) / nSim;
pCloser  = (1 + sum(Tsim >= Tobs)) / (nSim + 1);
pFurther = (1 + sum(Tsim <= Tobs)) / (nSim + 1);

cdfNullMean = mean(cdfSim, 1);
cdfNullLo   = prctile(cdfSim, 100 * alpha / 2, 1);
cdfNullHi   = prctile(cdfSim, 100 * (1 - alpha / 2), 1);

vals = [sdi, pCloser, pFurther, Tobs, ...
    mean(dObs), mean(mean(dSim, 2)), median(dObs), mean(median(dSim, 2)), ...
    mean(dObs <= touchDistance), mean(mean(dSim <= touchDistance, 2)), ...
    nFallback / (nSim * nObj)];
Tsum = summaryRow(meta, nObj, nSim, vals, windowArea);

n = nBins;
Tcurve = table(repmat({meta.code}, n, 1), repmat(meta.channel, n, 1), ...
    repmat(meta.section, n, 1), repmat(meta.frame, n, 1), ...
    repmat(meta.cellID, n, 1), rGrid(:), cdfObs(:), cdfNullMean(:), ...
    cdfNullLo(:), cdfNullHi(:), ...
    'VariableNames', {'filename','channel','section','frame','cellID','r', ...
    'cdfObs','cdfNullMean','cdfNullLo','cdfNullHi'});

end % cellShuffleTest


% =========================================================================
function T = summaryRow(meta, nObj, nSim, vals, windowArea)
names = {'sdi','pCloser','pFurther','maxDeviation','meanDistObs', ...
    'meanDistNull','medianDistObs','medianDistNull','fracTouchingObs', ...
    'fracTouchingNull','fallbackFraction'};
T = table({meta.code}, meta.channel, meta.section, meta.frame, meta.cellID, ...
    nObj, nSim, 'VariableNames', {'filename','channel','section','frame', ...
    'cellID','nObjects','nSimulations'});
for k = 1:numel(names)
    T.(names{k}) = vals(k);
end
T.windowAreaMicron2 = windowArea;
end % summaryRow


% =========================================================================
function v = getOpt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end % getOpt
