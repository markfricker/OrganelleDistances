function [d2ErStats, morphologyStats, plotLines, plotPoints] = organelleD2ErCompute( ...
    morphologyStats, erSkeleton, erCisternae, cellBoundary, edgelist, ...
    Cidx, calibration, span, radius)
%ORGANELLED2ERCOMPUTE  Radial nearest-distance from organelle membranes to the ER.
%
%   [d2ErStats, morphologyStats, plotLines, plotPoints] = organelleD2ErCompute( ...
%       morphologyStats, erSkeleton, erCisternae, cellBoundary, edgelist, ...
%       Cidx, calibration, span, radius)
%
% Casts outward rays from a spline-smoothed resampling of every object's
% perimeter and finds the nearest intersection with the ER polyline
% (cisternae boundaries + tubule skeleton edgelist), using the shared
% normalsToNearestIntersectionSpline ray-cast core (also used by
% organelleD2OrganelleCompute for mito-to-mito distance).
%
% INPUTS
%   morphologyStats – {nC_in x nZ x nT} cell array of organelle stats tables.
%                     Each table must have .organellePixelIdxList and
%                     .organellePerimeterIdxList.
%   erSkeleton      – [nY x nX x sC x sZ x sT] logical ER skeleton image.
%   erCisternae     – [nY x nX x cC x cZ x cT] logical, or [] if absent.
%   cellBoundary    – [nY x nX x bC x bZ x bT] logical cell mask, or [] if absent.
%   edgelist        – {eC x eZ x eT} cell array of tubule vertex lists.
%   Cidx            – [1 x nCh] channel indices into morphologyStats.
%   calibration     – scalar µm per pixel.
%   span            – perimeter sample-density control (px). Higher span
%                     -> fewer, smoother spline samples around each
%                     perimeter (numSamples = max(8, round(nPerim/span))).
%                     span<=1 keeps one sample per raw perimeter pixel.
%   radius          – outward ray search length (px) -- maps directly to
%                     the ray-cast core's opts.maxRange.
%
% OUTPUTS
%   d2ErStats       – {nCh x nZ x nT} cell array; each cell is statsIn
%                     augmented with organelleErOverlapArea/Idx,
%                     organelleRadialIdxList, organelleErDistancePix,
%                     organelleErDistance.
%   morphologyStats – input morphologyStats with the ER-distance columns
%                     written back into the processed channel slots, so the
%                     caller holds a single updated table (no GUI needed).
%   plotLines       – {nCh x nZ x nT} cell array of [* x 2] overlay lines.
%   plotPoints      – {nCh x nZ x nT} cell array of [* x 2] contact points.

[nC_in, nZ, nT] = size(morphologyStats);
nCh             = numel(Cidx);
[nY, nX, sC, sZ, sT] = size(erSkeleton);

if ~isempty(erCisternae)
    [~,~,cisC,cisZ,cisT] = size(erCisternae);
    hasCis = true;
else
    [cisC,cisZ,cisT] = deal(1,1,1);
    hasCis = false;
end

if ~isempty(cellBoundary)
    [~,~,bC,bZ,bT] = size(cellBoundary);
    hasBound = true;
else
    [bC,bZ,bT] = deal(1,1,1);
    hasBound = false;
end

[eC, eZ, eT] = size(edgelist);

d2ErStats  = cell(nCh, nZ, nT);
plotLines  = cell(nCh, nZ, nT);
plotPoints = cell(nCh, nZ, nT);

for iT = 1:nT
    for iZ = 1:nZ
        for iCh = 1:nCh
            iC = min(Cidx(iCh), nC_in);

            statsIn = morphologyStats{iC, iZ, iT};
            if isempty(statsIn)
                continue
            end

            % --- ER mask (logical lookup vector) -------------------------
            erSk = erSkeleton(:,:, min(iCh,sC), min(iZ,sZ), min(iT,sT));
            if hasCis
                erCis = erCisternae(:,:, min(iCh,cisC), min(iZ,cisZ), min(iT,cisT));
                erMaskVec = (erSk | erCis);
            else
                erCis     = [];
                erMaskVec = erSk;
            end
            erMaskVec = erMaskVec(:);

            % --- ER polyline (cisternae boundaries + skeleton edgelist) --
            if hasCis
                cis0         = bwboundaries(erCis);
                cis1         = cellfun(@(x) [double(x(:,1)) double(x(:,2)); nan nan], cis0, 'UniformOutput', false);
                cisternaeAll = cat(1, cis1{:});
            else
                cisternaeAll = zeros(0,2);
            end
            tub0       = edgelist{min(eC,iCh), min(eZ,iZ), min(eT,iT)};
            tub1       = cellfun(@(x) [double(x(:,1)) double(x(:,2)); nan nan], tub0, 'UniformOutput', false);
            tubulesAll = cat(1, tub1{:});
            ER = [cisternaeAll; nan nan; tubulesAll];
            ER = funcPolylineDropSingletons(ER);

            % --- pixels outside the cell boundary (col, row) -------------
            if hasBound
                [rB,cB] = find(~cellBoundary(:,:, min(iCh,bC), min(iZ,bZ), min(iT,bT)));
                cellBg  = [cB rB];
            else
                cellBg  = [];
            end

            % --- per-plane geometry (local subfunction) ------------------
            [stats, pLines, pPoints] = planeD2Er( ...
                statsIn, ER, erMaskVec, cellBg, nY, nX, calibration, span, radius);

            d2ErStats{iCh, iZ, iT}       = stats;
            morphologyStats{iC, iZ, iT}  = stats;   % write ER columns back
            plotLines{iCh, iZ, iT}       = pLines;
            plotPoints{iCh, iZ, iT}      = pPoints;
        end
    end
end

end % organelleD2ErCompute


% =========================================================================
function [stats, plotLines, plotPoints] = planeD2Er( ...
    statsIn, ER, erMaskVec, cellBg, nY, nX, calibration, span, radius)
%PLANED2ER  Single-plane worker — see organelleD2ErCompute for argument docs.

pixList   = statsIn.organellePixelIdxList;
perimList = statsIn.organellePerimeterIdxList;

% overlap: organelle pixels coincident with the ER
overlapIdx = cellfun(@(x) x(erMaskVec(x)), pixList, 'UniformOutput', false);

stats = statsIn;
stats.organelleErOverlapArea = cellfun(@(x) numel(x).*calibration.^2, overlapIdx);
stats.organelleErOverlapIdx  = overlapIdx;

nO = height(stats);
stats.organelleRadialIdxList = cell(nO,1);
stats.organelleErDistancePix = nan(nO,1);
stats.organelleErDistance    = nan(nO,1);

plotLinesCell  = cell(nO,1);
plotPointsCell = cell(nO,1);

rayOpts          = struct();
rayOpts.maxRange = radius;

for iO = 1:nO
    perim = perimList{iO,1};
    if numel(perim) < 3
        plotLinesCell{iO}  = zeros(0,2);
        plotPointsCell{iO} = zeros(0,2);
        continue
    end

    srcContour = contourFromPerimeterIdx(perim, [nY nX]);

    if span > 1
        rayOpts.numSamples = max(8, round(numel(perim) / span));
    else
        rayOpts.numSamples = numel(perim);
    end

    [dist, hitPts, sampled] = normalsToNearestIntersectionSpline(srcContour, ER, rayOpts);

    % points already coincident with the ER: zero distance, self as hit
    sampledLin = sub2ind([nY nX], ...
        min(max(round(sampled(:,1)),1),nY), min(max(round(sampled(:,2)),1),nX));
    onEr = erMaskVec(sampledLin);
    dist(onEr)      = 0;
    hitPts(onEr, :) = sampled(onEr, :);

    points = [hitPts(:,2), hitPts(:,1), dist];   % [x y distance], matches legacy convention

    % drop points outside the image or outside the cell boundary
    outside = points(:,1) < 1 | points(:,1) > nX | ...
              points(:,2) < 1 | points(:,2) > nY;
    points(outside,:) = NaN;
    if ~isempty(cellBg)
        outside = ismember(round(points(:,1:2)), cellBg, 'rows');
        points(outside,:) = NaN;
    end

    stats.organelleRadialIdxList{iO,1} = points;
    stats.organelleErDistancePix(iO,1) = min(points(:,3));
    stats.organelleErDistance(iO,1)    = min(points(:,3)) .* calibration;

    plSeg = arrayfun(@(x1,y1,x2,y2) [x1 y1; x2 y2; nan nan], ...
        sampled(:,2), sampled(:,1), points(:,1), points(:,2), 'UniformOutput', false);
    plotLinesCell{iO}  = [nan nan; cat(1, plSeg{:})];
    plotPointsCell{iO} = [nan nan; points(:,1:2)];
end

plotLines  = cat(1, plotLinesCell{:});
plotPoints = cat(1, plotPointsCell{:});

end % planeD2Er
