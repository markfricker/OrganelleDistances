function [d2ErStats, morphologyStats, plotLines, plotPoints] = organelleD2ErCompute( ...
    morphologyStats, erSkeleton, erCisternae, cellBoundary, edgelist, ...
    Cidx, calibration, span, radius, contactDistance)
%ORGANELLED2ERCOMPUTE  Radial nearest-distance from organelle membranes to the ER.
%
%   [d2ErStats, morphologyStats, plotLines, plotPoints] = organelleD2ErCompute( ...
%       morphologyStats, erSkeleton, erCisternae, cellBoundary, edgelist, ...
%       Cidx, calibration, span, radius, contactDistance)
%
% Casts outward rays from a spline-smoothed resampling of every object's
% perimeter and finds the nearest intersection with the ER polyline
% (cisternae boundaries + tubule skeleton edgelist), using the shared
% normalsToNearestIntersectionSpline ray-cast core (also used by
% organelleD2OrganelleCompute for mito-to-mito distance).
%
% A ray only finds ER that happens to lie along a sample point's own
% local outward normal -- ER sitting immediately adjacent but off-axis
% from that direction can be missed even though it's touching or nearly
% touching. A very small (1-2px), direction-agnostic distance-transform
% rescue catches exactly that case (see rescueRadiusPix in planeD2Er); it
% only ever improves on the ray-cast's own answer within that tiny
% tolerance, so every other case -- including the normal-to-surface
% appearance of the overlay -- is unaffected.
%
% Whichever candidate wins (ray-cast or rescue) is also checked against
% every OTHER organelle in the plane (segmentsCrossOtherOrganelles): a
% straight line that cuts through another organelle's own bulk on the way
% to the ER isn't a legitimate "unobstructed access" answer, even though
% the endpoint itself is genuinely the nearest ER pixel in that direction.
% Obstructed points report NaN, same as "nothing found".
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
%   span            – tangent/normal smoothing window (px). Rays are
%                     always cast at full density (one per raw perimeter
%                     pixel) -- this only controls how far along the arc
%                     the outward-normal direction is averaged, decoupled
%                     from ray count. Larger span -> smoother-looking
%                     (less pixel-jagged) normals.
%   radius          – outward ray search length (px) -- maps directly to
%                     the ray-cast core's opts.maxRange.
%   contactDistance – (optional) µm tolerance for ER contact length. When
%                     given, adds organelleErContactFraction (fraction of
%                     the organelle's perimeter pixels within
%                     contactDistance of the ER target, overlapping pixels
%                     count as 0) and organelleErContactLength (that
%                     fraction x organellePerimeter in µm, or x perimeter
%                     pixel count x calibration if organellePerimeter is
%                     absent). Direction-agnostic (Euclidean distance
%                     transform, like DiAna's dilation-based contact
%                     surface) rather than the normal ray-cast, so an ER
%                     tubule running alongside the organelle counts along
%                     its whole length. NOTE: the target is the ER
%                     skeleton centreline (+ cisternae), so the tolerance
%                     effectively includes the tubule half-width.
%
% OUTPUTS
%   d2ErStats       – {nCh x nZ x nT} cell array; each cell is statsIn
%                     augmented with organelleErOverlapArea/Idx,
%                     organelleRadialIdxList, organelleErDistancePix,
%                     organelleErDistance (+ organelleErContactFraction /
%                     organelleErContactLength when contactDistance given).
%                     NB organelleErOverlapArea counts organelle pixels ON
%                     the ER target, i.e. for tubules the 1-px centreline
%                     crossing the organelle, not a full-width footprint.
%   morphologyStats – input morphologyStats with the ER-distance columns
%                     written back into the processed channel slots, so the
%                     caller holds a single updated table (no GUI needed).
%   plotLines       – {nCh x nZ x nT} cell array of [* x 3] overlay lines
%                     [x y distance] (distance repeated at both ends of
%                     each ray, NaN-separated between rays).
%   plotPoints      – {nCh x nZ x nT} cell array of [* x 3] contact points
%                     [x y distance].

if nargin < 10
    contactDistance = [];
end

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
                statsIn, ER, erMaskVec, cellBg, nY, nX, calibration, span, radius, contactDistance);

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
    statsIn, ER, erMaskVec, cellBg, nY, nX, calibration, span, radius, contactDistance)
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

% Very small (1-2px), direction-agnostic rescue for the one specific case
% the ray-cast is known to miss: ER sitting immediately adjacent to a
% boundary point but off-axis from that point's local outward normal, so
% no ray happens to cross it even though it's touching or nearly
% touching. This only ever improves on the ray-cast's own answer (closer,
% and within this tiny tolerance) -- everything beyond it is left
% entirely to the unmodified ray-cast below, exactly as before.
rescueRadiusPix      = 2;
erMask2D             = reshape(erMaskVec, nY, nX);
[DRescue, idxRescue] = bwdist(erMask2D);

% Whole-plane label image (every organelle's own pixels marked with its
% own row index), so a winning hit can be checked for whether the
% straight line to it passes through some OTHER organelle's own bulk --
% not a legitimate "unobstructed" answer even though the endpoint itself
% is a valid ER pixel.
Lorg = zeros(nY, nX, 'int32');
for k = 1:nO
    Lorg(pixList{k}) = k;
end

for iO = 1:nO
    perim = perimList{iO,1};
    if numel(perim) < 3
        plotLinesCell{iO}  = zeros(0,3);
        plotPointsCell{iO} = zeros(0,3);
        continue
    end

    srcContour = contourFromPerimeterIdx(perim, [nY nX]);

    % Ray/plot density is now always full (one sample per raw perimeter
    % pixel) -- `span` no longer trades off density against tangent
    % stability, it controls only the smoothing window below.
    rayOpts.numSamples   = numel(perim);
    rayOpts.smoothSpanPx = max(span, 1);

    [dist, hitPts, sampled] = normalsToNearestIntersectionSpline(srcContour, ER, rayOpts);

    % points already coincident with the ER: zero distance, self as hit
    sampledLin = sub2ind([nY nX], ...
        min(max(round(sampled(:,1)),1),nY), min(max(round(sampled(:,2)),1),nX));
    onEr = erMaskVec(sampledLin);
    dist(onEr)      = 0;
    hitPts(onEr, :) = sampled(onEr, :);

    % rescue: ER within 1-2px that the ray-cast's normal direction missed
    rescueDist = double(DRescue(sampledLin));
    useRescue  = ~onEr & rescueDist <= rescueRadiusPix & (isnan(dist) | rescueDist < dist);
    if any(useRescue)
        [rescueRow, rescueCol] = ind2sub([nY nX], idxRescue(sampledLin(useRescue)));
        dist(useRescue)      = rescueDist(useRescue);
        hitPts(useRescue, 1) = rescueRow;
        hitPts(useRescue, 2) = rescueCol;
    end

    % other-organelle obstruction: the straight line to the winning hit
    % may cut through some OTHER organelle's own bulk on the way -- not a
    % legitimate "unobstructed" answer even though the endpoint itself is
    % a valid ER pixel.
    obstructed = ~onEr & segmentsCrossOtherOrganelles(sampled, hitPts, Lorg, iO);
    dist(obstructed)      = NaN;
    hitPts(obstructed, :) = NaN;

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

    % Distance is carried through as a third column on both plotLines and
    % plotPoints (same value at both ends of a ray) -- not rendered with
    % colour yet for the radial lines, but ready for it: see
    % funcRoiOverlay's 'd2Er'/'d2Organelle' scatter case for the points
    % side, which already does colour by this value.
    plSeg = arrayfun(@(x1,y1,x2,y2,d) [x1 y1 d; x2 y2 d; nan nan nan], ...
        sampled(:,2), sampled(:,1), points(:,1), points(:,2), points(:,3), 'UniformOutput', false);
    plotLinesCell{iO}  = [nan nan nan; cat(1, plSeg{:})];
    plotPointsCell{iO} = [nan nan nan; points];
end

% Overlap area is computed from the raw, unsmoothed pixel mask, while
% distance comes from a spline-smoothed/resampled ray-cast -- a small
% overlap can fall between the sampled perimeter points and be missed by
% the ray-cast, reporting a small positive distance despite genuinely
% touching the ER. Overlap implies zero distance by definition, so it
% overrides whatever the ray-cast found.
hasOverlap = stats.organelleErOverlapArea > 0;
stats.organelleErDistancePix(hasOverlap) = 0;
stats.organelleErDistance(hasOverlap)    = 0;

% ER contact length: perimeter pixels within contactDistance of the ER
% target, direction-agnostic (DRescue is the full-plane distance transform
% of the same ER target). NaN when there is no ER in the plane at all.
if ~isempty(contactDistance)
    hasEr = any(erMaskVec);
    contactFrac = nan(nO, 1);
    contactLen  = nan(nO, 1);
    hasPerimLen = ismember('organellePerimeter', stats.Properties.VariableNames);
    for iO = 1:nO
        perim = perimList{iO,1};
        if isempty(perim) || ~hasEr
            continue
        end
        contactFrac(iO) = mean(double(DRescue(perim)) .* calibration <= contactDistance);
        if hasPerimLen
            contactLen(iO) = contactFrac(iO) .* stats.organellePerimeter(iO);
        else
            contactLen(iO) = contactFrac(iO) .* numel(perim) .* calibration;
        end
    end
    stats.organelleErContactFraction = contactFrac;
    stats.organelleErContactLength   = contactLen;
end

plotLines  = cat(1, plotLinesCell{:});
plotPoints = cat(1, plotPointsCell{:});

end % planeD2Er
