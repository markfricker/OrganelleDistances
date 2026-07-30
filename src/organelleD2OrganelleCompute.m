function [morphologyStats, d2OrganelleStats, plotLines, plotPoints] = organelleD2OrganelleCompute( ...
    morphologyStats, imSize, Cidx, calibration, span, radius, neighborRadius)
%ORGANELLED2ORGANELLECOMPUTE  Per-object inter-organelle (mito-to-mito)
% nearest-neighbour distance and a radial "crowding" vector around each
% object's own perimeter.
%
%   [morphologyStats, d2OrganelleStats, plotLines, plotPoints] = ...
%       organelleD2OrganelleCompute(morphologyStats, imSize, Cidx, calibration, span, radius, neighborRadius)
%
% Replaces the earlier all-pairs organelleDistanceCompute (a complete
% graph of every object pair — O(n^2) edges, unusable to render or plot
% for realistic object counts). Instead, for every object this casts
% outward rays from a spline-smoothed resampling of its own perimeter
% (shared normalsToNearestIntersectionSpline core, also used by
% organelleD2ErCompute) against a target polyline built from every OTHER
% object's perimeter in the same plane (otherObjectsPolyline), and keeps
% only:
%   - the nearest-neighbour distance + which object it is,
%   - a full radial distance-to-nearest-other-object vector (so a mito
%     crowded on only one flank is visible, not just a single scalar),
%   - a simple local clustering readout (neighbour count / % of
%     perimeter within neighborRadius).
%
% All per-object results are written back into morphologyStats (same
% pattern as organelleD2ErCompute) so they become ordinary morphology
% columns — plottable, exportable, and colour-fillable exactly like any
% other metric, instead of living in a separate all-pairs graph object.
%
% A ray only finds a neighbour that happens to lie along a sample point's
% own local outward normal -- a neighbour sitting immediately adjacent but
% off-axis from that direction can be missed even though it's touching or
% nearly touching. A very small (1-2px), direction-agnostic rescue against
% a per-object label crop catches exactly that case (see rescueRadiusPix
% in planeD2Organelle); it only ever improves on the ray-cast's own answer
% within that tiny tolerance, so every other case is unaffected.
%
% Whichever candidate wins (ray-cast or rescue) is also checked against
% every OTHER organelle in the plane (segmentsCrossOtherOrganelles): a
% straight line that cuts through a third organelle's own bulk on the way
% to the chosen neighbour isn't a legitimate "unobstructed access" answer.
% Obstructed points report NaN, same as "nothing found".
%
% INPUTS
%   morphologyStats – {nC_in x nZ x nT} cell array of organelle stats
%                     tables. Each table must have
%                     .organellePerimeterIdxList and .organelleCentroid.
%   imSize          – [nY nX]; size of the 2-D image plane the perimeter
%                     linear indices were traced on.
%   Cidx            – [1 x nCh] channel indices into morphologyStats.
%   calibration     – scalar µm per pixel.
%   span            – tangent/normal smoothing window (px), same meaning
%                     as in organelleD2ErCompute -- decoupled from ray
%                     density, which is always full (one per raw
%                     perimeter pixel).
%   radius          – outward ray search length (px). Also used, with a
%                     margin, as a coarse centroid-distance prefilter so
%                     only plausibly-reachable neighbours are considered
%                     when building each object's target polyline (keeps
%                     cost roughly linear in object count instead of
%                     quadratic for large fields of objects).
%   neighborRadius  – µm; radius used for the neighbour-count / percent-
%                     perimeter-near clustering readout.
%
% OUTPUTS
%   morphologyStats     – input with new columns written back per object:
%       organelleNnDistancePix    – nearest-neighbour distance (px)
%       organelleNnDistance       – nearest-neighbour distance (µm)
%       organelleNnID             – organelleID of the nearest neighbour
%       organelleNnRadialIdxList  – [x y distance] per perimeter sample,
%                                    distance to the nearest OTHER object
%       organelleNeighborCount    – # distinct other objects with any
%                                    perimeter sample within neighborRadius
%       organelleNnPctNear        – % of perimeter samples within
%                                    neighborRadius of another object
%   d2OrganelleStats    – {nCh x nZ x nT} cell array; same tables as
%                         written into morphologyStats, for callers that
%                         want the per-step result without re-deriving it.
%   plotLines           – {nCh x nZ x nT} cell array of [* x 2] overlay
%                         lines (source point -> nearest-neighbour point).
%   plotPoints          – {nCh x nZ x nT} cell array of [* x 2] contact points.

[nC_in, nZ, nT] = size(morphologyStats);
nCh = numel(Cidx);

d2OrganelleStats = cell(nCh, nZ, nT);
plotLines        = cell(nCh, nZ, nT);
plotPoints       = cell(nCh, nZ, nT);

neighborRadiusPix = neighborRadius / calibration;
reachPix          = radius + 2 * neighborRadiusPix;   % generous coarse prefilter bound

for iT = 1:nT
    for iZ = 1:nZ
        for iCh = 1:nCh
            iC = min(Cidx(iCh), nC_in);

            statsIn = morphologyStats{iC, iZ, iT};
            if isempty(statsIn)
                continue
            end

            [stats, pLines, pPoints] = planeD2Organelle( ...
                statsIn, imSize, calibration, span, radius, neighborRadiusPix, reachPix);

            d2OrganelleStats{iCh, iZ, iT} = stats;
            morphologyStats{iC, iZ, iT}   = stats;
            plotLines{iCh, iZ, iT}        = pLines;
            plotPoints{iCh, iZ, iT}       = pPoints;
        end
    end
end

end % organelleD2OrganelleCompute


% =========================================================================
function [stats, plotLines, plotPoints] = planeD2Organelle( ...
    statsIn, imSize, calibration, span, radius, neighborRadiusPix, reachPix)
%PLANED2ORGANELLE  Single-plane worker — see organelleD2OrganelleCompute for docs.

perimList = statsIn.organellePerimeterIdxList;
centroids = statsIn.organelleCentroid;   % [x y]
pixList   = statsIn.organellePixelIdxList;
nO        = height(statsIn);

stats = statsIn;
stats.organelleNnDistancePix   = nan(nO,1);
stats.organelleNnDistance      = nan(nO,1);
stats.organelleNnID            = nan(nO,1);
stats.organelleNnRadialIdxList = cell(nO,1);
stats.organelleNeighborCount   = zeros(nO,1);
stats.organelleNnPctNear       = zeros(nO,1);

plotLinesCell  = cell(nO,1);
plotPointsCell = cell(nO,1);

rayOpts          = struct();
rayOpts.maxRange = radius;

% Whole-plane label image (each organelle's own pixels marked with its own
% row index), for the same small (1-2px), direction-agnostic rescue used
% in organelleD2ErCompute: the ray-cast only finds a neighbour that
% happens to lie along a boundary point's own local outward normal, so a
% genuinely adjacent neighbour sitting off-axis from that direction can be
% missed even though it's touching or nearly touching. Cropped to a small
% local window per object (not a full-image bwdist) since the rescue
% radius is tiny -- cheap, and the crop itself excludes the querying
% object's own label, so there's no self-match risk.
rescueRadiusPix = 2;
L = zeros(imSize, 'int32');
for k = 1:nO
    L(pixList{k}) = k;
end

for iO = 1:nO
    perim = perimList{iO,1};
    if numel(perim) < 3 || nO < 2
        plotLinesCell{iO}  = zeros(0,2);
        plotPointsCell{iO} = zeros(0,2);
        continue
    end

    [targetPoly, segOwner] = otherObjectsPolyline(perimList, imSize, iO, centroids, reachPix);
    if isempty(targetPoly)
        plotLinesCell{iO}  = zeros(0,2);
        plotPointsCell{iO} = zeros(0,2);
        continue
    end

    srcContour = contourFromPerimeterIdx(perim, imSize);

    % Pinpoint exactly which object and which raw perimeter indices are
    % responsible before this reaches computeSplineNormals's generic
    % (caller-blind) non-finite check further down the call chain.
    if any(~isfinite(srcContour(:)))
        objID = NaN;
        if ismember('organelleID', statsIn.Properties.VariableNames)
            objID = statsIn.organelleID(iO);
        end
        perimD = double(perim(:));
        outOfRange = ~isfinite(perimD) | perimD < 1 | perimD > prod(imSize);
        badPerimVals = perimD(outOfRange);
        error('organelleD2OrganelleCompute:nonFiniteSrcContour', ...
            ['Object iO=%d (organelleID=%s) decoded to a non-finite contour. ' ...
            'numel(perim)=%d, min(perim)=%s, max(perim)=%s, prod(imSize)=%d, imSize=[%d %d], ' ...
            '%d perim value(s) outside [1, prod(imSize)] or non-finite (first few: %s). ' ...
            'Such a value decodes to NaN via ind2sub without erroring -- check ' ...
            'organellePerimeterIdxList for this object against the imSize ' ...
            'actually passed in.'], ...
            iO, num2str(objID), numel(perim), num2str(min(perimD)), num2str(max(perimD)), ...
            prod(imSize), imSize(1), imSize(2), ...
            numel(badPerimVals), mat2str(badPerimVals(1:min(10,end))'));
    end

    % Ray/plot density is now always full (one sample per raw perimeter
    % pixel) -- `span` no longer trades off density against tangent
    % stability, it controls only the smoothing window below.
    rayOpts.numSamples   = numel(perim);
    rayOpts.smoothSpanPx = max(span, 1);
    rayOpts.segOwner     = segOwner;

    [dist, hitPts, sampled, ~, hitOwner] = normalsToNearestIntersectionSpline(srcContour, targetPoly, rayOpts);

    % rescue: nearest OTHER labeled organelle within a tiny local crop, for
    % the case the ray-cast's normal direction missed a genuinely close
    % neighbour that's touching or nearly touching.
    sampledLin = sub2ind(imSize, ...
        min(max(round(sampled(:,1)),1),imSize(1)), min(max(round(sampled(:,2)),1),imSize(2)));
    [qr, qc] = ind2sub(imSize, sampledLin);
    rMin = max(1, min(qr) - rescueRadiusPix - 1);
    rMax = min(imSize(1), max(qr) + rescueRadiusPix + 1);
    cMin = max(1, min(qc) - rescueRadiusPix - 1);
    cMax = min(imSize(2), max(qc) + rescueRadiusPix + 1);
    Lcrop = L(rMin:rMax, cMin:cMax);
    otherMaskCrop = (Lcrop > 0) & (Lcrop ~= iO);
    if any(otherMaskCrop(:))
        [Dcrop, idxCrop] = bwdist(otherMaskCrop);
        cropLin    = sub2ind(size(Lcrop), qr - rMin + 1, qc - cMin + 1);
        rescueDist = double(Dcrop(cropLin));
        useRescue  = rescueDist <= rescueRadiusPix & (isnan(dist) | rescueDist < dist);
        if any(useRescue)
            nearestLin = idxCrop(cropLin(useRescue));
            [rRow, rCol] = ind2sub(size(Lcrop), nearestLin);
            dist(useRescue)      = rescueDist(useRescue);
            hitPts(useRescue, 1) = rRow + rMin - 1;
            hitPts(useRescue, 2) = rCol + cMin - 1;
            hitOwner(useRescue)  = double(Lcrop(nearestLin));
        end
    end

    % other-organelle obstruction: the straight line to the winning
    % neighbour may cut through some THIRD organelle's own bulk on the
    % way -- not a legitimate "unobstructed" answer. Excludes both the
    % source (self) and the target being reached (arriving at it
    % necessarily touches its own pixels).
    obstructed = segmentsCrossOtherOrganelles(sampled, hitPts, L, [repmat(iO, numel(dist), 1), hitOwner]);
    dist(obstructed)      = NaN;
    hitPts(obstructed, :) = NaN;
    hitOwner(obstructed)  = NaN;

    points = [hitPts(:,2), hitPts(:,1), dist];   % [x y distance]

    stats.organelleNnRadialIdxList{iO,1} = points;

    valid = isfinite(dist);
    if any(valid)
        [minD, iMin] = min(dist(valid));
        validIdx = find(valid);
        ownerAtMin = hitOwner(validIdx(iMin));

        stats.organelleNnDistancePix(iO,1) = minD;
        stats.organelleNnDistance(iO,1)    = minD .* calibration;
        stats.organelleNnID(iO,1)          = statsIn.organelleID(ownerAtMin);

        nearMask = dist < neighborRadiusPix;
        stats.organelleNeighborCount(iO,1) = numel(unique(hitOwner(nearMask & isfinite(hitOwner))));
        stats.organelleNnPctNear(iO,1)     = 100 * sum(nearMask) / numel(dist);
    end

    plSeg = arrayfun(@(x1,y1,x2,y2) [x1 y1; x2 y2; nan nan], ...
        sampled(:,2), sampled(:,1), points(:,1), points(:,2), 'UniformOutput', false);
    plotLinesCell{iO}  = [nan nan; cat(1, plSeg{:})];
    plotPointsCell{iO} = [nan nan; points(:,1:2)];
end

plotLines  = cat(1, plotLinesCell{:});
plotPoints = cat(1, plotPointsCell{:});

end % planeD2Organelle
