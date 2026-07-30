function [morphologyStats, d2BoundaryStats, plotLines, plotPoints] = organelleD2BoundaryCompute( ...
    morphologyStats, imSize, cellBoundary, Cidx, calibration, span, radius)
%ORGANELLED2BOUNDARYCOMPUTE  Radial nearest-distance from organelle
% membranes to the cell boundary.
%
%   [morphologyStats, d2BoundaryStats, plotLines, plotPoints] = ...
%       organelleD2BoundaryCompute(morphologyStats, imSize, cellBoundary, Cidx, calibration, span, radius)
%
% Casts outward rays from a spline-smoothed resampling of every object's
% perimeter and finds the nearest intersection with the cell boundary
% (traced from the cellBoundary mask), using the same
% normalsToNearestIntersectionSpline ray-cast core as
% organelleD2ErCompute/organelleD2OrganelleCompute. Complements the
% existing organelleBoundary touch flag (analyzerFeatureAnalysis) with an
% actual distance for objects that are NOT touching -- e.g. to separate
% cortical vs. perinuclear mitochondria by how far they sit from the
% plasma membrane.
%
% A ray only finds boundary that happens to lie along a sample point's own
% local outward normal -- boundary sitting immediately adjacent but
% off-axis from that direction can be missed even though it's touching or
% nearly touching. A very small (1-2px), direction-agnostic distance-
% transform rescue catches exactly that case (see rescueRadiusPix in
% planeD2Boundary); it only ever improves on the ray-cast's own answer
% within that tiny tolerance, so every other case is unaffected.
%
% Whichever candidate wins (ray-cast or rescue) is also checked against
% every OTHER organelle in the plane (segmentsCrossOtherOrganelles): a
% straight line that cuts through another organelle's own bulk on the way
% to the cell wall isn't a legitimate "unobstructed access" answer, even
% though the endpoint itself is genuinely the nearest boundary pixel in
% that direction. Obstructed points report NaN, same as "nothing found".
%
% INPUTS
%   morphologyStats – {nC_in x nZ x nT} cell array of organelle stats
%                     tables. Each table must have
%                     .organellePerimeterIdxList and .organellePixelIdxList.
%   imSize          – [nY nX]; size of the 2-D organelle image plane the
%                     PerimeterIdxList entries are linear into (matches
%                     organelleD2OrganelleCompute's imSize argument --
%                     NOT necessarily the same size as cellBoundary; see
%                     the size-mismatch check below).
%   cellBoundary    – [nYb x nXb x bC x bZ x bT] logical; TRUE = inside the
%                     cell (same convention as organelleBoundaryFlag /
%                     app.images.cellBoundary). May be static (bT==1).
%   Cidx            – [1 x nCh] channel indices into morphologyStats.
%   calibration     – scalar µm per pixel.
%   span            – tangent/normal smoothing window (px), same meaning
%                     as in organelleD2ErCompute -- decoupled from ray
%                     density, which is always full (one per raw
%                     perimeter pixel).
%   radius          – outward ray search length (px). Objects farther
%                     from the boundary than this get NaN -- raise it for
%                     large cells / deeply interior mitochondria.
%
% OUTPUTS
%   morphologyStats     – input with new columns written back per object:
%       organelleBoundaryDistancePix – nearest distance to cell boundary (px)
%       organelleBoundaryDistance    – nearest distance to cell boundary (µm)
%       organelleBoundaryRadialIdxList – [x y distance] per perimeter
%                                        sample point (0 where already
%                                        outside the cell mask)
%   d2BoundaryStats     – {nCh x nZ x nT} cell array; same tables as
%                         written into morphologyStats.
%   plotLines           – {nCh x nZ x nT} cell array of [* x 2] overlay lines.
%   plotPoints          – {nCh x nZ x nT} cell array of [* x 2] contact points.

[nC_in, nZ, nT] = size(morphologyStats);
nCh = numel(Cidx);
nY  = imSize(1);
nX  = imSize(2);
[nYb, nXb, bC, bZ, bT] = size(cellBoundary);

% cellBoundary must share the organelle image's pixel grid -- ray-casting
% from an organelle's perimeter to the cell boundary is only meaningful
% if both are the same crop/resample. A mismatch here means cellBoundary
% is stale or was computed from a different processing stage (e.g. not
% re-run after a re-crop); previously this silently corrupted every
% organelle's decoded perimeter via ind2sub against the wrong image size,
% scattering points and reporting distance 0 for many unrelated objects.
if nYb ~= nY || nXb ~= nX
    warning('organelleD2BoundaryCompute:sizeMismatch', ...
        ['cellBoundary is [%d x %d] but the organelle image is [%d x %d] -- ' ...
        'cellBoundary looks stale or was computed from a different ' ...
        'processing stage. organelleBoundaryDistance results will be ' ...
        'wrong until cellBoundary is recomputed to match the current ' ...
        'organelle segmentation.'], nYb, nXb, nY, nX);
end

d2BoundaryStats = cell(nCh, nZ, nT);
plotLines       = cell(nCh, nZ, nT);
plotPoints      = cell(nCh, nZ, nT);

for iT = 1:nT
    for iZ = 1:nZ
        for iCh = 1:nCh
            iC = min(Cidx(iCh), nC_in);

            statsIn = morphologyStats{iC, iZ, iT};
            if isempty(statsIn)
                continue
            end

            % Fill internal holes (nucleus, vacuole, any segmentation gap)
            % before tracing -- this metric means distance to the outer
            % plasma membrane. Without filling, bwboundaries also traces
            % hole edges, and any organelle sitting against a hole (e.g.
            % perinuclear mitochondria, which cluster centrally) gets
            % flagged as touching "outside" and reports 0 even though it
            % may be tens of pixels from the true cell edge.
            mask = imfill(cellBoundary(:,:, min(iCh,bC), min(iZ,bZ), min(iT,bT)), 'holes');
            if ~any(mask(:)) || all(mask(:))
                % no boundary in frame (mask absent, or nothing outside it)
                stats = statsIn;
                nO = height(stats);
                stats.organelleBoundaryDistancePix   = nan(nO,1);
                stats.organelleBoundaryDistance      = nan(nO,1);
                stats.organelleBoundaryRadialIdxList = cell(nO,1);
                d2BoundaryStats{iCh, iZ, iT} = stats;
                morphologyStats{iC, iZ, iT}  = stats;
                plotLines{iCh, iZ, iT}       = zeros(0,2);
                plotPoints{iCh, iZ, iT}      = zeros(0,2);
                continue
            end

            B0 = bwboundaries(mask);
            B1 = cellfun(@(x) [double(x(:,1)) double(x(:,2)); nan nan], B0, 'UniformOutput', false);
            boundaryPoly = cat(1, B1{:});
            boundaryPoly = funcPolylineDropSingletons(boundaryPoly);
            maskVec = mask(:);

            [stats, pLines, pPoints] = planeD2Boundary( ...
                statsIn, boundaryPoly, maskVec, [nY nX], [nYb nXb], calibration, span, radius);

            d2BoundaryStats{iCh, iZ, iT} = stats;
            morphologyStats{iC, iZ, iT}  = stats;
            plotLines{iCh, iZ, iT}       = pLines;
            plotPoints{iCh, iZ, iT}      = pPoints;
        end
    end
end

end % organelleD2BoundaryCompute


% =========================================================================
function [stats, plotLines, plotPoints] = planeD2Boundary( ...
    statsIn, boundaryPoly, maskVec, orgSize, bSize, calibration, span, radius)
%PLANED2BOUNDARY  Single-plane worker — see organelleD2BoundaryCompute for docs.
%
% orgSize decodes the organelle's OWN PerimeterIdxList (must match the
% image the organelle segmentation/regionprops actually ran on); bSize
% indexes into maskVec, which is linearised from the cellBoundary-derived
% mask and so must use cellBoundary's own size. These are deliberately
% kept separate -- collapsing them back into one shared size is exactly
% the bug this function was fixed for.

perimList = statsIn.organellePerimeterIdxList;
pixList   = statsIn.organellePixelIdxList;
nO        = height(statsIn);

stats = statsIn;
stats.organelleBoundaryDistancePix   = nan(nO,1);
stats.organelleBoundaryDistance      = nan(nO,1);
stats.organelleBoundaryRadialIdxList = cell(nO,1);

plotLinesCell  = cell(nO,1);
plotPointsCell = cell(nO,1);

rayOpts          = struct();
rayOpts.maxRange = radius;

% Very small (1-2px), direction-agnostic rescue for the one specific case
% the ray-cast is known to miss: boundary sitting immediately adjacent to
% a boundary point but off-axis from that point's local outward normal
% (same mechanism as organelleD2ErCompute's ER rescue, applied to the
% cell-boundary mask instead). bwdist(~mask2D) gives, for any point
% inside the cell, the distance to the nearest OUTSIDE pixel -- exactly
% the boundary distance. Only ever improves on the ray-cast's own answer
% within this tiny tolerance.
rescueRadiusPix       = 2;
mask2D                = reshape(maskVec, bSize);
[DRescue, idxRescue]  = bwdist(~mask2D);

% Whole-plane label image (every organelle's own pixels marked with its
% own row index, in the organelle image's own coordinate space -- orgSize,
% same space perimeters/pixels are decoded in), for the other-organelle
% obstruction check below.
Lorg = zeros(orgSize, 'int32');
for k = 1:nO
    Lorg(pixList{k}) = k;
end

for iO = 1:nO
    perim = perimList{iO,1};
    if numel(perim) < 3
        plotLinesCell{iO}  = zeros(0,2);
        plotPointsCell{iO} = zeros(0,2);
        continue
    end

    srcContour = contourFromPerimeterIdx(perim, orgSize);

    % Ray/plot density is now always full (one sample per raw perimeter
    % pixel) -- `span` no longer trades off density against tangent
    % stability, it controls only the smoothing window.
    rayOpts.numSamples   = numel(perim);
    rayOpts.smoothSpanPx = max(span, 1);

    [dist, hitPts, sampled] = normalsToNearestIntersectionSpline(srcContour, boundaryPoly, rayOpts);

    sampledLin = sub2ind(bSize, ...
        min(max(round(sampled(:,1)),1),bSize(1)), min(max(round(sampled(:,2)),1),bSize(2)));
    outside = ~maskVec(sampledLin);

    % rescue: boundary within 1-2px that the ray-cast's normal direction
    % missed (only applies to points still inside the cell -- points
    % already outside get distance 0 below regardless).
    rescueDist = double(DRescue(sampledLin));
    useRescue  = ~outside & rescueDist <= rescueRadiusPix & (isnan(dist) | rescueDist < dist);
    if any(useRescue)
        [rescueRow, rescueCol] = ind2sub(bSize, idxRescue(sampledLin(useRescue)));
        dist(useRescue)      = rescueDist(useRescue);
        hitPts(useRescue, 1) = rescueRow;
        hitPts(useRescue, 2) = rescueCol;
    end

    % other-organelle obstruction: the straight line to the winning hit
    % may cut through some OTHER organelle's own bulk on the way -- not a
    % legitimate "unobstructed" answer (skip points already outside the
    % cell, which get distance 0 below regardless).
    obstructedByOrganelle = ~outside & segmentsCrossOtherOrganelles(sampled, hitPts, Lorg, iO);
    dist(obstructedByOrganelle)      = NaN;
    hitPts(obstructedByOrganelle, :) = NaN;

    % points already outside the cell mask: zero distance, self as hit
    dist(outside)      = 0;
    hitPts(outside, :) = sampled(outside, :);

    points = [hitPts(:,2), hitPts(:,1), dist];   % [x y distance]

    stats.organelleBoundaryRadialIdxList{iO,1} = points;
    stats.organelleBoundaryDistancePix(iO,1)   = min(points(:,3));
    stats.organelleBoundaryDistance(iO,1)      = min(points(:,3)) .* calibration;

    plSeg = arrayfun(@(x1,y1,x2,y2) [x1 y1; x2 y2; nan nan], ...
        sampled(:,2), sampled(:,1), points(:,1), points(:,2), 'UniformOutput', false);
    plotLinesCell{iO}  = [nan nan; cat(1, plSeg{:})];
    plotPointsCell{iO} = [nan nan; points(:,1:2)];
end

plotLines  = cat(1, plotLinesCell{:});
plotPoints = cat(1, plotPointsCell{:});

end % planeD2Boundary
