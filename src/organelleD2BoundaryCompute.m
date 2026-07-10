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
% INPUTS
%   morphologyStats – {nC_in x nZ x nT} cell array of organelle stats
%                     tables. Each table must have
%                     .organellePerimeterIdxList.
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
%   span            – perimeter sample-density control (px), same
%                     meaning as in organelleD2ErCompute.
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
nO        = height(statsIn);

stats = statsIn;
stats.organelleBoundaryDistancePix   = nan(nO,1);
stats.organelleBoundaryDistance      = nan(nO,1);
stats.organelleBoundaryRadialIdxList = cell(nO,1);

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

    srcContour = contourFromPerimeterIdx(perim, orgSize);

    if span > 1
        rayOpts.numSamples = max(8, round(numel(perim) / span));
    else
        rayOpts.numSamples = numel(perim);
    end

    [dist, hitPts, sampled] = normalsToNearestIntersectionSpline(srcContour, boundaryPoly, rayOpts);

    % points already outside the cell mask: zero distance, self as hit
    sampledLin = sub2ind(bSize, ...
        min(max(round(sampled(:,1)),1),bSize(1)), min(max(round(sampled(:,2)),1),bSize(2)));
    outside = ~maskVec(sampledLin);
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
