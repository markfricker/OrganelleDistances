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
% INPUTS
%   morphologyStats – {nC_in x nZ x nT} cell array of organelle stats
%                     tables. Each table must have
%                     .organellePerimeterIdxList and .organelleCentroid.
%   imSize          – [nY nX]; size of the 2-D image plane the perimeter
%                     linear indices were traced on.
%   Cidx            – [1 x nCh] channel indices into morphologyStats.
%   calibration     – scalar µm per pixel.
%   span            – perimeter sample-density control (px), same
%                     meaning as in organelleD2ErCompute.
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

    if span > 1
        rayOpts.numSamples = max(8, round(numel(perim) / span));
    else
        rayOpts.numSamples = numel(perim);
    end
    rayOpts.segOwner = segOwner;

    [dist, hitPts, sampled, ~, hitOwner] = normalsToNearestIntersectionSpline(srcContour, targetPoly, rayOpts);

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
