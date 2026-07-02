function [poly, segOwner] = otherObjectsPolyline(perimList, imSize, excludeIdx, centroids, maxReach)
%OTHEROBJECTSPOLYLINE  Concatenated [row col] target polyline of every
% object's perimeter except excludeIdx, for use as the target contour in
% normalsToNearestIntersectionSpline (self-exclusion nearest-neighbour
% search).
%
%   [poly, segOwner] = otherObjectsPolyline(perimList, imSize, excludeIdx)
%   [poly, segOwner] = otherObjectsPolyline(perimList, imSize, excludeIdx, centroids, maxReach)
%
% INPUTS
%   perimList  – {nObj x 1} cell of ordered PerimeterIdxList linear
%                indices (as stored in morphologyStats).
%   imSize     – [nY nX] size of the frame perimList indices are linear into.
%   excludeIdx – scalar row index into perimList to exclude (the source
%                object itself).
%   centroids  – (optional) nObj x 2 [x y] centroids, used with maxReach
%                to skip objects that cannot possibly be hit (keeps the
%                target polyline small when there are many objects).
%   maxReach   – (optional) scalar px; an object is dropped from the
%                target if its centroid is farther than maxReach from
%                excludeIdx's centroid. Pass a generous bound (ray
%                maxRange + a couple of object diameters) -- this is a
%                cheap coarse prefilter, not an exact cutoff.
%
% OUTPUTS
%   poly      – Px2 [row col] polyline, NaN-separated between objects.
%   segOwner  – (P-1)x1 owner id (row index into perimList) per segment
%               of poly; NaN for the degenerate segments that touch a
%               NaN separator row (never matched by the caller).

nO = numel(perimList);

useReach = nargin >= 5 && ~isempty(centroids) && ~isempty(maxReach);
if useReach
    d = hypot(centroids(:,1) - centroids(excludeIdx,1), ...
              centroids(:,2) - centroids(excludeIdx,2));
    keep = d <= maxReach;
else
    keep = true(nO, 1);
end
keep(excludeIdx) = false;

polyParts  = cell(nO, 1);
ownerParts = cell(nO, 1);
for k = 1:nO
    if ~keep(k) || isempty(perimList{k})
        continue
    end
    [r, c] = ind2sub(imSize, perimList{k});
    pk = [double(r(:)), double(c(:))];
    if any(pk(1,:) ~= pk(end,:))
        pk = [pk; pk(1,:)]; %#ok<AGROW>
    end
    polyParts{k}  = [pk; nan nan];
    ownerParts{k} = [repmat(k, size(pk,1), 1); nan];
end

poly    = cat(1, polyParts{:});
ptOwner = cat(1, ownerParts{:});

if isempty(poly)
    poly     = zeros(0, 2);
    segOwner = zeros(0, 1);
else
    segOwner = ptOwner(1:end-1);
end
end
