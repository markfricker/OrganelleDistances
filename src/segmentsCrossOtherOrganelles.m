function obstructed = segmentsCrossOtherOrganelles(srcPts, dstPts, L, excludeLabels)
%SEGMENTSCROSSOTHERORGANELLES  Flag straight-line paths that pass through
% another organelle's own pixels along the way.
%
%   obstructed = segmentsCrossOtherOrganelles(srcPts, dstPts, L, excludeLabels)
%
% A ray-cast (or the small local rescue) only checks whether its endpoint
% is a valid target -- it has no notion of what lies *between* the source
% and that endpoint. For a "distance to nearest ER/organelle/boundary"
% metric meant to reflect unobstructed access, a straight line that
% happens to cut through some other organelle's own bulk on the way isn't
% a legitimate answer, even though the endpoint itself is fine.
%
% INPUTS
%   srcPts, dstPts – Nx2 [row col]; the candidate segment for each sample
%                    point (source boundary point -> the winning hit).
%   L              – label image, 0 = background, k = organelle k's own
%                    pixels (every organelle in the plane, not just the
%                    one being queried).
%   excludeLabels  – Nx1 or Nx2; label(s) to NOT count as an obstruction
%                    for each row. Always include the querying organelle's
%                    own label (a ray starts right at its own boundary).
%                    For organelle-to-organelle distance, also include the
%                    per-point target organelle's own label, since
%                    reaching it necessarily touches its own pixels at the
%                    endpoint. Pass 0 in an unused column (background is
%                    already excluded regardless).
%
% OUTPUT
%   obstructed     – Nx1 logical; true where the straight line crosses a
%                    labeled pixel not in excludeLabels(i,:).

[nY, nX] = size(L);
N = size(srcPts, 1);
obstructed = false(N, 1);

valid = all(isfinite(srcPts), 2) & all(isfinite(dstPts), 2);
if ~any(valid)
    return
end

segLen = hypot(dstPts(:,1) - srcPts(:,1), dstPts(:,2) - srcPts(:,2));
maxLen = max(segLen(valid));
if isempty(maxLen) || ~isfinite(maxLen) || maxLen < 1
    return
end

% One sample per pixel of the longest segment in this batch -- shorter
% segments are oversampled (harmless, just some repeated/coincident
% points) rather than undersampled (which could miss a thin obstruction).
nSteps = max(2, ceil(maxLen));
tt = linspace(0, 1, nSteps);   % 1 x nSteps

rr = round(srcPts(:,1) + (dstPts(:,1) - srcPts(:,1)) .* tt);   % N x nSteps
cc = round(srcPts(:,2) + (dstPts(:,2) - srcPts(:,2)) .* tt);

inBounds = rr >= 1 & rr <= nY & cc >= 1 & cc <= nX;
rrClamped = min(max(rr, 1), nY);
ccClamped = min(max(cc, 1), nX);
lin = sub2ind([nY nX], rrClamped, ccClamped);
labelsAlong = reshape(L(lin), size(lin));

notExcluded = labelsAlong ~= 0;
for k = 1:size(excludeLabels, 2)
    notExcluded = notExcluded & labelsAlong ~= excludeLabels(:, k);
end

crosses = inBounds & notExcluded;
obstructed = valid & any(crosses, 2);

end
