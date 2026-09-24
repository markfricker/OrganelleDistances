function [placedIdx, isFallback] = shuffleObjectsInWindow(pixIdx, windowMask, opts)
%SHUFFLEOBJECTSINWINDOW  Randomly re-place a set of pixel objects inside a
% window, preserving each object's shape and preventing mutual overlap.
%
%   [placedIdx, isFallback] = shuffleObjectsInWindow(pixIdx, windowMask, opts)
%
% The null model behind organelleShuffleTest (and DiAna's "shuffle",
% Gilles et al. 2017, Methods 115:55-64): each object keeps its exact
% pixel shape but is translated to a uniformly random position, accepted
% only if every one of its pixels lands inside WINDOWMASK and on no pixel
% already claimed by a previously placed object. Objects are placed in a
% random order each call, so no object systematically gets first pick of
% the free space.
%
% Each candidate position is drawn by picking a random window pixel for
% the object's (rounded) centroid, so placements are uniform over the
% window for small objects; large objects are necessarily restricted to
% where they physically fit, exactly as they would be in the real cell.
%
% INPUTS
%   pixIdx      - {nObj x 1} cell of linear pixel indices into
%                 windowMask's array (the objects' observed positions).
%   windowMask  - [nY x nX] logical; the region objects may occupy.
%   opts        - (optional) struct:
%     .dihedral  - also apply one of the 8 exact 90-degree rotations /
%                  mirror flips to each object (default false). Pixel-exact,
%                  so object area and shape are preserved perfectly --
%                  unlike an arbitrary-angle rotation, which would need
%                  resampling.
%     .maxTries  - candidate positions tried per object before giving up
%                  (default 1000).
%     .batchSize - candidate positions tested per vectorised batch
%                  (default 64; performance only).
%
% OUTPUTS
%   placedIdx  - {nObj x 1} cell of linear indices of each object's new
%                position (same order as pixIdx).
%   isFallback - [nObj x 1] logical; true where no valid position was
%                found within maxTries and the object was left at its
%                observed position instead (typically only for an object
%                too large to fit anywhere else in a crowded window).

if nargin < 3 || isempty(opts)
    opts = struct();
end
dihedral  = getOpt(opts, 'dihedral', false);
maxTries  = getOpt(opts, 'maxTries', 1000);
batchSize = getOpt(opts, 'batchSize', 64);

[nY, nX]  = size(windowMask);
winIdx    = find(windowMask);
nWin      = numel(winIdx);
nObj      = numel(pixIdx);

placedIdx  = pixIdx(:);
isFallback = false(nObj, 1);
occupied   = false(nY, nX);

if nWin == 0
    isFallback(:) = true;
    return
end

for k = randperm(nObj)
    idx0 = pixIdx{k}(:);
    if isempty(idx0)
        continue
    end
    [r, c] = ind2sub([nY nX], idx0);
    offs   = [r c] - round(mean([r c], 1));
    if dihedral
        offs = dihedralTransform(offs, randi(8));
    end

    placed = false;
    tries  = 0;
    while ~placed && tries < maxTries
        nb = min(batchSize, maxTries - tries);
        tries = tries + nb;

        [ar, ac] = ind2sub([nY nX], winIdx(randi(nWin, nb, 1)));
        R = offs(:,1) + ar';            % [nPix x nb]
        C = offs(:,2) + ac';
        inBounds = all(R >= 1 & R <= nY & C >= 1 & C <= nX, 1);
        R(:, ~inBounds) = 1;            % dummy, masked out below
        C(:, ~inBounds) = 1;
        cand = R + (C - 1) * nY;
        ok = inBounds & all(windowMask(cand) & ~occupied(cand), 1);

        j = find(ok, 1);
        if ~isempty(j)
            placedIdx{k} = cand(:, j);
            placed = true;
        end
    end

    if ~placed
        isFallback(k) = true;           % leave at observed position
    end
    occupied(placedIdx{k}) = true;
end

end % shuffleObjectsInWindow


% =========================================================================
function offs = dihedralTransform(offs, t)
%DIHEDRALTRANSFORM  One of the 8 pixel-exact rotations/reflections.
if t > 4
    offs = offs(:, [2 1]);              % transpose (mirror about diagonal)
end
switch mod(t - 1, 4)
    case 1, offs(:,1) = -offs(:,1);
    case 2, offs(:,2) = -offs(:,2);
    case 3, offs = -offs;
end
end % dihedralTransform


% =========================================================================
function v = getOpt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end % getOpt
