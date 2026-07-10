classdef testOrganelleDistance < matlab.unittest.TestCase
%TESTORGANELLEDISTANCE  Unit tests for the OrganelleDistances_sandbox library.
%
% USAGE
%   Run all tests from the MATLAB command window:
%       results = runtests('tests/testOrganelleDistance');
%       table(results)
%
% COVERAGE
%   computeSplineNormals / normalsToNearestIntersectionSpline — 4 tests
%   otherObjectsPolyline / contourFromPerimeterIdx            — 2 tests
%   organelleD2ErCompute                                      — 3 tests
%   organelleD2OrganelleCompute                                — 4 tests
%   organelleD2BoundaryCompute                                 — 4 tests
%
% REQUIREMENTS
%   MATLAB R2019b+ (matlab.unittest framework)
%   Image Processing Toolbox (regionprops, bwtraceboundary, bwboundaries)
%   All src/ and utils/ on the path (added automatically by TestClassSetup).
%   core/utils/funcPolylineDropSingletons.m from the sibling
%   AnalyzERproject_sandbox repo (organelleD2ErCompute/organelleD2BoundaryCompute's
%   only external dependency) -- also added automatically if that sibling
%   repo is found next to this one.

    properties (Constant)
        ImSize = [200, 200]
        DistTol = 3.0   % px -- spline-resampling numerical tolerance
    end

    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'utils'));
            addpath(fullfile(rootDir, 'demos'));

            siblingUtils = fullfile(rootDir, '..', 'AnalyzERproject_sandbox', 'core', 'utils');
            if isfolder(siblingUtils)
                addpath(siblingUtils);
            end
        end
    end

    % =====================================================================
    % Shared synthetic helpers
    % =====================================================================
    methods (Static, Access = private)
        function [labelIm, centers, radii] = threeCircles()
            [X, Y] = meshgrid(1:testOrganelleDistance.ImSize(2), 1:testOrganelleDistance.ImSize(1));
            centers = [60 100; 140 100; 100 20];   % [x y]
            radii   = [15 15 10];
            labelIm = zeros(testOrganelleDistance.ImSize);
            for k = 1:3
                labelIm((X-centers(k,1)).^2 + (Y-centers(k,2)).^2 <= radii(k)^2) = k;
            end
        end

        function statsIn = statsFromLabelIm(labelIm)
            nO = max(labelIm(:));
            S  = regionprops('table', labelIm, 'Centroid');
            perimList = cell(nO, 1);
            pixList   = cell(nO, 1);
            for k = 1:nO
                mask = labelIm == k;
                [r0, c0] = find(mask, 1);
                B = bwtraceboundary(mask, [r0 c0], 'N', 8);
                perimList{k} = sub2ind(size(labelIm), B(:,1), B(:,2));
                pixList{k}   = find(mask);
            end
            statsIn = table();
            statsIn.organelleCentroid         = S.Centroid;
            statsIn.organellePerimeterIdxList = perimList;
            statsIn.organellePixelIdxList     = pixList;
            statsIn.organelleID               = (1:nO)';
        end
    end

    % =====================================================================
    % Ray-cast core
    % =====================================================================
    methods (Test)
        function testSingleTargetKnownGap(tc)
            [labelIm, centers, radii] = tc.threeCircles();
            mask1 = labelIm == 1;
            mask2 = labelIm == 2;
            c1 = bwboundaries(mask1, 8, 'noholes'); c1 = c1{1};
            c2 = bwboundaries(mask2, 8, 'noholes'); c2 = c2{1};

            gapExpected = (centers(2,1) - radii(2)) - (centers(1,1) + radii(1));
            opts = struct('numSamples', 200, 'maxRange', 100, 'chunkSize', 50, 'pad', 5);
            [dist] = normalsToNearestIntersectionSpline(c1, c2, opts);

            tc.verifyEqual(min(dist), gapExpected, 'AbsTol', tc.DistTol);
        end

        function testMultiTargetOwnerId(tc)
            [labelIm, centers, radii] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            perimList = statsIn.organellePerimeterIdxList;
            centroids = statsIn.organelleCentroid;

            [poly, segOwner] = otherObjectsPolyline(perimList, size(labelIm), 1, centroids, 200);
            srcContour = contourFromPerimeterIdx(perimList{1}, size(labelIm));
            opts = struct('numSamples', 200, 'maxRange', 200, 'chunkSize', 50, 'pad', 5, 'segOwner', segOwner);
            [dist, ~, ~, ~, hitOwner] = normalsToNearestIntersectionSpline(srcContour, poly, opts);

            [minD, iMin] = min(dist);
            gapExpected = (centers(2,1) - radii(2)) - (centers(1,1) + radii(1));

            tc.verifyEqual(minD, gapExpected, 'AbsTol', tc.DistTol);
            tc.verifyEqual(hitOwner(iMin), 2);
        end

        function testNoTargetInRangeGivesNaN(tc)
            [labelIm] = tc.threeCircles();
            mask1 = labelIm == 1;
            c1 = bwboundaries(mask1, 8, 'noholes'); c1 = c1{1};
            farAway = [10 10; 10 11; 11 11; 11 10; 10 10];  % tiny closed square, far from c1

            opts = struct('numSamples', 50, 'maxRange', 5, 'chunkSize', 50, 'pad', 1);
            dist = normalsToNearestIntersectionSpline(c1, farAway, opts);
            tc.verifyTrue(all(isnan(dist)));
        end

        function testContourRoundTrip(tc)
            [labelIm] = tc.threeCircles();
            mask1 = labelIm == 1;
            [r0, c0] = find(mask1, 1);
            B = bwtraceboundary(mask1, [r0 c0], 'N', 8);
            lin = sub2ind(size(labelIm), B(:,1), B(:,2));
            c1b = contourFromPerimeterIdx(lin, size(labelIm));
            tc.verifyEqual(c1b, double(B));
        end
    end

    % =====================================================================
    % otherObjectsPolyline
    % =====================================================================
    methods (Test)
        function testOtherObjectsExcludesSelf(tc)
            [labelIm] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            [poly, segOwner] = otherObjectsPolyline(statsIn.organellePerimeterIdxList, size(labelIm), 1, statsIn.organelleCentroid, 1000);
            tc.verifyTrue(all(segOwner(~isnan(segOwner)) ~= 1));
            tc.verifyTrue(any(segOwner == 2));
            tc.verifyTrue(any(segOwner == 3));
        end

        function testOtherObjectsReachPrefilterDrops(tc)
            [labelIm] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            % centroid-to-centroid: obj1->obj2 = 80px, obj1->obj3 = 89.4px.
            % A reach of 85 keeps object 2 but drops object 3.
            [poly, segOwner] = otherObjectsPolyline(statsIn.organellePerimeterIdxList, size(labelIm), 1, statsIn.organelleCentroid, 85); %#ok<ASGLU>
            tc.verifyTrue(any(segOwner == 2));
            tc.verifyFalse(any(segOwner == 3));
        end
    end

    % =====================================================================
    % organelleD2ErCompute
    % =====================================================================
    methods (Test)
        function testD2ErKnownDistance(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;

            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            erStrip = false(imSize);
            erStrip(:, 139:141) = true;
            erImg = reshape(erStrip, imSize(1), imSize(2), 1, 1, 1);

            [~, morphOut] = organelleD2ErCompute(morphologyStats, erImg, erImg, [], {cell(0,1)}, 1, 1, 1, 100);
            stats = morphOut{1};

            expectedDist = 139 - (cx + r);
            tc.verifyEqual(stats.organelleErDistance(1), expectedDist, 'AbsTol', tc.DistTol);
            tc.verifyEqual(stats.organelleErOverlapArea(1), 0);
        end

        function testD2ErOverlap(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;

            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            erStrip = false(imSize);
            erStrip(:, cx-1:cx+1) = true;   % strip runs straight through the object
            erImg = reshape(erStrip, imSize(1), imSize(2), 1, 1, 1);

            [~, morphOut] = organelleD2ErCompute(morphologyStats, erImg, erImg, [], {cell(0,1)}, 1, 1, 1, 100);
            stats = morphOut{1};

            tc.verifyEqual(stats.organelleErDistance(1), 0);
            tc.verifyGreaterThan(stats.organelleErOverlapArea(1), 0);
        end

        function testD2ErRadialVectorMatchesScalar(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;

            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            erStrip = false(imSize);
            erStrip(:, 139:141) = true;
            erImg = reshape(erStrip, imSize(1), imSize(2), 1, 1, 1);

            [~, morphOut] = organelleD2ErCompute(morphologyStats, erImg, erImg, [], {cell(0,1)}, 1, 1, 1, 100);
            stats = morphOut{1};
            radial = stats.organelleRadialIdxList{1};

            tc.verifyEqual(min(radial(:,3)), stats.organelleErDistancePix(1), 'AbsTol', 1e-9);
        end
    end

    % =====================================================================
    % organelleD2BoundaryCompute
    % =====================================================================
    methods (Test)
        function testD2BoundaryKnownDistance(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            % compact cell mask wholly inside the frame (not touching the
            % image border, unlike a naive half-plane mask -- otherwise the
            % image edge itself becomes a spurious "wall").
            cellMask = X >= 10 & X <= 150 & Y >= 10 & Y <= 190;
            cellImg  = reshape(cellMask, imSize(1), imSize(2), 1, 1, 1);

            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            morphOut = organelleD2BoundaryCompute(morphologyStats, imSize, cellImg, 1, 1, 1, 150);
            stats = morphOut{1};

            expectedDist = cx - r - 10;   % nearest wall is the left edge at col=10
            tc.verifyEqual(stats.organelleBoundaryDistance(1), expectedDist, 'AbsTol', tc.DistTol);
        end

        function testD2BoundaryTouching(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cellMask = X >= 10 & X <= 150 & Y >= 10 & Y <= 190;
            cellImg  = reshape(cellMask, imSize(1), imSize(2), 1, 1, 1);

            cx = 145; cy = 100; r = 15;   % spans col 130-160, crosses the right wall at col=150
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            morphOut = organelleD2BoundaryCompute(morphologyStats, imSize, cellImg, 1, 1, 1, 150);
            stats = morphOut{1};

            tc.verifyEqual(stats.organelleBoundaryDistance(1), 0);
        end

        function testD2BoundaryRadialVectorMatchesScalar(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cellMask = X >= 10 & X <= 150 & Y >= 10 & Y <= 190;
            cellImg  = reshape(cellMask, imSize(1), imSize(2), 1, 1, 1);

            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            morphOut = organelleD2BoundaryCompute(morphologyStats, imSize, cellImg, 1, 1, 1, 150);
            stats = morphOut{1};
            radial = stats.organelleBoundaryRadialIdxList{1};

            tc.verifyEqual(min(radial(:,3)), stats.organelleBoundaryDistancePix(1), 'AbsTol', 1e-9);
        end

        function testD2BoundaryStaleCellBoundaryWarnsAndDecodesCorrectly(tc)
            % Regression test for the ind2sub-decode bug: cellBoundary's
            % own array size must NOT be used to decode the organelle's
            % own PerimeterIdxList (that's imSize's job). Pads cellImg to
            % a DIFFERENT array size than imSize while keeping the actual
            % wall at the same absolute pixel coordinates as
            % testD2BoundaryKnownDistance -- before the fix, using
            % cellBoundary's (larger) size to decode the organelle's own
            % perimeter indices via ind2sub scrambled the organelle's
            % contour coordinates and gave a wrong/garbled distance;
            % after the fix, imSize is used for that decode and the
            % result is unaffected by cellBoundary's own array size,
            % matching the known geometric distance exactly.
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            cellMask = X >= 10 & X <= 150 & Y >= 10 & Y <= 190;

            padCols = 50;
            cellMaskPadded = [cellMask, false(imSize(1), padCols)];
            cellImgPadded  = reshape(cellMaskPadded, imSize(1), imSize(2)+padCols, 1, 1, 1);

            cx = 60; cy = 100; r = 15;
            labelIm = zeros(imSize);
            labelIm((X-cx).^2 + (Y-cy).^2 <= r^2) = 1;
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            tc.verifyWarning( ...
                @() organelleD2BoundaryCompute(morphologyStats, imSize, cellImgPadded, 1, 1, 1, 150), ...
                'organelleD2BoundaryCompute:sizeMismatch');

            warnState = warning('off', 'organelleD2BoundaryCompute:sizeMismatch');
            cleanupObj = onCleanup(@() warning(warnState)); %#ok<NASGU>
            morphOut = organelleD2BoundaryCompute(morphologyStats, imSize, cellImgPadded, 1, 1, 1, 150);
            stats = morphOut{1};

            expectedDist = cx - r - 10;
            tc.verifyEqual(stats.organelleBoundaryDistance(1), expectedDist, 'AbsTol', tc.DistTol);
        end
    end

    % =====================================================================
    % organelleD2OrganelleCompute
    % =====================================================================
    methods (Test)
        function testD2OrganelleNearestNeighbour(tc)
            [labelIm, centers, radii] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            [morphOut] = organelleD2OrganelleCompute(morphologyStats, size(labelIm), 1, 1, 1, 100, 20);
            stats = morphOut{1};

            gapExpected = (centers(2,1) - radii(2)) - (centers(1,1) + radii(1));
            tc.verifyEqual(stats.organelleNnDistance(1), gapExpected, 'AbsTol', tc.DistTol);
            tc.verifyEqual(stats.organelleNnID(1), 2);
            tc.verifyEqual(stats.organelleNnID(2), 1);
        end

        function testD2OrganelleRadialVectorMatchesScalar(tc)
            [labelIm] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            [morphOut] = organelleD2OrganelleCompute(morphologyStats, size(labelIm), 1, 1, 1, 100, 20);
            stats = morphOut{1};
            radial = stats.organelleNnRadialIdxList{1};

            tc.verifyEqual(min(radial(:,3)), stats.organelleNnDistancePix(1), 'AbsTol', 1e-9);
        end

        function testD2OrganelleNeighborCountRespondsToRadius(tc)
            [labelIm] = tc.threeCircles();
            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            % gap between obj1/obj2 is ~51px; a neighborRadius well above that
            % should mark them as neighbours, a tiny one should not.
            morphWide  = organelleD2OrganelleCompute(morphologyStats, size(labelIm), 1, 1, 1, 100, 60);
            morphTight = organelleD2OrganelleCompute(morphologyStats, size(labelIm), 1, 1, 1, 100, 5);

            tc.verifyGreaterThan(morphWide{1}.organelleNeighborCount(1), 0);
            tc.verifyEqual(morphTight{1}.organelleNeighborCount(1), 0);
        end

        function testD2OrganelleSingleObjectIsAllNaN(tc)
            imSize = tc.ImSize;
            [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));
            labelIm = zeros(imSize);
            labelIm((X-60).^2 + (Y-100).^2 <= 15^2) = 1;

            statsIn = tc.statsFromLabelIm(labelIm);
            morphologyStats = {statsIn};

            morphOut = organelleD2OrganelleCompute(morphologyStats, imSize, 1, 1, 1, 100, 20);
            stats = morphOut{1};

            tc.verifyTrue(isnan(stats.organelleNnDistance(1)));
            tc.verifyEqual(stats.organelleNeighborCount(1), 0);
        end
    end
end
