classdef testOrganelleShuffleTest < matlab.unittest.TestCase
%TESTORGANELLESHUFFLETEST  Unit tests for organelleShuffleTest and
% shuffleObjectsInWindow (object-shuffle Monte Carlo test of organelle-to-
% target proximity).
%
% USAGE
%   results = runtests('tests/testOrganelleShuffleTest');
%   table(results)
%
% COVERAGE
%   shuffleObjectsInWindow  - shape/area preserved, inside window, no
%                             overlap, dihedral, exclusion respected
%   organelleShuffleTest    - attraction detected, avoidance detected,
%                             null calibration (false-positive rate),
%                             empty target, per-cell separation +
%                             restrictTargetToCell, rngSeed reproducibility

    properties (Constant)
        ImSize = [120 120]
        Cal    = 0.1      % um/px
    end

    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'utils'));
        end
    end

    % =====================================================================
    % Synthetic helpers
    % =====================================================================
    methods (Static, Access = private)
        function cellMask = diskCell(center, radius)
            sz = testOrganelleShuffleTest.ImSize;
            [X, Y] = meshgrid(1:sz(2), 1:sz(1));
            cellMask = (X-center(1)).^2 + (Y-center(2)).^2 <= radius^2;
        end

        function tgt = verticalLine(col, rows)
            tgt = false(testOrganelleShuffleTest.ImSize);
            tgt(rows, col) = true;
        end

        function pix = squares(topLefts, side)
            % topLefts: [nObj x 2] [row col]
            sz  = testOrganelleShuffleTest.ImSize;
            pix = cell(size(topLefts,1), 1);
            for k = 1:size(topLefts,1)
                [R, C] = ndgrid(topLefts(k,1) + (0:side-1), topLefts(k,2) + (0:side-1));
                pix{k} = sub2ind(sz, R(:), C(:));
            end
        end

        function stats = statsTable(pix, cellIDs)
            stats = table(pix(:), cellIDs(:), ...
                'VariableNames', {'organellePixelIdxList','cellID'});
        end

        function [summary, curves] = runOne(pix, cellMask, tgt, nSim, opts)
            if nargin < 5, opts = struct(); end
            stats = testOrganelleShuffleTest.statsTable(pix, ones(numel(pix),1));
            [curvesOut, summaryOut] = organelleShuffleTest({stats}, tgt, ...
                double(cellMask), 1, testOrganelleShuffleTest.Cal, nSim, 'synthetic', opts);
            summary = summaryOut{1};
            curves  = curvesOut{1};
        end
    end

    % =====================================================================
    % shuffleObjectsInWindow
    % =====================================================================
    methods (Test)
        function shuffleKeepsShapeInsideWindowNoOverlap(tc)
            win = testOrganelleShuffleTest.diskCell([60 60], 50);
            pix = testOrganelleShuffleTest.squares([30 30; 40 70; 70 50; 80 80], 4);
            pix{end+1} = pix{1}(1:6);   % an L-ish irregular object
            for dihedral = [false true]
                rng(1);
                [placed, fb] = shuffleObjectsInWindow(pix, win, struct('dihedral', dihedral));
                tc.verifyFalse(any(fb));
                tc.verifyEqual(cellfun(@numel, placed), cellfun(@numel, pix));
                all_ = cat(1, placed{:});
                tc.verifyTrue(all(win(all_)), 'placed pixel outside window');
                tc.verifyEqual(numel(unique(all_)), numel(all_), 'objects overlap');
                if ~dihedral
                    % translation only: same shape up to an offset
                    for k = 1:numel(pix)
                        [r0, c0] = ind2sub(size(win), sort(pix{k}));
                        [r1, c1] = ind2sub(size(win), sort(placed{k}));
                        d = [r1 c1] - [r0 c0];
                        tc.verifyEqual(size(unique(d, 'rows'), 1), 1);
                    end
                end
            end
        end

        function shuffleRespectsExclusion(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 50);
            excl = false(size(cellMask));
            excl(:, 1:60) = true;                     % left half excluded
            pix = testOrganelleShuffleTest.squares([30 30; 50 30; 70 30], 3);
            rng(2);
            for s = 1:20
                placed = shuffleObjectsInWindow(pix, cellMask & ~excl);
                all_ = cat(1, placed{:});
                tc.verifyFalse(any(excl(all_)));
            end
        end

        % =================================================================
        % organelleShuffleTest
        % =================================================================
        function attractionDetected(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = testOrganelleShuffleTest.verticalLine(60, 10:110);
            rows = (20:8:100)';
            pix = testOrganelleShuffleTest.squares([rows, 62*ones(size(rows))], 3);
            [S, C] = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 99, struct('rngSeed', 3));
            tc.verifyLessThanOrEqual(S.sdi, 0.025);
            tc.verifyLessThanOrEqual(S.pCloser, 0.02);
            tc.verifyGreaterThan(S.maxDeviation, 0);
            tc.verifyLessThan(S.meanDistObs, S.meanDistNull);
            tc.verifyEqual(height(C), 200);
            tc.verifyTrue(all(C.cdfNullLo <= C.cdfNullHi));
        end

        function avoidanceDetected(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = testOrganelleShuffleTest.verticalLine(60, 10:110);
            rows = (45:6:75)';
            pix = testOrganelleShuffleTest.squares( ...
                [rows, 12*ones(size(rows)); rows, 106*ones(size(rows))], 3);
            S = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 99, struct('rngSeed', 4));
            tc.verifyGreaterThanOrEqual(S.sdi, 0.975);
            tc.verifyLessThanOrEqual(S.pFurther, 0.02);
            tc.verifyLessThan(S.maxDeviation, 0);
        end

        function touchingFractionReported(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = testOrganelleShuffleTest.verticalLine(60, 10:110);
            rows = (20:8:100)';
            pix = testOrganelleShuffleTest.squares([rows, 59*ones(size(rows))], 3); % straddle line
            S = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 39, struct('rngSeed', 5));
            tc.verifyEqual(S.fracTouchingObs, 1);
            tc.verifyEqual(S.meanDistObs, 0);
            tc.verifyLessThan(S.fracTouchingNull, 0.5);
        end

        function nullFalsePositiveRateControlled(tc)
            % Objects placed by the null model itself: p-values should be
            % ~uniform, so rejections at 0.05 should be rare.
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = testOrganelleShuffleTest.verticalLine(60, 10:110);
            template = testOrganelleShuffleTest.squares([30 30; 40 70; 70 50; 80 80; 50 40; 60 90; 90 60; 35 55], 3);
            rng(6);
            nRep = 40;
            rejCloser = 0; rejFurther = 0; sdis = zeros(nRep,1);
            for rep = 1:nRep
                pix = shuffleObjectsInWindow(template, cellMask);
                S = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 39);
                rejCloser  = rejCloser  + (S.pCloser  <= 0.05);
                rejFurther = rejFurther + (S.pFurther <= 0.05);
                sdis(rep)  = S.sdi;
            end
            tc.verifyLessThanOrEqual(rejCloser  / nRep, 0.2);
            tc.verifyLessThanOrEqual(rejFurther / nRep, 0.2);
            tc.verifyGreaterThan(mean(sdis), 0.3);
            tc.verifyLessThan(mean(sdis), 0.7);
        end

        function emptyTargetGivesNaNRow(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = false(testOrganelleShuffleTest.ImSize);
            pix = testOrganelleShuffleTest.squares([30 30; 60 60], 3);
            [S, C] = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 19);
            tc.verifyEqual(height(S), 1);
            tc.verifyEqual(S.nObjects, 2);
            tc.verifyTrue(isnan(S.sdi));
            tc.verifyEmpty(C);
        end

        function perCellSeparationAndTargetRestriction(tc)
            c1 = testOrganelleShuffleTest.diskCell([30 60], 25);
            c2 = testOrganelleShuffleTest.diskCell([90 60], 25);
            label = double(c1) + 2*double(c2);
            tgt = testOrganelleShuffleTest.verticalLine(30, 40:80);   % in cell 1 only
            pix = [testOrganelleShuffleTest.squares([50 32; 60 32; 70 32], 3); ...
                   testOrganelleShuffleTest.squares([50 90; 60 90; 70 90], 3)];
            stats = testOrganelleShuffleTest.statsTable(pix, [1 1 1 2 2 2]);

            [curvesOut, summaryOut] = organelleShuffleTest({stats}, tgt, label, 1, ...
                testOrganelleShuffleTest.Cal, 49, 'f', struct('rngSeed', 7));
            S = summaryOut{1};
            tc.verifyEqual(S.cellID, [1; 2]);
            tc.verifyEqual(S.nObjects, [3; 3]);
            tc.verifyFalse(isnan(S.sdi(1)));
            tc.verifyTrue(isnan(S.sdi(2)), 'cell 2 must not see cell 1''s target');
            tc.verifyEqual(unique(curvesOut{1}.cellID), 1);

            % unrestricted: cell 2 now measures against cell 1's line
            [~, summaryOut] = organelleShuffleTest({stats}, tgt, label, 1, ...
                testOrganelleShuffleTest.Cal, 49, 'f', ...
                struct('rngSeed', 7, 'restrictTargetToCell', false));
            tc.verifyFalse(isnan(summaryOut{1}.sdi(2)));
        end

        function rngSeedReproducible(tc)
            cellMask = testOrganelleShuffleTest.diskCell([60 60], 55);
            tgt = testOrganelleShuffleTest.verticalLine(60, 10:110);
            pix = testOrganelleShuffleTest.squares([30 30; 40 70; 70 50], 3);
            S1 = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 29, struct('rngSeed', 8, 'dihedral', true));
            S2 = testOrganelleShuffleTest.runOne(pix, cellMask, tgt, 29, struct('rngSeed', 8, 'dihedral', true));
            tc.verifyEqual(S1, S2);
        end

        function progressCallbackCanCancel(tc)
            c1 = testOrganelleShuffleTest.diskCell([30 60], 25);
            c2 = testOrganelleShuffleTest.diskCell([90 60], 25);
            label = double(c1) + 2*double(c2);
            tgt = testOrganelleShuffleTest.verticalLine(30, 40:80) | ...
                  testOrganelleShuffleTest.verticalLine(90, 40:80);
            pix = [testOrganelleShuffleTest.squares([50 32; 60 32], 3); ...
                   testOrganelleShuffleTest.squares([50 92; 60 92], 3)];
            stats = testOrganelleShuffleTest.statsTable(pix, [1 1 2 2]);
            calls = 0;
            function cancel = cb(~, ~)
                calls = calls + 1;
                cancel = true;
            end
            [~, summaryOut] = organelleShuffleTest({stats}, tgt, label, 1, ...
                testOrganelleShuffleTest.Cal, 9, 'f', struct('progressFcn', @cb));
            tc.verifyEqual(calls, 1);
            tc.verifyEqual(height(summaryOut{1}), 1);
        end
    end
end
