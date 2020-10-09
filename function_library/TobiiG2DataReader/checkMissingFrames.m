function checkMissingFrames(data, toleranceRatio, reportThresholdSec)
%CHECKMISSINGFRAMES looks for missing frames in the eye and scene video
%   This function first tries to estimate the intersample period for the
%   video using the median of intersample differences. It then searches for
%   samples that violate this expectation (+- the toleranceRatio), reporting:
%   a) data for the whole video
%   b) data for faulty periods >= reportThresholdSec (given in seconds)


if isfield(data.video,'eye')
    checkMissingFramesImpl(data.video.eye.fts, 'eye', toleranceRatio, reportThresholdSec);
end
checkMissingFramesImpl(data.video.scene.fts, 'scene', toleranceRatio, reportThresholdSec);




function checkMissingFramesImpl(ts, id, toleranceRatio, reportThresholdSec)
%toleranceRatio = 0.05;
%reportThresholdMs = 0.1;

dt = diff(ts);
ifp = median(dt); % estimate for the interframe period
tolerance = toleranceRatio * ifp;
faulty = dt > ifp + tolerance | dt < ifp - tolerance;
dtFaulty = dt(faulty);
totalFaulty = sum(dtFaulty);
totalVideo = ts(end) - ts(1);
percFaulty = 100 * totalFaulty / totalVideo;
if totalFaulty > 0
    warning('Missing %.2f/%.2f seconds of the %s video (%.2f%%)\n', ...
        totalFaulty, totalVideo, id, percFaulty);
    if any(dtFaulty>reportThresholdSec)
        fprintf('Worst offenders:\n');
        for v = sort(dtFaulty, 'descend')
            if v >= reportThresholdSec
                fprintf('%.3fs (%.2f%%)\n', ...
                    v, 100 * v / totalVideo );
            end
        end
    end
end