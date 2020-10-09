function data = fillMissingSamples(data,expectedFs)
% NB!! we have data for each eye and binocular for matched gidx (gidx by
% now removed from data as no longer needed). Data for same gidx doesn't
% however always come with exactly the same timestamps. As differences are
% (should be) tiny, we can ignore this. check differences are tiny though.
u=unique(data.left.ts-data.right.ts);
assert(all(isnan(u) | u<1000))  % arbitrarily decide that less than one ms is small
u=unique(data.left.ts-data.binocular.ts);
assert(all(isnan(u) | u<1000))  % arbitrarily decide that less than one ms is small

% fill gaps with nans (just split time interval between two samples
% into equally sized bits). You may think we can do this based on gidx.
% every non-existent gidx just needs a ts (or better yet, just keep those
% ts). Can't do that as the gidx are a bit messy. there are some spurious
% ones in there with non-zero s (crap). And there are genuine gaps in the
% data where gidx doesn't increase. So best we can do is just look at the
% data and go from there. I see there are some very few samples where
% there is e.g. 13 ms between samples instead of 10. we'll just have to
% live with that

% gaps are places where ISI is 1.333 (4/3) times longer than expected given
% sampling frequency expectedFs (arbitrarily chosen). Plug 'em, filling the
% signals up with nan data
thr = round(1000*1000/expectedFs*4/3); % consider a gap as more than thr time elapsed between consequtive samples (4/3 means gap when ISI is 1.333 times longer than expected)
thr2= round(1000*1000/expectedFs*3/4);
for c=1:3
    switch c
        case 1
            ch = 'left';
        case 2
            ch = 'right';
        case 3
            ch = 'binocular';
    end
    
    % find all gaps
    dt      = diff(data.(ch).ts);
    iGap    = find(dt>thr);    % index indicates sample _after_ which there is a gap
    
    if any(dt<thr2)
        warning('got some samples closer together than expected.')
    end
    
    % determine how long the new signal will be with gaps filled
    gapSzs          = dt(iGap);
    nSampMissing    = round(gapSzs/(1000*1000/expectedFs))-1;     % round instead of ceil or floor gives smallest deviation from nominal framerate: 5.4->5, 5.6->6
    
    % place samples in the right places (effectively inserts the missing
    % samples)
    % indicate where real samples should be in timeline
    idxs            = ones(1,length(data.(ch).ts));
    idxs(iGap+1)    = nSampMissing+1;
    idxs            = cumsum(idxs);
    % put them there, couched in nan
    fields = fieldnames(data.(ch));
    for f=1:length(fields)
        temp                    = nan(idxs(end),size(data.(ch).(fields{f}),2));
        temp(idxs,:)            = data.(ch).(fields{f});
        data.(ch).(fields{f})   = temp;
    end
    
    % fill gaps in time with faked equally intersecting intervals
    data.(ch).ts = round(interp1(idxs,data.(ch).ts(idxs),1:idxs(end),'linear')).';
end
% tiny differences between timestamps for left and right eye have been
% spotted in the wild, and are not caused by this code: they were indeed in
% the original json file.
% arbitrarily decide that differences less than 50 microseconds are not
% worth reporting. If anyone spots anything larger, i'd love to have a look
% at it.
if (~isequal(size(data.left.ts),size(data.right.ts)) || any(abs(data.left.ts-data.right.ts)>50)) || (~isequal(size(data.left.ts),size(data.binocular.ts)) || any(abs(data.left.ts-data.binocular.ts)>50))
    warning('timestamps for the different eyes are not the same. please contact dcnieho@gmail.com if you are willing to share this recording, I would love to see it and check if this is a problem')
end