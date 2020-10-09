function data = organizeTobiiGlassesEyeData(pc,pd,gd,gp,gp3)
% pc : pupil center (per eye)
% pd : pupil diameter (per eye)
% gd : gaze direction (per eye)
% gp : gaze position on scene video (binocular)
% gp3: gaze convergence position in 3D space (binocular)
assert(  all(ismember( gp.ts,pc.ts))) % binocular signal should not contain timestamps not in the monocular data, as need data from at least one eye for a binocular signal to be ejected
assert(isempty(setxor(gp3.ts,gp.ts))) % both binocular channels should contain the same timestamps
assert(isempty(setxor( pd.ts,pc.ts))) % monocular channels should contain the same timestamps
assert(isempty(setxor( gd.ts,pc.ts))) % monocular channels should contain the same timestamps

gidx = unique(pc.gidx);
data.left       = struct('gidx',gidx);
data.right      = struct('gidx',gidx);
data.binocular  = struct('gidx',gidx);
nSamp           = length(gidx);

monoc = struct('pc',pc,'pd',pd,'gd',gd);
binoc = struct('gp',gp,'gp3',gp3);

% deal with monocular data
vars = {'pc','pd','gd'};    % both field in monoc struct, and the field in the input, so will index actual data by applying field twice in a row. Looks weird but is correct.
for v=1:length(vars)
    % setup
    qCopyTs = v==1;   % else check
    qLeftEye = monoc.(vars{v}).eye=='l';
    
    for e=1:2
        if e==1
            eye = 'left';
            qEye=  qLeftEye;
        else
            eye = 'right';
            qEye= ~qLeftEye;
        end
        
        % get gidx and match to reference gidx already in output struct
        gidx    = monoc.(vars{v}).gidx(qEye);
        qIdxs   = ismember(data.(eye).gidx,gidx);   % gidx is already ordered
        assert(sum(qIdxs)==sum(qEye))
        
        % deal with ts
        if qCopyTs
            % copy ts field
            data.(eye).ts = nan(nSamp,1);
            data.(eye).ts(qIdxs) = monoc.(vars{v}).ts(qEye);
        else
            % check tss match
            tss = data.(eye).ts(qIdxs) - monoc.(vars{v}).ts(qEye);
            assert(all(tss==0 | isnan(tss)))
            if any(isnan(tss))
                warning('check for matching nans in tss failed, we have at least one sample where not all monocular packets are present for an eye');
            end
        end
        
        % allocate array
        data.(eye).(vars{v}) = nan(nSamp,size(monoc.(vars{v}).(vars{v}),2));
        % copy values to correct positions
        data.(eye).(vars{v})(qIdxs,:) = monoc.(vars{v}).(vars{v})(qEye,:);
        
        % check gidxs match
        assert(all(data.(eye).gidx(qIdxs) == monoc.(vars{v}).gidx(qEye)))
    end
end

% deal with binocular data
vars = {'gp','gp3'};    % both field in monoc struct, and the field in the input, so will index actual data by applying field twice in a row. Looks weird but is correct.
for v=1:length(vars)
    % setup
    qCopyTs = v==1;   % else check
    
    % get gidx and match to reference gidx already in output struct
    gidx    = binoc.(vars{v}).gidx;
    qIdxs   = ismember(data.binocular.gidx,gidx);   % gidx is already ordered
    assert(sum(qIdxs)==length(gidx))
    
    % deal with ts
    if qCopyTs
        % copy ts field
        data.binocular.ts = nan(nSamp,1);
        data.binocular.ts(qIdxs) = binoc.(vars{v}).ts;
    else
        % check tss match
        assert(all(data.binocular.ts(qIdxs) == binoc.(vars{v}).ts))
    end
    
    % allocate array
    data.binocular.(vars{v}) = nan(nSamp,size(binoc.(vars{v}).(vars{v}),2));
    % copy values to correct positions
    data.binocular.(vars{v})(qIdxs,:) = binoc.(vars{v}).(vars{v});
    
    % check gidxs match
    assert(all(data.(eye).gidx(qIdxs) == binoc.(vars{v}).gidx))
end

% for each binocular gidx, check whether monocular data for left and right
% eye was available
qLeftEye = monoc.pc.eye=='l';
il = find( qLeftEye);
ir = find(~qLeftEye);
% check if we have monocular sample at all for each binocular
[hasLeft ,ifl] = ismember(data.binocular.gidx,monoc.pc.gidx(il));
[hasRight,ifr] = ismember(data.binocular.gidx,monoc.pc.gidx(ir));
% for those for which we have a sample, check if there was gaze
hasLeft  = hasLeft  & ~any(isnan(monoc.gd.gd(il(ifl),:)),2);
hasRight = hasRight & ~any(isnan(monoc.gd.gd(ir(ifr),:)),2);
% now count on how many monocular samples each binocular sample is based
data.binocular.nEye = sum([hasLeft hasRight],2);

% check gidx is still monotonically increasing and then remove, it has
% served its purpose
for eye={'left','right','binocular'}
    assert(issorted(data.(eye{1}).gidx,'monotonic'))
    data.(eye{1}) = rmfield(data.(eye{1}),'gidx');
end