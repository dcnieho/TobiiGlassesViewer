function data = readG2DataFiles(recordingDir,qDEBUG)

% Cite as: Niehorster, D.C., Hessels, R.S., and Benjamins, J.S. (2020).
% GlassesViewer: Open-source software for viewing and analyzing data from
% the Tobii Pro Glasses 2 eye tracker. Behavior Research Methods. doi:
% 10.3758/s13428-019-01314-1

% set file format version. cache files older than this are overwritten with
% a newly generated cache file
fileVersion = 10;

if ~isempty(which('matlab.internal.webservices.fromJSON'))
    jsondecoder = @matlab.internal.webservices.fromJSON;
elseif ~isempty(which('jsondecode'))
    jsondecoder = @jsondecode;
else
    error('Your MATLAB version does not provide a way to decode json (which means its really old), upgrade to something newer');
end

qGenCacheFile = ~exist(fullfile(recordingDir,'livedata.mat'),'file');
if ~qGenCacheFile
    % we have a cache file, check its file version
    cache = load(fullfile(recordingDir,'livedata.mat'),'fileVersion');
    qGenCacheFile = cache.fileVersion~=fileVersion;
end

if qGenCacheFile
    % 0 get info about participant and recording
    fid = fopen(fullfile(recordingDir,'recording.json'),'rt');
    recording = jsondecoder(fread(fid,inf,'*char').');
    fclose(fid);
    fid = fopen(fullfile(recordingDir,'participant.json'),'rt');
    participant = jsondecoder(fread(fid,inf,'*char').');
    fclose(fid);
    expectedFs = round(recording.rec_et_samples/recording.rec_length/50)*50;    % find nearest 50Hz
    if qDEBUG
        fprintf('determined fs: %d Hz\n',expectedFs);
    end
    % 1 per segment, read data
    segments = FolderFromFolder(fullfile(recordingDir,'segments'));
    txt = '';
    for s=1:length(segments)
        % read in all segments, just concat
        % 1.1 unpack gz file, if doesn't exist
        if ~exist(fullfile(recordingDir,'segments',segments(s).name,'livedata.json'),'file')
            gunzip(fullfile(recordingDir,'segments',segments(s).name,'livedata.json.gz'));
        end
        % 1.2 read all lines
        fprintf('reading: %s\n',fullfile(recordingDir,'segments',segments(s).name,'livedata.json'));
        fid = fopen(fullfile(recordingDir,'segments',segments(s).name,'livedata.json'),'rt');
        txt = [txt fread(fid,inf,'*char').']; %#ok<AGROW>
        fclose(fid);
    end
    % 2 get ts, s and type per packet
    types = regexp(txt,'(?:"ts":)(\d+,).*?(?:"s":)(\d+,).*?(?:")(pc|pd|gd|gp|gp3|gy|ac|vts|evts|pts|epts|sig|type)(?:")','tokens');   % only get first item of status code, we just care if 0 or not. get delimiter of ts on purpose, so we can concat all matches and just sscanf
    types = cat(1,types{:});
    % 3 transform packets into struct of arrays we do our further
    % processing on
    % fields:
    % - ts
    % - type
    % - dat (whatever is in the type field)
    % - gidx or nan
    % - eye (first letter) or x
    % corresponding packet
    % in the process, remove all packets with non-zero status field (-means
    % its crap, failed or otherwise shouldn't be used), split off and
    % remove fields requiring custom processing (sig, type) and remove
    % unwanted fields (pts, epts)
    qSplit  = strcmp(types(:,3),'sig')|strcmp(types(:,3),'type');
    qUnwanted = strcmp(types(:,3),'pts')|strcmp(types(:,3),'epts');
    qKeep   = ~qSplit & ~qUnwanted;  % remove non-zero s, split off fields and unwanted fields
    nKeep   = sum(qKeep);
    dat     = struct('ts',[],'dat',[],'gidx',[],'eye',[],'s',[]);
    dat.ts  = sscanf(cat(2,types{qKeep,1}),'%f,');
    dat.qValid = sscanf(cat(2,types{qKeep,2}),'%f,')==0;
    dat.dat = nan(nKeep,3);
    dat.gidx= nan(nKeep,1);
    dat.s   = nan(nKeep,1);
    dat.eye = repmat('x',nKeep,1);
    % check where each packet begins and ends
    iNewLines = [1 find(txt==newline)];
    if iNewLines(end)~=length(txt)
        iNewLines = [iNewLines length(txt)];
    end
    assert(length(iNewLines)==length(qKeep)+1,'must have missed a type in the above regexp, or the assumption that each element has only one of the above strings is no longer true')
    % split off the ones we process separately
    q = bounds2bool(iNewLines([qSplit; false])+1,iNewLines([false; qSplit]));
    syncPortAPISyncTxt  = txt(q);
    syncPortAPISyncType = types(qSplit,3);
    % get types for main packets
    types   = types(qKeep,3);
    % now remove these packets from string
    qRem    = ~qKeep;
    q       = bounds2bool(iNewLines([qRem; false])+1,iNewLines([false; qRem]));
    txt(q)  = [];
    clear qSplit qUnwanted qKeep qRem iNewLines q
    % find where each type is in this data
    qData    = [strcmp(types,'pc')  strcmp(types,'pd') strcmp(types,'gd')...
                strcmp(types,'gp')  strcmp(types,'gp3')...
                strcmp(types,'gy')  strcmp(types,'ac')...
                strcmp(types,'vts') strcmp(types,'evts')];
    qHasEye  = any(qData(:,1:3),2);
    qHasGidx = any([qHasEye qData(:,4:5)],2);
    qHasEyeVideo = any(qData(:,9));
    % parse in gidx, eye and s
    gidx = regexp(txt,'(?<="gidx":)\d+,','match');
    assert(~isempty(gidx),'The data file does not fulfill the requirements of this code, the ''gidx'' field is missing. Possibly the firmware of the recording unit was too old.')
    dat.gidx(qHasGidx) = sscanf(cat(2,gidx{:}),'%f,');
    eye  = regexp(txt,'(?<="eye":")[lr]','match');
    dat.eye(qHasEye) = cat(1,eye{:});
    qLeftEye = qHasEye & dat.eye=='l';
    s    = regexp(txt,'(?<="s":)\d+,','match');
    dat.s   = sscanf(cat(2,s{:}),'%d,');
    clear gidx eye s types qHasEye
    % parse in scalar data (only pd)
    dat.dat(qData(:,2),1  )   = parseTobiiGlassesData(txt,  'pd',1);
    % parse in 2-vector (only gp)
    dat.dat(qData(:,4),1:2)   = parseTobiiGlassesData(txt,  'gp',3);
    % parse in 3-vector (other three)
    dat.dat(qData(:,1), : )   = parseTobiiGlassesData(txt,  'pc',4);
    dat.dat(qData(:,3), : )   = parseTobiiGlassesData(txt,  'gd',4);
    dat.dat(qData(:,5), : )   = parseTobiiGlassesData(txt, 'gp3',4);
    % parse in gy and ac
    dat.dat(qData(:,6), : )   = parseTobiiGlassesData(txt,  'gy',4);
    dat.dat(qData(:,7), : )   = parseTobiiGlassesData(txt,  'ac',4);
    % parse in vts and evts
    dat.dat(qData(:,8),1  )   = parseTobiiGlassesData(txt, 'vts',2);
    if qHasEyeVideo
        dat.dat(qData(:,9),1  )   = parseTobiiGlassesData(txt,'evts',2);
    end
    % set data to nan if status code is non-zero
    dat.dat(~~dat.s,:) = nan;
    clear txt
    % 4 organize into types
    % the overall strategy to deal with crap in the files is to:
    % 1. completely remove all gidx for which we have an unexpected number
    %    of packets (~=8).
    % 2. remove all gidx with more than 8 valid packets or less than 3
    % 3. remove binocular data if data for both eyes is incomplete
    % 4. check if any monocular gidx with 2 packets for single eye,
    %    remove
    % 5. remove data where (when sorted by gidx) time apparently went
    %    backward (check per eye if monocular). This mostly gets a few
    %    packets of monocular data where only a single eye is available
    %    for a given gidx (this one is done in
    %    getDataTypeFromTobiiArray)
    minGidx     = min(dat.gidx(qHasGidx));
    maxGidx     = max(dat.gidx(qHasGidx));
    allGidx     = [minGidx:maxGidx].';
    gidxCount   = accumarray(dat.gidx(qHasGidx)-minGidx+1,true(sum(qHasGidx),1));
    validCount  = accumarray(dat.gidx(qHasGidx)-minGidx+1,dat.qValid(qHasGidx));
    % 4.1 completely remove gidx for which there are not the right numbers
    % of packets
    qRemove     = ismember(dat.gidx,allGidx(gidxCount~=8));
    qData   (qRemove,:) = [];
    qLeftEye(qRemove,:) = [];
    dat = replaceElementsInStruct(dat,qRemove,[],[],1);
    % 4.2 first some global processing
    % 4.2.1 remove all gidx with more than 8 valid packets or less than 3 (should have 3 monocular for each eye and 2 binocular if recording succeeded for both eyes. due to bug, sometimes multiple monocular packets for single eye, resp. calculated based on single camera and stereo views, without knowing which is the right one)
    qBad    = validCount>8 | validCount<3;    % less than 3 packets, we probably have a lonely pc while even pd and such failed, better ignore it, must be crap
    qRemove = ismember(dat.gidx,allGidx(qBad,1));
    dat.dat   (qRemove,:) = nan;
    dat.qValid(qRemove  ) = false;
    % 4.2.2 now for each gidx, build a table flagging what we have, so we
    % can throw out all broken samples
    % per gidx, see which valid samples we have
    qValidPerGidx   = [qData(:,1:3)&qLeftEye&dat.qValid qData(:,1:3)&~qLeftEye&dat.qValid qData(:,4:5)&dat.qValid]; % three fields left eye, three fields right, two binocular
    [i,i2]          = find(qValidPerGidx.');   % get column
    qValidPerGidx   = false(maxGidx-minGidx+1,size(qValidPerGidx,2));
    qValidPerGidx(sub2ind(size(qValidPerGidx),dat.gidx(i2)-minGidx+1,i)) = true;
    % 4.2.3 remove binocular data for gidx for which monocular data for
    % both eyes is incomplete (one eye incomplete is possible, binocular
    % data is then computed by the eye tracker based on assumption of
    % unchanged vergence distance since last available true binocular data)
    qBad    = ~all(qValidPerGidx(:,1:3),2) & ~all(qValidPerGidx(:,4:6),2) & any(qValidPerGidx(:,7:8),2);
    qRemove = ismember(dat.gidx,allGidx(qBad))&any(qData(:,4:5),2);   % remove binocular data for these
    % 4.2.4 remove monocular data for incomplete eyes
    qBad    = ~all(qValidPerGidx(:,1:3),2);    % left eye
    qRemove = qRemove | ismember(dat.gidx,allGidx(qBad))&any(qData(:,1:3),2)& qLeftEye;
    qBad    = ~all(qValidPerGidx(:,4:6),2);    % right eye
    qRemove = qRemove | ismember(dat.gidx,allGidx(qBad))&any(qData(:,1:3),2)&~qLeftEye;
    % now remove all these flagged data
    dat.dat   (qRemove,:) = nan;
    dat.qValid(qRemove  ) = false;
    % 4.2.5 check for case with monocular data twice from same eye for
    % given gidx. like for gidx with more than 8 packets, we don't know
    % which is the right one, so remove both
    % To check: sort samples by e, use that to sort gidx and see if any
    % same numbers in a row. those are an issue as that means same gidx
    % multiple times for same eye
    qSel  = qData(:,1)&dat.qValid;
    [~,i] = sort(dat.eye(qSel));
    gidx  = dat.gidx(qSel);
    gidx  = gidx(i);
    qDupl = diff(gidx)==0;
    qRemove = ismember(dat.gidx,gidx(qDupl));
    % now remove all these flagged data
    dat.dat   (qRemove,:) = nan;
    dat.qValid(qRemove  ) = false;
    clear qValidPerGidx qSel qBad qDupl qHasGidx qLeftEye qRemove allGidx gidx gidxCount i i2 validCount

    % 4.3 pupil data: pupil center (3D position w.r.t. scene camera) and pupil diameter
    pc  = getDataTypeFromTobiiArray(dat,qData(:,1), 'pc',3,{'gidx','eye'},2,qDEBUG);        % pupil center
    pd  = getDataTypeFromTobiiArray(dat,qData(:,2), 'pd',1,{'gidx','eye'},2,qDEBUG);        % pupil diameter
    % 4.4 gaze data: gaze direction vector. gaze position on scene video. 3D gaze position (where eyes intersect, in camera coordinate system?)
    gd  = getDataTypeFromTobiiArray(dat,qData(:,3), 'gd',3,{'gidx','eye'},2,qDEBUG);        % gaze direction
    gp  = getDataTypeFromTobiiArray(dat,qData(:,4), 'gp',2,{'gidx'      },1,qDEBUG);        % gaze position on scene video
    gp3 = getDataTypeFromTobiiArray(dat,qData(:,5),'gp3',3,{'gidx'      },1,qDEBUG);        % gaze convergence position in 3D space
    % 4.5 gyroscope and accelerometer data
    gy  = getDataTypeFromTobiiArray(dat,qData(:,6), 'gy',3,{'ts'        },0,false );        % gyroscope
    ac  = getDataTypeFromTobiiArray(dat,qData(:,7), 'ac',3,{'ts'        },0,false );        % accelerometer
    % 4.6 get video sync info. (e)vts package contains what video
    % timestamp corresponds to a given data timestamp, these occur once
    % in a while so we can see how the two clocks have progressed.
    vts = getDataTypeFromTobiiArray(dat,qData(:,8),'vts',1,{'ts'        },0,false );       % scene video
    if qHasEyeVideo
        evts = getDataTypeFromTobiiArray(dat,qData(:,9),'evts',1,{'ts'        },0,false );       % eye video
    end
    % clean up
    clear dat qData
    % 4.7 parse sync signal (parse directly as it has extra field of interest)
    temp = regexp(syncPortAPISyncTxt,'(?:"ts":)(\d+,).*?(?:"dir":")(in|out)(?:","sig":)(\d})','tokens');   % only get first item of status code, we just care if 0 or not. get delimiter of ts on purpose, so we can concat all matches and just sscanf
    temp = cat(1,temp{:});
    qOut = strcmp(temp(:,2),'out');
    if any(qOut)
        sig.out.ts    = sscanf(cat(2,temp{qOut,1}),'%f,');
        sig.out.state = sscanf(cat(2,temp{qOut,3}),'%f}');
    else
        [sig.out.ts,sig.out.state] = deal([]);
    end
    if any(~qOut)
        sig.in.ts    = sscanf(cat(2,temp{~qOut,1}),'%f,');
        sig.in.state = sscanf(cat(2,temp{~qOut,3}),'%f}');
    else
        [sig.in.ts,sig.in.state] = deal([]);
    end
    
    % 4.8 TODO: parse API sync-packages
    % qAPISync = strcmp(syncPortAPISyncType,'type');
    % clean up
    clear temp syncPortAPISyncTxt syncPortAPISyncType;
    
    
    % 5 reorganize eye data into binocular data, left eye data and right eye data
    data.eye = organizeTobiiGlassesEyeData(pc,pd,gd,gp,gp3);
    clear pc pd gd gp gp3
    data.eye.fs = expectedFs;
    % 5.1 convert gaze vectors to azimuth elevation
    [la,le] = cart2sph(data.eye. left.gd(:,1),data.eye. left.gd(:,3),data.eye. left.gd(:,2));   % matlab's Z and Y are reversed w.r.t. ours
    [ra,re] = cart2sph(data.eye.right.gd(:,1),data.eye.right.gd(:,3),data.eye.right.gd(:,2));
    data.eye. left.azi  =  la*180/pi-90;    % I have checked sign and offset of azi and ele so that things match the gaze position on the scene video in the data file (gp)
    data.eye.right.azi  =  ra*180/pi-90;
    data.eye. left.ele  = -le*180/pi;
    data.eye.right.ele  = -re*180/pi;
    
    % 6 fill up missed samples with nan
    data.eye = fillMissingSamples(data.eye,data.eye.fs);
    
    % 7 add gyroscope and accelerometer data to output file
    assert(issorted(gy.ts,'monotonic'))
    assert(issorted(ac.ts,'monotonic'))
    data.gyroscope      = gy;
    data.accelerometer  = ac;
    clear gy ac
    
    % 8 add video sync data to output file
    assert(issorted( vts.ts,'monotonic'))
    assert(isscalar(unique(vts.ts-vts.vts)))    % this is an assumption of the fts calculation code below
    data.video.scene.sync   =  vts;
    if qHasEyeVideo
        assert(issorted(evts.ts,'monotonic'))
        assert(isscalar(unique(evts.ts-evts.evts)))     % this is an assumption of the fts calculation code below
        data.video.eye.sync     = evts;
    end
    clear vts evts
    
    % 9 add sync port data to output file
    data.syncPort = sig;
    clear sig
    
    % 10 add API sync data to output file (TODO)
    
    % 11 determine t0, convert all timestamps to s
    % set t0 as start point of latest video
    t0s = min(data.video.scene.sync.ts);
    if qHasEyeVideo
        t0s = [t0s min(data.video.eye.sync.ts)];
    end
    t0 = max(t0s);
    data.eye.left.ts        = (data.eye.left.ts-t0)./1000000;
    data.eye.right.ts       = (data.eye.right.ts-t0)./1000000;
    data.eye.binocular.ts   = (data.eye.binocular.ts-t0)./1000000;
    data.gyroscope.ts       = (data.gyroscope.ts-t0)./1000000;
    data.accelerometer.ts   = (data.accelerometer.ts-t0)./1000000;
    data.video.scene.sync.ts= (data.video.scene.sync.ts-t0)./1000000;
    if qHasEyeVideo
        data.video.eye.sync.ts  = (data.video.eye.sync.ts-t0)./1000000;
    end
    data.syncPort.out.ts    = (data.syncPort.out.ts-t0)./1000000;
    data.syncPort. in.ts    = (data.syncPort. in.ts-t0)./1000000;
    
    % 12 open video files for each segment, check how many frames, and make
    % frame timestamps
    data.video.scene.fts        = [];
    data.video.scene.segframes  = [];
    if qHasEyeVideo
        data.video.eye.fts          = [];
        data.video.eye.segframes    = [];
    end
    for s=1:length(segments)
        for p=1:1+qHasEyeVideo
            switch p
                case 1
                    file = 'fullstream.mp4';
                    field= 'scene';
                    tsoff= data.video.scene.sync.ts(data.video.scene.sync.vts==0);
                case 2
                    file = 'eyesstream.mp4';
                    field= 'eye';
                    tsoff= data.video.  eye.sync.ts(data.video.  eye.sync.evts==0);
            end
            fname = fullfile(recordingDir,'segments',segments(s).name,file);
            % get frame timestamps and such from info stored in the mp4
            % file's atoms
            [timeInfo,sttsEntries,atoms,videoTrack] = getMP4VideoInfo(fname);
            % 1. timeInfo (from mdhd atom) contains info about timescale,
            % duration in those units and duration in ms
            % 2. stts table, contains the info needed to determine
            % timestamp for each frame. Use entries in stts to determine
            % frame timestamps. Use formulae described here:
            % https://developer.apple.com/library/content/documentation/QuickTime/QTFF/QTFFChap2/qtff2.html#//apple_ref/doc/uid/TP40000939-CH204-25696
            fIdxs = SmartVec(sttsEntries(:,2),sttsEntries(:,1),'flat');
            timeStamps = cumsum([0 fIdxs]);
            timeStamps = timeStamps/timeInfo.time_scale;
            % last is timestamp for end of last frame, should be equal to
            % length of video
            assert(floor(timeStamps(end)*1000)==timeInfo.duration_ms,'these should match')
            % 3. determine number of frames in file that matlab can read by
            % direct indexing. It seems the Tobii files sometimes have a
            % few frames at the end erroneously marked as keyframes. All
            % those cannot be read by matlab (when using read for a
            % specific time or frame number), so take number of frames as
            % last real keyframe. If not a problem, just take number of
            % frames as last for which we have timeStamp
            lastFrame = atoms.tracks(videoTrack).stss.table(find(diff(atoms.tracks(videoTrack).stss.table)==1,1));
            if isempty(lastFrame)
                lastFrame = sum(sttsEntries(:,1));
            end
            % now that we know number of readable frames, we may have more
            % timestamps than actually readable frames, throw away ones we
            % don't need as we can't read those frames
            assert(length(timeStamps)>=lastFrame)
            timeStamps(lastFrame+1:end) = [];
            % Sync video frames with data by offsetting the timelines for
            % each based on timesync info in tobii data file
            data.video.(field).fts = [data.video.(field).fts timeStamps+tsoff(s)];
            data.video.(field).segframes = [data.video.(field).segframes lastFrame];
            
            % resolution sanity check
            assert(atoms.tracks(videoTrack).tkhd.width ==atoms.tracks(videoTrack).stsd.width , 'mp4 file weird: video widths in tkhd and stsd atoms do not match')
            assert(atoms.tracks(videoTrack).tkhd.height==atoms.tracks(videoTrack).stsd.height,'mp4 file weird: video heights in tkhd and stsd atoms do not match')
            data.video.(field).width  = atoms.tracks(videoTrack).tkhd.width;
            data.video.(field).height = atoms.tracks(videoTrack).tkhd.height;
        end
    end
    % clean up unneeded fields
    for p=1:1+qHasEyeVideo
        switch p
            case 1
                field= 'scene';
            case 2
                field= 'eye';
        end
        data.video.(field) = rmfield(data.video.(field),'sync');
    end
    
    % 13 scale binocular gaze point on video data to pixels
    % we can do so now that we know how big the scene video is
    data.eye.binocular.gp(:,1) = data.eye.binocular.gp(:,1)*data.video.scene.width;
    data.eye.binocular.gp(:,2) = data.eye.binocular.gp(:,2)*data.video.scene.height;
    
    % 14 add time information -- data interval to be used
    % use data from last start of video (scene or eye, whichever is later)
    % to first end of video.
    % 14.1 start time: timestamps are already relative to last video start
    % time, so just get time of first sample at 0 or just before
    data.time.startTime   = data.eye.left.ts(find(data.eye.left.ts<=0,1,'last'));
    % 14.2 end time
    if qHasEyeVideo
        data.time.endTime = min([data.video.scene.fts(end) data.video.eye.fts(end)]);
    else
        data.time.endTime = data.video.scene.fts(end);
    end
    
    % 15 read scene camera calibration info from tslv file
    tslv = readTSLV(fullfile(recordingDir,'segments',segments(s).name,'et.tslv.gz'),'camera',true);
    idxs = find(strcmp(tslv(:,2),'camera'));
    if isempty(idxs)
        % no camera calibration info found, weird but ok
        data.video.scene.calibration = [];
    else
        assert(isscalar(idxs),'more than one camera calibration info found in tslv, contact dcnieho@gmail.com')
        data.video.scene.calibration = rmfield(tslv{idxs,3},'status');  % remove status field, 0==ok is unimportant info for user
    end
    
    % 16 store to cache file
    data.subjName       = participant.pa_info.Name;
    data.recName        = recording.rec_info.Name;
    data.fileVersion    = fileVersion;
    save(fullfile(recordingDir,'livedata.mat'),'-struct','data');
else
    data = load(fullfile(recordingDir,'livedata.mat'));
end

checkMissingFrames(data, 0.05, 0.1);
