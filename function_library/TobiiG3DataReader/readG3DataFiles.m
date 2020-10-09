function data = readG3DataFiles(recordingDir,qDEBUG)

% Cite as: Niehorster, D.C., Hessels, R.S., and Benjamins, J.S. (2020).
% GlassesViewer: Open-source software for viewing and analyzing data from
% the Tobii Pro Glasses 2 eye tracker. Behavior Research Methods. doi:
% 10.3758/s13428-019-01314-1

% set file format version. cache files older than this are overwritten with
% a newly generated cache file
fileVersion = 1;

if ~isempty(which('matlab.internal.webservices.fromJSON'))
    jsondecoder = @matlab.internal.webservices.fromJSON;
elseif ~isempty(which('jsondecode'))
    jsondecoder = @jsondecode;
else
    error('Your MATLAB version does not provide a way to decode json (which means its really old), upgrade to something newer');
end

qGenCacheFile = ~exist(fullfile(recordingDir,'gazedata.mat'),'file');
if ~qGenCacheFile
    % we have a cache file, check its file version
    cache = load(fullfile(recordingDir,'gazedata.mat'),'fileVersion');
    qGenCacheFile = cache.fileVersion~=fileVersion;
end

if qGenCacheFile
    % 0 get info about participant and recording
    fid = fopen(fullfile(recordingDir,'recording.g3'),'rt');
    recording = jsondecoder(fread(fid,inf,'*char').');
    fclose(fid);
    expectedFs = round(recording.gaze.samples/recording.duration/50)*50;    % find nearest 50Hz
    if qDEBUG
        fprintf('determined fs: %d Hz\n',expectedFs);
    end
    fid = fopen(fullfile(recordingDir,recording.meta_folder,'participant'),'rt');
    participant = jsondecoder(fread(fid,inf,'*char').');
    fclose(fid);
    
    % 1 read in gaze data
    % 1.1 unpack gz file, if doesn't exist
    gzFile = fullfile(recordingDir,recording.gaze.file);
    [~,gazeFile,~] = fileparts(gzFile);
    if ~exist(fullfile(recordingDir,gazeFile),'file')
        gunzip(fullfile(recordingDir,recording.gaze.file));
    end
    % 1.2 read in gaze data
    fid = fopen(fullfile(recordingDir,gazeFile),'rt');
    gazeData = fread(fid,inf,'*char').';
    fclose(fid);
    % turn into something we can read
    gazeData(gazeData==10) = ',';
    gazeData = jsondecoder(['[' gazeData ']']);
    % do quick checks
    types = unique({gazeData.type});
    assert(isscalar(types) && strcmp(types{1},'gaze'),'Data not as expected')
    
    % 2 turn into our data format
    % 2.1 prep storage
    data.eye.fs             = expectedFs;
    data.eye.left.ts        = cat(1,gazeData.timestamp);
    data.eye.right.ts       = data.eye.left.ts;
    data.eye.binocular.ts   = data.eye.left.ts;
    nSamp                   = length(gazeData);
    [data.eye.left.pc, data.eye.left.gd]    = deal(nan(nSamp,3));
    data.eye.left.pd                        =      nan(nSamp,1) ;
    [data.eye.right.pc, data.eye.right.gd]  = deal(nan(nSamp,3));
    data.eye.right.pd                       =      nan(nSamp,1) ;
    data.eye.binocular.gp                   =      nan(nSamp,2) ;
    data.eye.binocular.gp3                  =      nan(nSamp,3) ;
    % 2.2 throw data into storage
    qNotMissing                             = arrayfun(@(x) isfield(x.data,'gaze2d'),gazeData); % struct is empty (and thus doesn't have this field) if there is no data
    gazeData                                = cat(1,gazeData(qNotMissing).data);
    data.eye.binocular.gp (qNotMissing,:)   = cat(2,gazeData.gaze2d).';
    data.eye.binocular.gp3(qNotMissing,:)   = cat(2,gazeData.gaze3d).';
    [qNotMissingLA,qNotMissingRA]           = deal(qNotMissing);
    qNotMissingL                            = arrayfun(@(x) isfield(x.eyeleft ,'gazeorigin'),gazeData);
    qNotMissingR                            = arrayfun(@(x) isfield(x.eyeright,'gazeorigin'),gazeData);
    qNotMissingLA(qNotMissing)              = qNotMissingL;
    qNotMissingRA(qNotMissing)              = qNotMissingR;
    left                                    = cat(1,gazeData(qNotMissingL).eyeleft);
    right                                   = cat(1,gazeData(qNotMissingR).eyeright);
    data.eye.left .pc(qNotMissingLA,:)      = cat(2,left.gazeorigin).';
    data.eye.left .pd(qNotMissingLA)        = cat(2,left.pupildiameter).';
    data.eye.left .gd(qNotMissingLA,:)      = cat(2,left.gazedirection).';
    data.eye.right.pc(qNotMissingRA,:)      = cat(2,right.gazeorigin).';
    data.eye.right.pd(qNotMissingRA)        = cat(2,right.pupildiameter).';
    data.eye.right.gd(qNotMissingRA,:)      = cat(2,right.gazedirection).';
    % 2.3 for each binocular sample, see on how many eyes its based
    data.eye.binocular.nEye                 = sum([qNotMissingLA qNotMissingRA],2);
    % clean up
    clear gazeData left right qNotMissing qNotMissingL qNotMissingLA qNotMissingR qNotMissingRA
    % 2.4 convert gaze vectors to azimuth elevation
    [la,le] = cart2sph(data.eye. left.gd(:,1),data.eye. left.gd(:,3),data.eye. left.gd(:,2));   % matlab's Z and Y are reversed w.r.t. ours
    [ra,re] = cart2sph(data.eye.right.gd(:,1),data.eye.right.gd(:,3),data.eye.right.gd(:,2));
    data.eye. left.azi  =  la*180/pi-90;    % I have checked sign and offset of azi and ele so that things match the gaze position on the scene video in the data file (gp)
    data.eye.right.azi  =  ra*180/pi-90;
    data.eye. left.ele  = -le*180/pi;
    data.eye.right.ele  = -re*180/pi;
    % clean up
    clear la le ra re
    
    % 3 read in event data
    % 3.1 unpack gz file, if doesn't exist
    gzFile = fullfile(recordingDir,recording.events.file);
    [~,eventFile,~] = fileparts(gzFile);
    if ~exist(fullfile(recordingDir,eventFile),'file')
        gunzip(fullfile(recordingDir,recording.events.file));
    end
    % 3.2 read in event data
    fid = fopen(fullfile(recordingDir,eventFile),'rt');
    eventData = fread(fid,inf,'*char').';
    fclose(fid);
    % turn into something we can read
    eventData(eventData==10) = ',';
    eventData = jsondecoder(['[' eventData ']']);
    % 3.3 sync signal
    qSync= strcmp({eventData.type},'syncport');
    sync = cat(1,eventData(qSync).data);
    ts   = cat(1,eventData(qSync).timestamp);
    qOut = strcmp({sync.direction},'out');
    if any(qOut)
        data.syncPort.out.ts    = ts(qOut);
        data.syncPort.out.state = cat(1,sync(qOut).value);
    else
        [data.syncPort.out.ts,data.syncPort.out.state] = deal([]);
    end
    if any(~qOut)
        data.syncPort.in.ts    = ts(~qOut);
        data.syncPort.in.state = cat(1,sync(~qOut).value);
    else
        [data.syncPort.in.ts,data.syncPort.in.state] = deal([]);
    end
    % clean up
    clear eventData sync ts qSync qOut
    
    % 4 open scene video file, check how many frames, and make frame
    % timestamps
    data.video.scene.fts        = [];
    data.video.scene.segframes  = [];
    if false
        data.video.eye.fts          = [];
        data.video.eye.segframes    = [];
    end
    for p=1:1+false
        switch p
            case 1
                file = recording.scenecamera.file;
                field= 'scene';
                %tsoff= data.video.scene.sync.ts(data.video.scene.sync.vts==0);
            case 2
                error todo
                file = 'eyesstream.mp4';
                field= 'eye';
                tsoff= data.video.  eye.sync.ts(data.video.  eye.sync.evts==0);
        end
        fname = fullfile(recordingDir,file);
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
        data.video.(field).fts = [data.video.(field).fts timeStamps];
        data.video.(field).segframes = [data.video.(field).segframes lastFrame];
        
        % resolution sanity check (doesn't work, tkhd width and height
        % appear to be broken)
%         assert(atoms.tracks(videoTrack).tkhd.width ==atoms.tracks(videoTrack).stsd.width , 'mp4 file weird: video widths in tkhd and stsd atoms do not match')
%         assert(atoms.tracks(videoTrack).tkhd.height==atoms.tracks(videoTrack).stsd.height,'mp4 file weird: video heights in tkhd and stsd atoms do not match')
        data.video.(field).width  = atoms.tracks(videoTrack).stsd.width;
        data.video.(field).height = atoms.tracks(videoTrack).stsd.height;
        
        % store name of video file
        data.video.(field).file   = file;
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
    if isempty(data.time.startTime)
        data.time.startTime = data.eye.left.ts(1);
    end
    % 14.2 end time
    if false
        data.time.endTime = min([data.video.scene.fts(end) data.video.eye.fts(end)]);
    else
        data.time.endTime = data.video.scene.fts(end);
    end
    
    % 15 read scene camera calibration info from tslv file
    data.video.scene.calibration = recording.scenecamera.camera_calibration;
    
    % 16 store to cache file
    data.subjName       = participant.name;
    data.fileVersion    = fileVersion;
    save(fullfile(recordingDir,'gazedata.mat'),'-struct','data');
else
    data = load(fullfile(recordingDir,'gazedata.mat'));
end

checkMissingFrames(data, 0.05, 0.1);
