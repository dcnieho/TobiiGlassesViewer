function parsedData = readTSLV(filename,typeList,exitAfterFirst)
% NB: This reads TSLV files until an "End" package is found. According to
% docs, such an End package may be missing e.g. if the device is suddenly
% switched off. This here will then crash upon trying to read after the end
% of the file. If you have such files and need this to work, contact the
% developer.

% if gz file passed and filename and unpacked version doesn't exist yet, unpack it
if strcmp(filename(max(1,end-2):end),'.gz')
    if ~exist(filename(1:max(1,end-3)),'file')
        gunzip(filename);
    end
    filename = filename(1:max(1,end-3));
end

% parse input
if nargin>1 && ~isempty(typeList)
    if ~iscell(typeList)
        typeList = {typeList};
    end
    typeList = cellfun(@typeIDStringConversion,typeList,'uni',false);
    assert(~any(cellfun(@isempty,typeList)),'One of the types you specified in typeList is not recognized')
    typeList = [typeList{:}];
else
    typeList = [];
end

% init
parsedData  = cell(1024,3);
i           = 0;

% open tslv
fid = fopen(filename,'rb');
% loop
endType = typeIDStringConversion('end');
while true
    % 1. read header of data package
    [type,status,payloadLength,payloadLengthPadded] = readHeader(fid);
    
    if isempty(typeList) || ismember(type,typeList)
        % 2. read data package
        out = readDataPackage(fid,type,status,payloadLength);
        if ~isempty(out)
            i=i+1;
            if i>size(parsedData,1)
                % grow to double the size
                parsedData(size(parsedData,1)*2,3) = {[]};
            end
            parsedData{i,1} = out.typeID;
            parsedData{i,2} = out.type;
            parsedData{i,3} = out.payload;
            if nargin>2 && exitAfterFirst
                break;
            end
        end
        
        % 3. go to next data package
        skipBytes(fid,payloadLengthPadded-payloadLength);
    else
        % 2. skip whole package
        skipBytes(fid,payloadLengthPadded);
    end
    % 4. check we're at end
    if type==endType
        break;
    end
end
fclose(fid);

% trim off excess from output cell
parsedData(i+1:end,:) = [];


end


%%% helpers
function [type,status,payloadLength,payloadLengthPadded] = readHeader(fid)

type   = readInt(fid, 2);
status = readInt(fid, 2);
payloadLength = readInt(fid, 4, []);
payloadLengthPadded = ceil(payloadLength/4)*4; % payloads are padded to multiples of 4 bytes, calculate how much we have to skip later
end


function out = readDataPackage(fid, type, status, payloadLength)

outputBuilder = @(varargin) buildOutput(type,status,typeIDStringConversion(type),varargin{:});
switch type
    case 0
        % unknown, skip
        skipBytes(fid,payloadLength);
        out = [];
        
    case 21332
        % FORMAT ID: depicts the stream data format. Always "TSLV1"
        format = char(readInt(fid, 1, payloadLength));
        assert(strcmp(format,'TSLV1'))
        out = outputBuilder('TSLVformat',format);
        
    case 2
        % STREAM TYPE: data stream content type: "ET" or "Mems".
        streamType = char(readInt(fid, 1, payloadLength));
        assert(ismember(streamType,{'ET','MEMS'}))
        out = outputBuilder('streamType',streamType);
        
    case 3
        % FREQUENCY: indicative (but not guaranteed or strict) data point frequency.
        frequency = readInt(fid, 4);
        out = outputBuilder('frequency',frequency);
        
    case 4
        % VIDEO_INFO: Video information for the video.
        videoId = readInt(fid, 1);
        skipBytes(fid,1);   % padding
        width = readInt(fid, 2);
        height = readInt(fid, 2);
        frequency = readInt(fid, 2);
        name = char(readInt(fid, 1, payloadLength-8));
        name(name==0) = [];
        out = outputBuilder('videoId',videoId,'width',width,'height',height,'frequency',frequency,'name',name);
        
    case 5
        % END: Data stream end marker. Will be the last TSLV of the data stream.
        %      See status bits and data value for end reason.
        %      If a stream is missing the trailing END TSLV it means the stream
        %      is truncated for whatever reason (out of storage space, device shutdown,
        %      etc).
        out = outputBuilder();
        
    case 6
        % CHECKSUM: Rolling checksum of data stream.
        %           The CHECKSUM TSLV itself is not included in the stream.
        %           Each CHECKSUM resets the checksum counting.
        %           Checksums are SHA1.
        % not implemented, skip
        skipBytes(fid,payloadLength);
        out = [];
        
    case 10
        % FRAME_ID: Marks the beginning of a new frame. The frame extends to
        %           the next FRAME_ID (non-inclusive).
        %           The 'seq' is an ever-increasing (by one) frame sequence id.
        %           'sid' is the internal session id and should be ignored.
        frameID   = readInt(fid, 8, [], @int64);
        sessionID = readInt(fid, 4);                      %#ok<NASGU>
        out = outputBuilder('frameID',frameID);
        
    case 11
        % SESSION_START: Indicates that a new internal session is starting.
        %                This functions as an internal state synchronization and should be
        %                ignored by any external readers.
        %                The session id is in the FRAME_ID TSLV.
        % not implemented, skip
        skipBytes(fid,payloadLength);
        out = [];
        
    case 12
        % SESSION_STOP: Indicates that the current session is ending.
        %               This functions as an internal state synchronization and should be
        %               ignored by any external readers.
        %               The session id is in the FRAME_ID TSLV.
        % not implemented, skip
        skipBytes(fid,payloadLength);
        out = [];
        
    case 50
        % SYSTEM_TIMESTAMP: Monotonic system timestamp in microseconds (64-bits).
        %                   Included in every frame.
        timestamp   = readInt(fid, 8, [], @int64);
        out = outputBuilder('timestamp',timestamp);
        
    case 51
        % WALLCLOCK_TIMESTAMP: The current system wall clock time.
        %                      Also indicates whether the time is synced by
        %                      external time source or not.
        %                      Included every 10 seconds.
        timestamp   = readInt(fid, 8, [], @int64);
        formattedTime = char(readInt(fid, 1, 8));
        formattedTime(formattedTime==0) = [];
        out = outputBuilder('timestamp',timestamp,'formattedTime',formattedTime);
        
    case 53
        % VIDEOCLOCK_TIMESTAMP: The current video clock timestamp.
        videoId = readInt(fid, 1);
        skipBytes(fid,3);   % padding
        timestamp   = readInt(fid, 8, [], @int64);
        out = outputBuilder('videoId',videoId,'timestamp',timestamp);
        
    case 54
        % VIDEOFILE_TIMESTAMP: The current video file timestamp.
        videoId = readInt(fid, 1);
        skipBytes(fid,3);   % padding
        timestamp   = readInt(fid, 8, [], @int64);
        out = outputBuilder('videoId',videoId,'timestamp',timestamp);
        
    case 57
        % GazeFrameCounter
        frameCounter   = readInt(fid, 8, [], @int64);
        out = outputBuilder('counter',frameCounter);
        
    case 103
        % PUPIL_CENTER: Position in glasses 3D coordinate frame.
        %               (Origin at scene camera center of projection, z-axis pointing forward,
        %               y-axis pointing down?).
        eye = eyeIDToEye(readInt(fid, 1));
        skipBytes(fid,3);   % padding
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        out = outputBuilder('eye',eye,'pos',[x y z]);
        
    case 104
        % GAZE_DIRECTION: Direction in glasses 3D coordinate frame (unit vector).
        eye = eyeIDToEye(readInt(fid, 1));
        skipBytes(fid,3);   % padding
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        out = outputBuilder('eye',eye,'vec',[x y z]);
        
    case 105
        % pupil diameter
        eye = eyeIDToEye(readInt(fid, 1));
        skipBytes(fid,3);   % padding
        diam = readFloat(fid, 'single');
        out = outputBuilder('eye',eye,'diam',diam);
        
    case 110
        % 3D GAZE POINT: Position in glasses 3D coordinate frame.
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        out = outputBuilder('vec',[x y z]);
        
    case 112
        % 2D GAZE POINT: Normalized x,y axis (0-1,0-1) according to the  full stream aspect ratio.
        %                May be outside 0-1 if point falls outside scene camera view scope.
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        out = outputBuilder('vec',[x y]);
        
    case 200
        % Gyro data
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        out = outputBuilder('vec',[x y z]);
        
    case 201
        % Accelerometer data
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        out = outputBuilder('vec',[x y z]);
        
    case 250
        % logged event
        % not implemented, skip
        warning('todo: implement');
        skipBytes(fid,payloadLength);
        out = [];
        
    case 251
        % sync signal
        % not implemented, skip
        signal = readInt(fid, 1);
        direction = directionToString(readInt(fid, 1));
        out = outputBuilder('signal',signal,'direction',direction);
        
    case 300
        % CAMERA: Camera positioning, manufacturing calibration, etc.
        %         Position and rotation in glasses 3D coordinate frame.
        id = readInt(fid, 1);
        location = readInt(fid, 1);
        
        skipBytes(fid,2);   % padding
        
        x = readFloat(fid, 'single');
        y = readFloat(fid, 'single');
        z = readFloat(fid, 'single');
        
        r11 = readFloat(fid, 'single');
        r12 = readFloat(fid, 'single');
        r13 = readFloat(fid, 'single');
        
        % the camera parameters contain a 3x3 matrix for the rotation. we have to skip the next 6 floats since we don't use them.
        readFloat(fid, 'single'); % r21
        readFloat(fid, 'single'); % r22
        readFloat(fid, 'single'); % r23
        readFloat(fid, 'single'); % r31
        readFloat(fid, 'single'); % r32
        readFloat(fid, 'single'); % r33
        
        fx = readFloat(fid, 'single');
        fy = readFloat(fid, 'single');
        
        skew = readFloat(fid, 'single');
        
        px = readFloat(fid, 'single');
        py = readFloat(fid, 'single');
        
        rd1 = readFloat(fid, 'single');
        rd2 = readFloat(fid, 'single');
        rd3 = readFloat(fid, 'single');
        
        t1 = readFloat(fid, 'single');
        t2 = readFloat(fid, 'single');
        t3 = readFloat(fid, 'single');
        sx = readInt(fid, 2);
        sy = readInt(fid, 2);
        out = outputBuilder('id',id,'location',location,'position',[x y z],'rodriguesRotation',[r11 r12 r13],'focalLength',[fx fy],'skew',skew,'principalPoint',[px py],'radialDistortion',[rd1 rd2 rd3],'tangentialDistortion',[t1 t2 t3],'sensorDimensions',[sx sy]);
        
    otherwise
        % not implemented, skip
        warning('todo: implement');
        skipBytes(fid,payloadLength);
        out = [];
end
end

function out = buildOutput(typeID,status,typeString,varargin)
out.typeID  = typeID;
out.type    = typeString;
out.payload = struct(varargin{:});
if ~isempty(varargin)
    out.payload.status = status;
end
end

function eye = eyeIDToEye(eyeID)
if eyeID==0
    eye = 'L';
else
    eye = 'R';
end
end

function eye = directionToString(direction)
if direction==0
    eye = 'in';
else
    eye = 'out';
end
end

function out = readInt(fid, nBit, nElem, classfun)
if nargin<3 || isempty(nElem)
    nElem = 1;
end
if nargin<4
    classfun = @double;
end
[temp,count] = fread(fid, nBit*nElem, '*uchar');
if count<nBit*nElem
    out = nan;
    return
end
temp = classfun(reshape(temp, nBit, nElem));
p = repmat(classfun(2).^classfun(8*(0:nBit-1))', 1, nElem);
out = sum(temp .* p, 1, 'native');
end

function out = readFloat(fid, type, nElem)
if nargin<3 || isempty(nElem)
    nElem = 1;
end
[temp,count] = fread(fid, nElem, type);
if count<nElem
    out = nan;
    return
end
out = double(temp);
end

function skipBytes(fid,nSkip)
if nSkip~=0
    fread(fid, nSkip, '*uchar');
end
end

function typeIDOrString = typeIDStringConversion(typeIDOrString)
table = {
    21332,'format'
    2,'streamType'
    3,'frequency'
    4,'videoInfo'
    5,'end'
    10,'frameID'
    50,'systemTimestamp'
    51,'wallClockTimestamp'
    53,'videoClockTimestamp'
    54,'videoFileTimestamp'
    57,'gazePackageCounter'
    103,'pupilCenter'
    104,'gazeDirection'
    105,'pupilDiameter'
    110,'gazePoint3D'
    112,'gazePoint2D'
    200,'gyro'
    201,'accelerometer'
    251,'syncSignal'
    300,'camera'
    };

if ischar(typeIDOrString)
    idx = find(strcmp(table(:,2),typeIDOrString),1);
    if ~isempty(idx)
        typeIDOrString = table{idx,1};
    else
        typeIDOrString = [];
    end
else
    idx = find(cat(1,table{:,1})==typeIDOrString,1);
    if ~isempty(idx)
        typeIDOrString = table{idx,2};
    else
        typeIDOrString = '';
    end
end
end
