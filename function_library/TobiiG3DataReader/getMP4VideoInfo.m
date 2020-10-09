function [time_info,sttsEntries,data,trackIdx] = getMP4VideoInfo(filename)

data.fid = fopen(filename,'rb');

% init
data.eof = 0;
data.total_tracks = 0;

% read first atom header
[headerSize,atomSize,~,typeIdx,data] = read_atom_header(data);

% loop
while ~data.eof
    % MOOV atom - info about contents of mp4 file
    if typeIdx==1
        data.moov.offset = ftell(data.fid) - headerSize;
        data.moov.size = atomSize;
    end
    % MDAT atom - the actual video/audio data
    if typeIdx==130
        data.mdat.offset = ftell(data.fid) - headerSize;
        data.mdat.size = atomSize;
    end
    if typeIdx<128
        % parse subatoms
        data = parse_subatoms(data, atomSize-headerSize);
    else
        % go to next atom header
        status = fseek(data.fid, atomSize-headerSize, 'cof');
        if status~=0
            error 'fseek failed';
        end
    end
    
    % read next atom header
    [headerSize,atomSize,~,typeIdx,data] = read_atom_header(data);
end
fclose(data.fid);
data = rmfield(data,{'fid','eof'});

% get data to return
trackIdx = find(arrayfun(@(x) strcmpi('avc1',x.stsd.type),data.tracks));
time_info = data.tracks(trackIdx).mdhd;
sttsEntries = [data.tracks(trackIdx).stts.sample_count; data.tracks(trackIdx).stts.sample_delta].';
if sttsEntries(end,2)==0
    sttsEntries(end,:) = [];
end
assert(sum(sttsEntries(:,1).*sttsEntries(:,2)) == time_info.duration)
end


%%% helpers
function data = parse_subatoms(data, totalsize)

% init
counted_size = 0;

while counted_size < totalsize
    % read atom header
    [headerSize,atomSize,typeTxt,typeIdx,data] = read_atom_header(data);
    counted_size = counted_size + atomSize;
    
    % track atom
    if typeIdx==2
        data.total_tracks = data.total_tracks + 1;
    end
    
    if typeIdx<128
        % parse subatoms
        data = parse_subatoms(data, atomSize-headerSize);
    else
        % read atom
        data = read_atom(data, atomSize-headerSize, typeTxt);
    end
end
end

function [headerSize,atomSize,typeTxt,typeIdx,data] = read_atom_header(data)

atomSize = read_int(data.fid, 4);
[bytes,count] = fread(data.fid, 4, 'uchar');
if isnan(atomSize) || count<4
    % end of file
    data.eof = 1;
    headerSize = 0;
    atomSize = 0;
    typeTxt = 0;
    typeIdx = 0;
else
    % header and atom size
    headerSize = 8;
    
    % long atom: size coded on 64 bits instead of 32 bits
    if atomSize==1
        headerSize = 16;
        atomSize = read_int(data.fid, 8, [], @int64);
    end
    
    % atom type (text)
    typeTxt = char(bytes)';
    
    % atom type (integer index)
    switch typeTxt
        case 'moov'
            typeIdx = 1;
        case 'trak'
            typeIdx = 2;
        case 'edts'
            typeIdx = 3;
        case 'mdia'
            typeIdx = 4;
        case 'minf'
            typeIdx = 5;
        case 'stbl'
            typeIdx = 6;
        case 'udta'
            typeIdx = 7;
        case 'ilst'
            typeIdx = 8;
        case 'ftyp'
            typeIdx = 129;
        case 'mdat'
            typeIdx = 130;
        otherwise
            typeIdx = 255;
    end
end
end

function data = read_atom(data, size, type)

destination = ftell(data.fid) + size;

switch type
    case 'mdhd'
        version = read_int(data.fid, 4);    % TODO: this should really be 1 8bit for version and 24bit for flags
        if version==1
            fread(data.fid,2,'int64');
            data.tracks(data.total_tracks).mdhd.time_scale = read_int(data.fid, 4);
            data.tracks(data.total_tracks).mdhd.duration = read_int(data.fid, 8);
        else
            fread(data.fid,2,'int32');
            data.tracks(data.total_tracks).mdhd.time_scale = read_int(data.fid, 4);
            data.tracks(data.total_tracks).mdhd.duration = read_int(data.fid, 4);
        end
        data.tracks(data.total_tracks).mdhd.duration_ms = floor(data.tracks(data.total_tracks).mdhd.duration*1000/data.tracks(data.total_tracks).mdhd.time_scale);
        
    case 'tkhd'
        fread(data.fid,1,'uchar');  % version (TODO: below creation time and such can also be 64 bit if version==1)
        fread(data.fid,3,'uchar');  % flags
        data.tracks(data.total_tracks).tkhd.creationTime    = read_int(data.fid, 4);
        data.tracks(data.total_tracks).tkhd.modificationTime= read_int(data.fid, 4);
        data.tracks(data.total_tracks).tkhd.trackId         = read_int(data.fid, 4);
        fread(data.fid,4,'uchar');  % reserved
        data.tracks(data.total_tracks).tkhd.duration        = read_int(data.fid, 4);
        fread(data.fid,8,'uchar');  % reserved
        fread(data.fid,2,'uchar');  % layer
        fread(data.fid,2,'uchar');  % alternative group
        fread(data.fid,2,'uchar');  % volume
        fread(data.fid,2,'uchar');  % reserved
        fread(data.fid,9*4,'uchar');  % matrix
        data.tracks(data.total_tracks).tkhd.width           = read_int(data.fid, 4)/2^16;
        data.tracks(data.total_tracks).tkhd.height          = read_int(data.fid, 4)/2^16;
        
        
    case 'stsd'
        fread(data.fid,1,'uchar');  % version
        fread(data.fid,3,'uchar');  % flags
        data.tracks(data.total_tracks).stsd.entry_count = read_int(data.fid, 4);
        % there's a box in this box
        [~,~,typeTxt,~,data] = read_atom_header(data);
        data.tracks(data.total_tracks).stsd.type = typeTxt;
        if strcmp(typeTxt,'avc1')
            % its an avc1 box, read some info about the encoded video data
            fread(data.fid,6,'uchar');  % reserved
            fread(data.fid,2,'uchar');  % data reference index
            fread(data.fid,2,'uchar');  % video encoding version
            fread(data.fid,2,'uchar');  % video encoding revision level
            fread(data.fid,4,'uchar');  % video encoding vendor
            fread(data.fid,4,'uchar');  % video temporal quality
            fread(data.fid,4,'uchar');  % video spatial quality
            data.tracks(data.total_tracks).stsd.width           = read_int(data.fid, 2);
            data.tracks(data.total_tracks).stsd.height          = read_int(data.fid, 2);
            data.tracks(data.total_tracks).stsd.dpiHori         = read_int(data.fid, 4)/2^16;
            data.tracks(data.total_tracks).stsd.dpiVert         = read_int(data.fid, 4)/2^16;
        end
        
    case 'stts'
        fread(data.fid,1,'uchar');
        fread(data.fid,3,'uchar');
        stts_entry_count = read_int(data.fid, 4);
        temp = read_int(data.fid, 4, stts_entry_count*2);
        data.tracks(data.total_tracks).stts.sample_count = temp(1:2:end);
        data.tracks(data.total_tracks).stts.sample_delta = temp(2:2:end);
        
    case 'stss'
        fread(data.fid,1,'uchar');
        fread(data.fid,3,'uchar');
        stss_entry_count = read_int(data.fid, 4);
        data.tracks(data.total_tracks).stss.table = read_int(data.fid, 4, stss_entry_count);
        
    case 'stsc'
        fread(data.fid,1,'uchar');
        fread(data.fid,3,'uchar');
        stsc_entry_count = read_int(data.fid, 4);
        temp = read_int(data.fid, 4, stsc_entry_count*3);
        data.tracks(data.total_tracks).stsc.first_chunk = temp(1:3:end);
        data.tracks(data.total_tracks).stsc.samples_per_chunk = temp(2:3:end);
        data.tracks(data.total_tracks).stsc.sample_desc_index = temp(3:3:end);
        
    case 'stsz'
        fread(data.fid,1,'uchar');
        fread(data.fid,3,'uchar');
        data.tracks(data.total_tracks).stsz.sample_size = read_int(data.fid, 4);
        data.tracks(data.total_tracks).stsz.sample_count = read_int(data.fid, 4);
        if data.tracks(data.total_tracks).stsz.sample_size==0
            data.tracks(data.total_tracks).stsz.table = read_int(data.fid, 4, data.tracks(data.total_tracks).stsz.sample_count);
        else
            data.tracks(data.total_tracks).stsz.table = data.tracks(data.total_tracks).stsz.sample_size*ones(1,data.tracks(data.total_tracks).stsz.sample_count);
        end
        
    case 'stco'
        fread(data.fid,1,'uchar');
        fread(data.fid,3,'uchar');
        stco_entry_count = read_int(data.fid, 4);
        data.tracks(data.total_tracks).stco.chunk_offset = read_int(data.fid, 4, stco_entry_count);
end

current_position = ftell(data.fid);
if destination ~= current_position
    status = fseek(data.fid, destination, 'bof');
    if status~=0
        error 'fseek failed';
    end
end
end

function out = read_int(fid, n, m, classfun)

if nargin<3 || isempty(m)
    m = 1;
end
if nargin<4
    classfun = @double;
end
[temp,count] = fread(fid, n*m, 'uchar');
if count<n*m
    out = nan;
    return
end
temp = classfun(reshape(temp, n, m));
p = repmat(classfun(2).^classfun(8*(n-1:-1:0))', 1, m);
out = sum(temp .* p, 'native');
end
