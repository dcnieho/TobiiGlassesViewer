function x = getDataTypeFromTobiiArray(dat,qType,type,nCol,sortFields,mode,qDEBUG)
% mode:
% 0: no gidx
% 1: has gidx, one value per gidx/timstamp expected
% 2: has gidx, two values per gidx/timstamp expected (one per eye)

% find data type we want in the array, and get those
x = struct('ts',{dat.ts(qType)},type,{dat.dat(qType,1:nCol)});
if mode==1 || mode==2
    x.gidx = dat.gidx(qType);
end
if mode==2
    x.eye = dat.eye(qType);
end

% sort by specified columns
x = sortByColumns(x,sortFields);

if mode>0
    % remove data where time apparently went backward
    iNeg = find(diff(x.ts)<0);
    iNeg(diff(x.gidx([iNeg iNeg+1]))==0) = []; % if negative timestep happened within same gidx, probably just slight time difference for different eyes. Don't worry about it
    x = replaceElementsInStruct(x,iNeg+1,[],[],1);
else
    % remove invalid data
    x = replaceElementsInStruct(x,~dat.qValid(qType),[],[],1);
end

% check we have expected number of values per gidx
if mode>0
    [~,n] = UniqueAndN(x.gidx);
    if qDEBUG
        if mode==1
            fprintf('% 3s: 1/gidx: % 5d, 2/gidx: % 5d\n',type,sum(n==1),sum(n==2));
        else
            fprintf('% 3s: 1/gidx: % 5d, 2/gidx: % 5d, 3/gidx: % 5d\n',type,sum(n==1),sum(n==2),sum(n==3));
        end
        
    end
end
end
