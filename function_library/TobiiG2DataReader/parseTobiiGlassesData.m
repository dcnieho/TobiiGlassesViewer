function dat = parseTobiiGlassesData(txt,name,type)
switch type
    case 1
        regFmt = ['(?<="' name '":)[-\d.]+,'];
        sscFmt = '%f,';
        numElem= 1;
    case 2
        regFmt = ['(?<="' name '":)[-\d.]+}'];
        sscFmt = '%f}';
        numElem= 1;
    case 3
        regFmt = ['(?<="' name '":)\[(?:[-\d.]+,{0,1})+\]'];
        sscFmt = '[%f,%f]';
        numElem= 2;
    case 4
        regFmt = ['(?<="' name '":)\[(?:[-\d.]+,{0,1})+\]'];
        sscFmt = '[%f,%f,%f]';
        numElem= 3;
end
dat = regexp(txt,regFmt,'match');
if numElem>1
    dat = reshape(sscanf(cat(2,dat{:}),sscFmt),numElem,[]).';
else
    dat =         sscanf(cat(2,dat{:}),sscFmt);
end
end