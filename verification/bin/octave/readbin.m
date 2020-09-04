function [dat, n, m] = readbin(fn)

    fid = fopen(fn);
    if(~fid)
        error('cannot to open file');
    end

    type = fread(fid, 1, 'int32');
    n    = fread(fid, 1, 'int32');
    m    = fread(fid, 1, 'int32');

    if(type==0)
        dat = fread(fid, [n*m, 1], 'double');
    end

    if(type==1)
        y = fread(fid, [n*m*2, 1], 'double');
        dat = y(1:2:end) + 1i * y(2:2:end);
    end

    dat = reshape(dat, n, m);

    fclose(fid);
end