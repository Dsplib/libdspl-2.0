function res = writebin(x, type, fn)

    if(type~=0 && type~=1)
        res = 2;
        return;
    end

    fid = fopen(fn, 'w');
    if(~fid)
        error('cannot to open file');
    end

    n = size(x, 1);
    m = size(x, 2);
        
    count = fwrite(fid,  type,  'int32');
    if(count ~= 1)
        res = 1;
        return;
    end
    count = fwrite(fid,  n,     'int32');
    if(count ~= 1)
        fclose(fid);
        res = 1;
        return;
    end
    count = fwrite(fid,  m,     'int32');
    if(count ~= 1)
        res = 1;
        fclose(fid);
        return;
    end
    
    flag = 0;
    if(type==0)
        
        count = fwrite(fid,  x,     'double');
        if(count ~= n*m)
            res = 1;
            fclose(fid);
            return;
        end
        flag = 1;
        
    else
        y = reshape(x, n*m, 1);
        z = zeros(2*n*m, 1);
        z(1:2:end) = real(y);
        z(2:2:end) = imag(y);
        count = fwrite(fid,  z,     'double');
        if(count ~= 2*n*m)
            res = 1;
            fclose(fid);
            return;
        end
        flag = 1;
    end
    if(flag == 0)
        res = 3;
    else 
        res = 0;
    end
    fclose(fid);
end