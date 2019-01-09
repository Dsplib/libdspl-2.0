function [dat, n, m] = dspl_readbin(fn)

fid = fopen(fn);
if(~fid)
	error('cannot to open file');
end
type = fread(fid, 1, 'int32');	
n    = fread(fid, 1, 'int32');
m    = fread(fid, 1, 'int32');

if(type==0)
	dat = fread(fid, [n*m, 1], 'double');
	dat = reshape(dat, n, m);
end
if(type==1)
	dat = fread(fid, [n*m*2, 1], 'double');
	dat = dat(1:2:end) + 1i * dat(2:2:end);
	dat = reshape(dat, n, m);	
end

fclose(fid);

end
