function dspl_writebin(dat, fn)

fid = fopen(fn, 'w');
if(~fid)
	error('cannot to open file');
end
if(isreal (dat))
  type = 0;
else
  type = 1;
end

n    = size(dat,1);
m    = size(dat,2);

fwrite (fid, type, 'int32');
fwrite (fid, n, 'int32');
fwrite (fid, m, 'int32');

if(type==0)
	dat = reshape(dat, n*m, 1);
  fwrite (fid, dat, 'double');
end
if(type==1)
  y = zeros(2*n*m, 1);
  dat = reshape(dat, n*m, 1);
  y(1:2:end) = real(dat);
  y(2:2:end) = imag(dat);
 
	fwrite (fid, y, 'double');	
end

fclose(fid);

end
