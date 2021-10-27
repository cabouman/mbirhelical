function sino=read_sino_float(filename)
%
% function sino=read_sino_float(filename)
%
% Reads helical sinogram into sino(Nr,Nc,Nv).
%

fp = fopen(filename,'r');
Nr = fscanf(fp,'%d',1);
Nc = fscanf(fp,'%d',1);
Nv = fscanf(fp,'%d',1);
temp = fread(fp,1,'int8');  % carriage return

% 3D memory allocation
sino = zeros(Nr,Nc,Nv);
for jv = 1:Nv
  sino(:,:,jv) = fread(fp,[Nr Nc],'float32');
end

fclose(fp);

