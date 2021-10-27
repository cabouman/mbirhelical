function [img res loc]=read_vjk(fin)
%
% function [img res loc]=read_vjk(fin)
%
% Read 3d volume using my made-up format.
% ASCII header should contain 3 lines
%   Nx Ny Nz
%   Dx Dy Dz  : voxel size (mm) (coming from res=[Dx Dy Dz])
%   x0 y0 z0  : location of 1st index (coming from loc=[x0 y0 z0])
%
% The rest is little endian 16-bit unsigned short representing Hounsfield
%
% Input
%   fin : filename
% Outputs
%   img : 3d image array in order (z,y,x)
%   res : voxel size [Dx Dy Dz]
%   loc : coords of 1st voxel [x0 y0 z0]

fp=fopen(fin,'r','l');

Nx=fscanf(fp,'%d',1);
Ny=fscanf(fp,'%d',1);
Nz=fscanf(fp,'%d',1);
temp=fread(fp,1,'int8');  % carriage return
Dx=fscanf(fp,'%f',1);
Dy=fscanf(fp,'%f',1);
Dz=fscanf(fp,'%f',1);
temp=fread(fp,1,'int8');  % carriage return
x0=fscanf(fp,'%f',1);
y0=fscanf(fp,'%f',1);
z0=fscanf(fp,'%f',1);
temp=fread(fp,1,'int8');  % carriage return

res=[Dx Dy Dz];
loc=[x0 y0 z0];

img=fread(fp,inf,'uint16');
%img(1:2)
%length(img)
%Nx*Ny*Nz

if(length(img) ~= Nx*Ny*Nz)
  disp('read_3dimg error : data length does not match header')
  img=0;
else
  img=reshape(img,Nz,Ny,Nx);
end

fclose(fp);

