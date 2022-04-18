function [img]=read_vjk_float(fin)
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
 
Nx=fread(fp,1,'int16');
temp=fread(fp,1,'int16');  % carriage return
Ny=fread(fp,1,'int16');
temp=fread(fp,1,'int16');  % carriage return
Nz=fread(fp,1,'int16');
temp=fread(fp,1,'int16');  % carriage return
 
img=fread(fp,inf,'float32');
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
