function write_vjk(img,fout,res,loc)
%
% function write_vjk(img,fout,res,loc)
%
% Write 3d volume to my made-up format.
% ASCII header will contain 3 lines
%   Nx Ny Nz
%   Dx Dy Dz  : voxel size (mm) (coming from res=[Dx Dy Dz])
%   x0 y0 z0  : location of 1st index (coming from loc=[x0 y0 z0])
%
% The rest is little endian 16-bit unsigned short representing Hounsfield

Nx=size(img,3);
Ny=size(img,2);
Nz=size(img,1);

if(exist('res','var')==0)
  res=[1 1 1];
end

if(exist('loc','var')==0)
  loc=[-(Nx-1)/2*res(1) -(Ny-1)/2*res(2) 0];
end

fp=fopen(fout,'w','l');

fprintf(fp,'%d %d %d\n',Nx,Ny,Nz);
fprintf(fp,'%f %f %f\n',res(1),res(2),res(3));
fprintf(fp,'%f %f %f\n',loc(1),loc(2),loc(3));

for jx = 1:Nx
for jy = 1:Ny
  fwrite(fp,img(:,jy,jx),'uint16');
end
end

fclose(fp);

