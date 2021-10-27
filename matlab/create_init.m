function create_init(fname, init, Nx, Ny, Nz)

img = init*ones(Nz,Ny,Nx);

xc = 0;
yc = 0;
zc = Nz/2;
Del_xy = 2;
Del_z = 2;

res=[Del_xy Del_xy Del_z];
loc=[ xc-(N-1)*Del_xy/2  yc-(N-1)*Del_xy/2  zc-(N-1)*Del_z/2 ];
img=1000*img/(.02-.0000226);
write_vjk(img,fname,res,loc);

return

% Output phantom as .img binary file 
fp = fopen(fname, 'w');
fprintf(fp, '%d %d %d \n', Nx, Ny, Nz);
fprintf(fp,'%f %f %f \n', xc, yc, zc);
fprintf(fp,'%f %f \n', Del_xy, Del_z);
fprintf(fp,'%f \n', rI);
for z = 1:Nz
	for y = 1:Ny
		for x = 1:Nx
			fprintf(fp,'%f ', img(z,y,x));
		end
	end
end
fclose(fp);
