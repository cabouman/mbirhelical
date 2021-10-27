function create_phantom(fname, N)

coord = linspace(-1,1,N);
[z,y,x] = ndgrid(coord, coord, coord);

img = 0.0002*ones(N,N,N);
zshift = -0.5;

% Generate "a" ellipse
x0 = 0;     % Centers
y0 = 0;
z0 = zshift + 0;
ax = 0.69;  % major and minor axes
ay = 0.92;
az = 0.43;
img(((x-x0)/ax).^2 + ((y-y0)/ay).^2 + ((z-z0)/az).^2 <= 1 ) = 0.06;

% Generate "b" ellipse
x0 = 0;
y0 = -0.0184;
z0 = zshift+0;
ax = 0.60;
ay = 0.82;
az = 0.33;
img(((x-x0)/ax).^2 + ((y-y0)/ay).^2 + ((z-z0)/az).^2 <= 1) = 0.04;

% Generate region "e"
x0 = 0;
y0 = 0.2500;
z0 = 0.1+zshift;
ax = 0.21;
ay = 0.20;
az = 0.20;
img(((x-x0)/ax).^2 + ((y-y0)/ay).^2 + ((z-z0)/az).^2 <= 1) = 0.058;

% Generate region "d" (tilted ellipse)
x0 = 0;
y0 = 0;
z0 = 0 + zshift;
ax = 0.16;
ay = 0.31;
az = 0.1;
relative_value = 0.03;
temp_img = ((x-x0)/ax).^2 + ((y-y0)/ay).^2 + ((z-z0)/az).^2 <= 1;
for iz=1:N,
  temp_img2=imrotate(squeeze(temp_img(iz,:,:)),36,'bilinear','crop');
  shift=round(0.22/2*N);  % number of pixels to shift ellipse
  temp_img2=[temp_img2(:,shift+1:N), zeros(N,shift)];  
  temp_img2=[temp_img2(shift+1:N,:); zeros(shift,N)];  
  img(iz,:,:) = squeeze(img(iz,:,:)) - relative_value*temp_img2(:,:);
end

% Generate line
z0 = zshift;
img((abs(y+0.3-(1.5*(x-0.15)))<=10/N) & (y>-0.3) & (y<0) & ((z-z0)>-0.1) & ((z-z0)<0)) = 0.01;

xc = 0;
yc = 0;
zc = N/2;
Del_xy = 1.0;
Del_z=1.0;

res=[Del_xy Del_xy Del_z];
loc=[ xc-(N-1)*Del_xy/2  yc-(N-1)*Del_xy/2  zc-(N-1)*Del_z/2 ];
img=1000*img/(.02-.0000226);
write_vjk(img,fname,res,loc);

return

% Output phantom as .img binary file 
fp = fopen(fname, 'w');
fprintf(fp, '%d %d %d \n', N, N, N);
fprintf(fp,'%f %f %f \n', xc, yc, zc);
fprintf(fp,'%f %f \n', Del_xy, Del_z);
fprintf(fp,'%f \n', rI);
for z = 1:N
	for y = 1:N
		for x = 1:N
			fprintf(fp,'%f ', img(z,y,x));
		end
	end
end
fclose(fp);
