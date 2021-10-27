function create_twoballs(fname, N)
% NxNxN phantom
% 500mm each dimension

coord = linspace(-1,1,N);
[z,y,x] = ndgrid(coord,coord,coord);

img = 0.0002*ones(N,N,N);
zshift = -0.5;

x0 = 0.3;
y0 = 0.3;
z0 = zshift + 0;
img(((x-x0).^2 + (y-y0).^2 + (z-z0).^2) <= 0.2) = 0.08;

x0 = -0.4;
y0 = -0.4;
z0 = zshift + 0;
img(((x-x0).^2 + (y-y0).^2 + (z-z0).^2) <= 0.2) = 0.01;

xc = 0;
yc = 0;
zc = N/2;
Del_xy = 1.0;
Del_z = 1.0;
rI=500;

res=[Del_xy Del_xy Del_z];
loc=[ xc-(N-1)*Del_xy/2  yc-(N-1)*Del_xy/2  zc-(N-1)*Del_z/2 ];
img=1000*img/(.02-.0000226);
write_vjk(img,fname,res,loc);

