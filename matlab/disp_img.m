function disp_img(fname)

% read file
[img res loc]=read_vjk(fname);

if(0)
fp = fopen(fname,'rb');
[Nx N] = fscanf(fp,'%d',1);
[Ny N] = fscanf(fp,'%d',1);
[Nz N] = fscanf(fp,'%d\n',1);
[xc N] = fscanf(fp,'%lf',1);
[yc N] = fscanf(fp,'%lf',1);
[zc N] = fscanf(fp,'%lf\n',1);
[Del_xy N] = fscanf(fp,'%lf',1);
[Del_z N] = fscanf(fp,'%lf\n',1);
[rI N] = fscanf(fp,'%lf\n',1);

% 3D memory allocation
img = zeros(Nx,Ny,Nz);
for jz = 1:Nz
	for jy = 1:Ny
		for jx = 1:Nx
			[img(jz,jy,jx) N] = fscanf(fp,'%f\n',1);
		end
	end
end
fclose(fp);
end

[slice_num, rows, cols] = size(img);

max_val = max(img(:));
min_val = min(img(:));
%max_val = 600;
%min_val = 0;


%initial display
sag = squeeze(img(:,:,floor(cols/2)));
ax = squeeze(img(floor(slice_num/2),:,:));
cor = squeeze(img(:,floor(rows/2),:));
mos = zeros((slice_num+rows),(rows+cols));
mos(0*slice_num+1:1*slice_num,0*rows+1:1*rows) = sag;
mos(0*slice_num+1:1*slice_num,1*rows+1:rows+cols) = cor;
mos(1*slice_num+1:slice_num+rows,1*rows+1:rows+cols) = ax;

figure(1)
colormap(gray(256))
imagesc(mos,[min_val,max_val])  % axial
colorbar
axis('image');

while(1)
	[c,r,button] = ginput(1);
	c = round(c);
	r = round(r);
	if ((c <= rows) & (r <= slice_num)) % input from sag view
		ax = squeeze(img(r,:,:));
		cor = squeeze(img(:,c,:));
	elseif ((r <= slice_num) & (c > rows)) % input from cor view
		c = c - rows;
		sag = squeeze(img(:,:,c));
		ax = squeeze(img(r,:,:));
	elseif ((r > slice_num) & (c > rows))
		r = r - slice_num;
		c = c - rows;
		sag = squeeze(img(:,:,c));
		cor = squeeze(img(:,r,:));
	end

	mos(0*slice_num+1:1*slice_num,0*rows+1:1*rows) = sag;
	mos(0*slice_num+1:1*slice_num,1*rows+1:rows+cols) = cor;
	mos(1*slice_num+1:slice_num+rows,1*rows+1:rows+cols) = ax;
	imagesc(mos,[min_val,max_val])
	colorbar
	axis('image');

	% save image
	%imwrite(uint8(mos), 'mos.tif', 'compression', 'lzw');
	%imwrite(mat2gray(sag, [min_val max_val]), 'sag.tif', 'compression', 'lzw');
	%imwrite(mat2gray(cor, [min_val max_val]), 'cor.tif', 'compression', 'lzw');
	%imwrite(mat2gray(ax, [min_val max_val]), 'ax.tif', 'compression', 'lzw');
end

return
