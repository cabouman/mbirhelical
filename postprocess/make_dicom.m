clear all;
cd ../matlab/


for data_idx=0:9
	if data_idx <67
		dataset =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/raw/dcmproj_copd/dcm_%03d/',data_idx)
	elseif data_idx >=67 && data_idx <134
		dataset =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/raw/dcmproj_lung_lesion/dcm_%03d/',data_idx)
	else
		dataset =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/raw/dcmproj_liver/dcm_%03d/',data_idx)
	end



	a = dir([dataset,'*.dcm']);
	NumViews = numel(a);
	fprintf('NumViews is %d ',NumViews);

	if NumViews <10000	
		proj=sprintf('proj_%04d.dcm',1);
	else
		proj=sprintf('proj_%05d.dcm',1);
	end
	proj_name =[dataset,proj];
	temp_info = dicominfo(proj_name,'dictionary','../preprocess/DICOM-CT-PD-dict_v10.txt');


	WaterAttenuationCoefficient = temp_info.WaterAttenuationCoefficient /10;

	disp(WaterAttenuationCoefficient)

	%dimension z x y
	reconname =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/recon/dcm%03d/reconmyid0.vjk',data_idx);

	disp(reconname)

	recon_myid0 = read_vjk_float(reconname);
	recon_myid0 = (recon_myid0-WaterAttenuationCoefficient)*1000/WaterAttenuationCoefficient;
	recon_myid0=recon_myid0(:,end:-1:1,:);


	zname =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/preprocess/dcm%03d_zPositionList.txt',data_idx);

	fid=fopen(zname); 	
	% set linenum to the desired line number that you want to import
	linenum = 501;
	% use '%s' if you want to read in the entire line or use '%f' if you want to read only the first numeric value
	recon_z_start = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	recon_z_start = recon_z_start{1};
	fseek(fid,0,'bof');
	linenum = NumViews-500;
	recon_z_end = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	fseek(fid,0,'bof');

	recon_z_end = recon_z_end{1};
	fprintf('recon_z_start %f, recon_z_end %f\n ',recon_z_start,recon_z_end);





	recon_parameter =sprintf('../data/aapm-parameters/dcm_%03d/info_recon.txt',data_idx);

	fid=fopen(recon_parameter); 	
	% set linenum to the desired line number that you want to import
	linenum = 20;
	% use '%s' if you want to read in the entire line or use '%f' if you want to read only the first numeric value
	z_center = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	z_center = z_center{1};
	fseek(fid,0,'bof');
	linenum = 26;
	z_spacing = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	z_spacing = z_spacing{1};
	fseek(fid,0,'bof');
	fprintf('z_center %f, z_spacing %f\n ',z_center,z_spacing);

	total_z_slices = size(recon_myid0,1);
	fprintf('total_z_slices %d\n',total_z_slices);



	if mod(total_z_slices,2)==0
		skipped_slices = fix((recon_z_start-(z_center-z_spacing/2-(total_z_slices/2-1)*z_spacing)+1e-4)/z_spacing);
	else
		skipped_slices = fix((recon_z_start-(z_center-fix(total_z_slices/2)*z_spacing)+1e-4)/z_spacing);

	end
	fprintf('skipped_slices %d\n',skipped_slices)



        useful_z_slices = ceil((recon_z_end - recon_z_start)/z_spacing)+1;

	fprintf('useful_z_slices %d\n',useful_z_slices);


	recon_myid0=recon_myid0((skipped_slices+1):(skipped_slices+useful_z_slices),:,:);
	
	fprintf('recon_myid0 number of slices after pruning %d\n',size(recon_myid0,1));

	for i=1:size(recon_myid0,1)
		temp_img = int16(squeeze(recon_myid0(i,:,:)));
		temp_info.slicelocation = recon_z_start +(i-1)*z_spacing;
		patient =sprintf('dcm_%03d',data_idx);

		temp_info.PatientName = patient;
		dicomname =sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/dicom/dcm%03d/%04d.dcm',data_idx,i);
	
	    	dicomwrite(temp_img, dicomname,temp_info,'SliceLocation',recon_z_start +(i-1)*z_spacing);

	end

end
