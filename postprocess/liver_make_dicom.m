clear all;
cd ../matlab/


for data_idx=134:199
	if data_idx <67
		dataset =sprintf('/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_%03d/',data_idx)
	elseif data_idx >=67 && data_idx <134
		dataset =sprintf('/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_lung_lesion/dcm_%03d/',data_idx)
	else
		dataset =sprintf('/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_liver/dcm_%03d/',data_idx)
	end



	a = dir([dataset,'*.dcm']);
	NumViews = numel(a);
	fprintf('NumViews is %d ',NumViews);


	WaterAttenuationCoefficient = 0.019922;

	disp(WaterAttenuationCoefficient)

	%dimension z x y
	reconname =sprintf('/gpfs/alpine/gen006/proj-shared/xf9/recon/dcm%03d/reconmyid0.vjk',data_idx);

	disp(reconname)

	recon_myid0 = read_vjk_float(reconname);
	recon_myid0(recon_myid0<0)=0;
	recon_myid0 = (recon_myid0-WaterAttenuationCoefficient)*1000/WaterAttenuationCoefficient;
	recon_myid0=recon_myid0(:,end:-1:1,:);
	%liver first, then lung


	%removes artifacts from negative sinogram entires
	mask = (recon_myid0>-1000);

	for j=1: size(recon_myid0,2)
		for k=1:size(recon_myid0,3)
			if (((j-256)*(j-256)/(230^2)+(k-256)*(k-256)/(230^2))<1)
				mask(:,j,k)=0;
			end
		end
	end
	recon_myid0(mask)=-1000;


	zname =sprintf('/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_zPositionList.txt',data_idx);

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





	recon_parameter =sprintf('/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/info_recon.txt',data_idx);

	fid=fopen(recon_parameter); 	
	% set linenum to the desired line number that you want to import
	linenum = 20;
	% use '%s' if you want to read in the entire line or use '%f' if you want to read only the first numeric value
	z_center = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	z_center = z_center{1};
	fseek(fid,0,'bof');
	linenum = 23;
	inplane_spacing = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	inplane_spacing = inplane_spacing{1};
	fseek(fid,0,'bof');

	linenum = 26;
	z_spacing = textscan(fid,'%f',1,'delimiter','\n', 'headerlines',linenum-1);
	z_spacing = z_spacing{1};
	fseek(fid,0,'bof');
	
	fprintf('z_center %f, z_spacing %f inplane_spacing %f\n ',z_center,z_spacing,inplane_spacing);

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


	uid = dicomuid;

	for i=1:size(recon_myid0,1)
		temp_img = int16(squeeze(recon_myid0(i,:,:)));
		patient =sprintf('dcm_%03d',data_idx);

		dicomname =sprintf('/gpfs/alpine/gen006/scratch/xf9/dicom/dcm%03d/%04d.dcm',data_idx,i);
	
	    	dicomwrite(temp_img, dicomname);

		temp_info=dicominfo(dicomname);

		temp_info.PatientName = patient;

		temp_info.(dicomlookup('0018', '0050'))=z_spacing;
		temp_info.(dicomlookup('0018', '1050'))=inplane_spacing;
		temp_info.(dicomlookup('0020', '1041'))=recon_z_start +(i-1)*z_spacing;
		temp_info.(dicomlookup('0020', '0011'))= double(0); %series number
		temp_info.(dicomlookup('0020', '000E'))= uid; %series instance number

	    	dicomwrite(temp_img, dicomname, temp_info,'CreateMode','copy');

	end

end
