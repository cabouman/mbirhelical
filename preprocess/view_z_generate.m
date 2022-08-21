clear all;

for data_idx=0:199
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


	z_center_list=zeros(NumViews,1);
	output_view_angles=zeros(NumViews,1);


	for iv=1:NumViews
		if NumViews <10000	
	    		proj=sprintf('proj_%04d.dcm',iv);
		else
			proj=sprintf('proj_%05d.dcm',iv);
		end
		proj_name =[dataset,proj];
		temp_info = dicominfo(proj_name,'dictionary','./DICOM-CT-PD-dict_v10.txt');
	    
		z_center_list(iv,1)=temp_info.DetectorFocalCenterAxialPosition;
	    
	    	output_view_angles(iv,1)=temp_info.DetectorFocalCenterAngularPosition;
		disp(iv)
	end

	for i=1:NumViews
		output_view_angles(i,1) = output_view_angles(i,1)-pi*3/2;
	end


	output_path_view=sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/preprocess/dcm%03d_viewAnglesList.txt',data_idx);

	output_path_z=sprintf('/global/cscratch1/sd/wang1698/AAPM_2022/preprocess/dcm%03d_zPositionList.txt',data_idx);

	dlmwrite(output_path_view,output_view_angles);
	dlmwrite(output_path_z,z_center_list);


	%dlmwrite('C:/Users/xiaow/Downloads/AAPM-2022/dcmproj_copd/dcm_000_preprocess/dcm000_viewAnglesList.txt',output_view_angles);
	%dlmwrite('C:/Users/xiaow/Downloads/AAPM-2022/dcmproj_copd/dcm_000_preprocess/dcm000_zPositionList.txt',z_s_list);
end
