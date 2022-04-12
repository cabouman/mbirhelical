clear all;


NumViews=9000;   
FramesPerRotation = 1000;

%NumColumns is the number of channels. NumRows is the number of detector
%rows
rho_0_list=zeros(NumViews,1);
z_center_list=zeros(NumViews,1);
phi_0_list=zeros(NumViews,1);

delta_rho_list=zeros(NumViews,1);
delta_z_list=zeros(NumViews,1);
delta_phi_list=zeros(NumViews,1);

for iv=1:NumViews
    dataset=sprintf('C:/Users/xiaow/Downloads/AAPM-2022/dcmproj_copd/dcm_000/proj_%04d.dcm',iv);
    temp_info = dicominfo(dataset,'dictionary','C:/Users/xiaow/Downloads/AAPM-2022/DICOM-CT-PD-dict_v10.txt');
    rho_0_list(iv,1) = temp_info.DetectorFocalCenterRadialDistance;
    delta_rho_list(iv,1)=0;
    
    
    z_center_list(iv,1)=temp_info.DetectorFocalCenterAxialPosition;
    delta_z_list(iv,1)=0;
    
    
    phi_0_list(iv,1)=temp_info.DetectorFocalCenterAngularPosition;
    delta_phi_list(iv,1)=0;
end

x_center_list = - rho_0_list .* sin(phi_0_list);
y_center_list = rho_0_list .* cos(phi_0_list);

x_s_list = - (rho_0_list + delta_rho_list) .* sin(phi_0_list + delta_phi_list);
y_s_list = (rho_0_list + + delta_rho_list) .* cos(phi_0_list + delta_phi_list);
z_s_list = z_center_list + delta_z_list;



%simulation now

phi_0_simulated = zeros(NumViews,1);
j = 0;
for i=1:NumViews
    if(i~=1 && phi_0_list(i,1)>phi_0_list(i-1,1))
        j=j+1;
    end
        phi_0_simulated(i,1) = phi_0_list(i,1)-pi*3/2;
end




output_view_angles =zeros(NumViews,1);
for i=1:NumViews
    output_view_angles(i) = phi_0_simulated(i) + delta_phi_list(i);
end

dlmwrite('C:/Users/xiaow/Downloads/AAPM-2022/dcmproj_copd/dcm_000_preprocess/dcm000_viewAnglesList.txt',output_view_angles);
dlmwrite('C:/Users/xiaow/Downloads/AAPM-2022/dcmproj_copd/dcm_000_preprocess/dcm000_zPositionList.txt',z_s_list);

