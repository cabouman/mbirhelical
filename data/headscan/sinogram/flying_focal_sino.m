dataFolder = '/run/media/ch194925/2C660A42660A0CF0/CT_ConeBeam/Head_Scan/BC_Under_3_Years_ptr/';

load([dataFolder, 'dataInfo'])  
load([dataFolder, 'dataInfo_UID20'])  
load([dataFolder, 'dataInfo_UID50'])  

flying_focal_spot = 4;
id_high_low = 1; % id_high_low: 1: high; 2: low
viewStart = 1;% 100;
viewEnd = dataInfo.ScanDescr.FramesPerRotation*4*2;% *4 is 360 degrees one rotation
Proj= zeros((viewEnd)/flying_focal_spot, dataInfo.ModeParXML.Type.ChnNum(id_high_low)*2, dataInfo.ModeParXML.ModePar.NoOfSlicesDMS*2);
sinogram_focal1= zeros((viewEnd)/flying_focal_spot, dataInfo.ModeParXML.Type.ChnNum(id_high_low), dataInfo.ModeParXML.ModePar.NoOfSlicesDMS*2);
j= 1;
for i = viewStart:flying_focal_spot:viewEnd
    projName = ['proj_', num2str(i)];
    headerName = ['header_', num2str(i)];
    load([dataFolder, projName], 'proj')
    sinogram_focal1(j, :, 1:2:end) = proj{id_high_low};
    
    projName = ['proj_', num2str(i+2)];
    headerName = ['header_', num2str(i+2)];
    load([dataFolder, projName], 'proj')
    sinogram_focal1(j, :, 2:2:end) = proj{id_high_low};    
    
    j = j+1;
end


viewStart = 2;% 100;
sinogram_focal2= zeros((viewEnd)/flying_focal_spot, dataInfo.ModeParXML.Type.ChnNum(id_high_low), dataInfo.ModeParXML.ModePar.NoOfSlicesDMS*2);
j= 1;
for i = viewStart:flying_focal_spot:viewEnd
    projName = ['proj_', num2str(i)];
    headerName = ['header_', num2str(i)];
    load([dataFolder, projName], 'proj')
    sinogram_focal2(j, :, 1:2:end) = proj{id_high_low};
    
    projName = ['proj_', num2str(i+2)];
    headerName = ['header_', num2str(i+2)];
    load([dataFolder, projName], 'proj')
    sinogram_focal2(j, :, 2:2:end) = proj{id_high_low};    
    
    j = j+1;
end



for i = 1:1:size(Proj,2)
    if(mod(i,2)==1)
        Proj(:,i,:)=sinogram_focal1(:,(i+1)/2,:);
    else
        Proj(:,i,:)=sinogram_focal2(:,i/2,:);
    end
end

clear sinogram_focal1;
clear sinogram_focal2;

Proj = permute(Proj,[3, 2, 1]);

%Proj = Proj - min(Proj(:));
temp = squeeze(Proj(32,1600:1800,1:300));
offset = mean(temp(:));

Proj = Proj - offset;

ratio_factor= max(Proj(:))/3;
Proj = Proj ./ ratio_factor;


fp=fopen('head_scan.sino','w');

fprintf(fp,'%d %d %d\n',dataInfo.ModeParXML.ModePar.NoOfSlicesDMS*2, dataInfo.ModeParXML.Type.ChnNum(id_high_low)*2,(viewEnd)/flying_focal_spot);


for i = 1:size(Proj,3)
    proj = squeeze(Proj(:,:,i)); 
    fwrite(fp,proj,'float32');
end    
fclose(fp);


% Proj=mat2gray(Proj);
% implay(Proj)
