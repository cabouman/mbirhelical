


dataFolder = '/run/media/ch194925/2C660A42660A0CF0/CT_ConeBeam/Head_Scan/BC_Under_3_Years_ptr/';

load([dataFolder, 'dataInfo'])  
load([dataFolder, 'dataInfo_UID20'])  
load([dataFolder, 'dataInfo_UID50'])  

flying_focal_spot = 4;
id_high_low = 1; % id_high_low: 1: high; 2: low
viewStart = 2;% 100;
%viewEnd = dataInfo.ScanDescr.FramesPerRotation*4*8;% *4 is 360 degrees one rotation
viewEnd = dataInfo.ScanDescr.FramesPerRotation*flying_focal_spot*2;
Proj= zeros((viewEnd)/flying_focal_spot, dataInfo.ModeParXML.Type.ChnNum(id_high_low), dataInfo.ModeParXML.ModePar.NoOfSlicesDMS);
j= 1;
for i = viewStart:flying_focal_spot:viewEnd
projName = ['proj_', num2str(i)];
    headerName = ['header_', num2str(i)];
    load([dataFolder, projName], 'proj')
    Proj(j, :, :) = proj{id_high_low};
    j = j+1;
end
clear temp;


Proj = permute(Proj,[3, 2, 1]);

%Proj = Proj - min(Proj(:));
temp = squeeze(Proj(32,1600:1800,1:300));
offset = mean(temp(:));

Proj = Proj - offset;

ratio_factor= max(Proj(:))/3;
Proj = Proj ./ ratio_factor;



fp=fopen('head_scan.sino','w');

fprintf(fp,'%d %d %d\n',dataInfo.ModeParXML.ModePar.NoOfSlicesDMS, dataInfo.ModeParXML.Type.ChnNum(id_high_low),(viewEnd)/flying_focal_spot);


for i = 1:size(Proj,3)
    proj = squeeze(Proj(:,:,i)); 
    fwrite(fp,proj,'float32');
end    
fclose(fp);


% Proj=mat2gray(Proj);
% implay(Proj)
