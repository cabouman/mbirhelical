D=595/(100/256);    %575 is the distance from source to object. 500/512 is the pixel size. D is the number of pixels from source to the center of the object
detector_angle = 1.023*920/1085.6;  % the detector angle in radius
dsensor = detector_angle/920;   % the detector spacing

Ifan3 = ifanbeam(circshift(sino,0),D, 'FanSensorSpacing',dsensor*180/pi,'FanSensorGeometry','arc','FanRotationIncrement',180/(1050),'OutputSize',920);
imagesc(Ifan3')
colormap(gray)

%For CT Head Scan, a rotation is 90 degrees. So 1050 is for 90 degrees