frameSizeX = 500;
frameSizeY = 500;
myImageCh1 = zeros(frameSizeY,frameSizeX,'int16');
myImageCh1(round(frameSizeX/2),round(frameSizeY/2)) = 255;
maxVoltageX = 3;
maxVoltageY = 3;
voltsPerPixel = 2*maxVoltageX/frameSizeX;
pixelRate = 1e6;%in Hz??





MpData.dt = 1.0/pixelRate;
MpData.startVoltageX = -1*maxVoltageX;
MpData.startVoltageY = -1*maxVoltageY;
MpData.stopVoltageX = maxVoltageX;
MpData.stopVoltageY = maxVoltageY;
MpData.voltsPerPixel = voltsPerPixel;
MpData.pixelDwellTime = MpData.dt;%redundant (historic)
MpData.pixelsToShiftX = 0;
MpData.pixelsToShiftY = 0;
MpData.zoomFactor = 1;
MpData.framesCh1 = myImageCh1;
%MpData.framesCh2 = ;

