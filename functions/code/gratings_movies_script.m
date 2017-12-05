% GENERATE GRATINGS

clc; 
clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % % --------------- spherical correction caluculations --------------- % % 
 
W=104.3;  % width of screen, in cm
H=58.8;
        
Eyepos=[0.5,0.5];
        
xEye = Eyepos(1) * W;
yEye = Eyepos(2) * H;
zEye = 30;

warpgridsizeW=1920;
warpgridsizeH=1080;

xgrid=warpgridsizeW;
ygrid=warpgridsizeH;
equalDistanceX = linspace(0, W, xgrid);
equalDistanceY = linspace(0, H, ygrid);
        
xc=zeros(warpgridsizeH,warpgridsizeW);
yc=zeros(warpgridsizeH,warpgridsizeW);
for j=1:warpgridsizeH
xc(j,:) = equalDistanceX - xEye;
end
for j=1:warpgridsizeW
yc(:,j) = -(equalDistanceY - yEye);
end
      
% pass in spherical coordinates
r = sqrt((xc).^2 + (yc).^2 + (zEye).^2);

azimuth = atan(xc ./ zEye);
altitude = asin(yc ./ r);

% calculate the texture coordinates: 
tx = zEye * (1 + xc ./ r)- zEye;
ty = zEye * (1 + yc ./ r) - zEye;
            
% calculate centralAngle
centralAngle=zeros(warpgridsizeH,warpgridsizeW);
for i=1:warpgridsizeH
    for j=1:warpgridsizeW
        
centralAngle(i,j) = acos (cos(altitude(i,j)) * cos(abs(azimuth(i,j))));   
            
    end
end

% distance from eyepoint to texture vertex
arcLength = centralAngle .* zEye;
% remap the texture coordinate
theta = atan2(ty, tx);
            
for i=1:warpgridsizeH
      for j=1:warpgridsizeW
        
tx(i,j) = arcLength(i,j) .* cos(theta(i,j));
ty(i,j) = arcLength(i,j) .* sin(theta(i,j));
            
     end
end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % % ----------------- gamma correction caluculations ----------------- % % 
 
% dati
pixmesureG=(0:8:256);
lummesureG=[0.306,0.317,0.299;0.303,0.269,0.284;0.297,0.301,0.305;0.301,0.282,0.284;0.304,0.295,0.292;0.358,0.358,0.352;0.968,0.980,0.950;1.928,1.933,1.921;3.338,3.309,3.269;4.803,4.822,4.867;6.632,6.697,6.681;9.444,9.521,9.464;13.28,13.58,13.41;18.38,18.28,18.94;24.75,24.78,24.18;31.44,31.95,31.91;39.05,39.07,39.10;46.93,47.27,47.54;55.12,56.25,56.02;65.30,65.65,65.70;76.04,77.81,76.35;88.06,86.02,86.99;97.48,97.15,98.60;109.0,109.8,110.4;120.9,120.2,120.5;133.1,133.1,133.1;143.7,143.5,143.9;156.6,156.9,156.4;170.8,169.3,169.5;180.9,181.8,181.7;192.2,191.9,191.0;203.0,204.1,204.0;212.9,213.7,214.3];
lumG=mean(lummesureG,2);

% gamma
gamma = @(p,xdata)p(1)*(xdata).^p(2);

% fitting
p0 = [1 1];
[optpar,resnorm,~,exitflag,output] = lsqcurvefit(gamma,p0,pixmesureG',lumG);

% inverse gamma
invgamma = @(p,y)(y/p(1)).^(inv(p(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % % ------------- initalize grating parameters and grid -------------- % % 

% choose folder
rootfold=[uigetdir,filesep];

% choose bicodes
bitcode  = (0:1:108)+100;

% spatial parameters
sf_deg = [0.02,0.04,0.08];
dir = (0:30:330);

% convert degs into pixels
struct.res         = [1920 1080]; % monitor resolution
struct.sz          = [104.3 58.8]; % monitor width and height
struct.vd          = 30; % viewing distance
[pix_deg, deg_pix] = VisAng(struct); % calculate pixels per deg

% frame size (fullscreen)
grat_size_pxl1     = 1920;
grat_size_pxl2     = 1080;

% timing parameters
RefreshRate          = 30; % in Hz
StimDurationSecs     = 0.9; % in s ... era 1
tfrequencieslabels = [2,4,8]; % label for temporal frequencies 2, 4, 8 Hz
NumFramesPerPeriod   = [16,8,4]; % corresponding to temporal frequencies 2, 4, 8 Hz

% create a grid 
[x,y] = meshgrid((-grat_size_pxl1/2):(grat_size_pxl1/2)-1, (-grat_size_pxl2/2):(grat_size_pxl2/2)-1);
% meshgrid y axis is inverted!!
y=-y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------- initalize grating parameters and grid --------------- % % 

counter=0;

for looptf=1:length(NumFramesPerPeriod)  % loop over temporal frequencies
   
StimDurationFrames = round(StimDurationSecs * RefreshRate);
StimFrameIndices   = mod(0:(StimDurationFrames-1), NumFramesPerPeriod(looptf)) + 1;
    
for loopsf=1:length(sf_deg) % loop over spatial frequencies

% current spatial frequency in pixels
sf_pxl  = (sf_deg(loopsf)/pix_deg(1))*2*pi;

message=['Start producing sf = ',num2str(sf_deg(loopsf)),'  tf = ',num2str(tfrequencieslabels(looptf)), ' ...\n'];
fprintf(message)

for i = 1:length(dir) % loop over directions 
tic
counter=counter+1;

% preallocation	of TextFrame
TextFrame=zeros(size(x,1),size(x,2),NumFramesPerPeriod(looptf));

% calculate wave vector corresponding to current direction
    Ang = (dir(i).*pi/180);
	K = [cos(Ang), sin(Ang)]*sf_pxl;
    
    for k=1:NumFramesPerPeriod(looptf) % loop on phase steps
		
        Phase = (k/NumFramesPerPeriod(looptf))*2*pi; % current phase
		TextFrame(:,:,k) = 0.5*sin(K(1)*x + K(2)*y - Phase)+0.5; % draw gratings
        TextFrame(:,:,k) = interp2(xc,yc , squeeze( TextFrame(:,:,k)),tx,ty); % spherically correct
	    linlum=(lumG(end).*squeeze( TextFrame(:,:,k))+lumG(1));
        TextFrame(:,:,k) = invgamma(optpar,linlum)./255; % gamma correct
        
    end
		
% Use repmat to reach desired duration 
HowManyReps = round(StimDurationFrames/NumFramesPerPeriod(looptf));
CompleteStimulus = repmat(TextFrame,[1,1,HowManyReps]);

% Per farlo da 0.9 secondi...
Fraction=StimDurationFrames/NumFramesPerPeriod(looptf)-round(StimDurationFrames/NumFramesPerPeriod(looptf));
NumFramesExcess=round(NumFramesPerPeriod(looptf)*Fraction);
if NumFramesExcess<0
CompleteStimulus =CompleteStimulus(:,:,1:NumFramesExcess+end);
else
CompleteStimulus =cat(3,CompleteStimulus,CompleteStimulus(:,:,1:NumFramesExcess));   
end

%% save each frame in its own folder

% create folder for the current grating and jump in
foldername=[rootfold,'0',num2str(bitcode(counter)),'_grating_',num2str(dir(i)),'_',num2str(sf_deg(loopsf)),'_',num2str(tfrequencieslabels(looptf))];
mkdir(foldername)
oldfolder=cd(foldername);

    for j = 1:size(CompleteStimulus,3) % loop over frames
        if j <= 9
            imname = sprintf('%sgrating_%d_frame_00%d.bmp', [foldername,filesep], dir(i), j);
        elseif j < 100
            imname = sprintf('%sgrating_%d_frame_0%d.bmp', [foldername,filesep], dir(i), j);
        elseif j >=100
            imname = sprintf('%sgrating_%d_frame_%d.bmp', [foldername,filesep], dir(i), j);
        end
        imwrite(CompleteStimulus(:,:,j),imname,'bmp'); % save current frame
    end
    
% jump out
cd(oldfolder);
	message=['direction ',num2str(i),' produced\n'];
fprintf(message)
toc	

end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
