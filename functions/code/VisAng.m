function [pixperdeg, degperpix]=VisAng(params)
%
% Takes as input a structure containing Monitor parameters:
%
% struct.res - the resolution of the monitor
% struct.sz - the size of the monitor in cm
% struct.vd - viewing distance in cm
%
% Es:
% Monitor Sony
% struct.res=[1024 768];
% struct.sz=[33 24];
% struct.vd=[57];


% Monitor Asus
% struct.res=[1280 800];
% struct.sz=[33.1 20.7];
% struct.vd=[114];

%
% (these values can either be along a single dimension or
% for both the width and height)
% Stimulus parameters
% struct.pixs - the size of the stimulus
% (these values can either be along a single dimension or
% for both the width and height)
% struct.vdist - the viewing distance in cm.
%
% Calculates the visual angle subtended by a single pixel
%
% Returns the pixels per degree
% and it's reciprocal - the degrees per pixel (in degrees,
% not radians)
%
% 1) Calculate the width of a single pixel in cm.
% (i) Measure the width and height of the screen (there is usually a black rim, don?t
% measure that) in cm
% (ii) Find out the number of pixels for that monitor horizontally and vertically
% width of single pix in cm=width screen/num horizontal pixels
% height single pixel in cm=height screen/num vertical pixels

% It?s possible that your pixels aren?t square. It may be that you don?t care. If you
% do, then adjust your monitor to make the pixels squarer. If your pixels still aren?t
% square then you may have to have different values for the number of pixels per
% degree in the horizontal and vertical directions.
% 2) Calculate the visual angle subtended by a single pixel
% 
% theta=atan(pix/vdist); 
% 
% This gives you two values theta=degperpix, and the reciprocal, 1/theta=pixperdeg
% written by IF 7/2000
vdist=params.vd;
pix=params.sz./params.res; %calculates the size of a pixel in cm
degperpix=(2*atan(pix./(2*vdist))).*(180/pi);
pixperdeg=1./degperpix;

cmperdeg=pix.*pixperdeg;
degpercm=1./cmperdeg;

return;

