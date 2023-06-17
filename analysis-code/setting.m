%--------------------------------------------------------------------------
% Parameter settings and initilization section
%
% Steven Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------


%% analysis file
rootname = ['/Users/bndsteven/Dropbox/ice-melting-AML/rollover/ice-flipping/' ...
    'experiments-data/for-data/'];
foldername = '2022-06-28-b/';
filename = 'record.MOV';
calibration = 'calibration.MOV';

disp(['=== Analysis Experiment name: ',foldername(1:end-1),' ==='])

videoName = strcat(rootname,foldername,filename);
calibrationName = strcat(rootname,foldername,calibration);

v = VideoReader(videoName);

global rfr
% framerate registeration
rfr = v.FrameRate;
% start frame, when the iceberg is steady
sf = 60*rfr;    
% ignore the last 120 second
endFrame = min(v.NumFrames,1200*rfr);
% interval between frame
% when interval=rfr, meaning that 1 frame per second
interval = 1*round(rfr);

%% instantize the function
fr = func_aux();
intp = func_itp();
bdyn = func_bodydynamics();
cf = func_curve();

%% parameter settings
global rho_water rho_air rho_ice g method dec bkt

% interpolation method
method = 'F';
% boundary tracking method
dec = 'A';

% interpolation degree
degmin = 5;
degmax = 32;
variationset = degmin:1:degmax;
% number of interpolation point
numpt = 2048;

% define density parameters
rho_water = 998; % kg/m3
rho_ice = 917; % kg/m3
% at 20 Celsius and 1 atmosphere pressure
rho_air = 0; % kg/cm3
% gravity constant
g = 9.8;

% melting rate (in simulation) registeration
vr = 2.8*10^(-5); % m/s

% overall filter threshold
areaopen = 40; % pixel [free parameter]
% water surface around area filter threshold
waterareaopen = 50; % pixel, 
% alpha shape parameter
alphanum = 140; % dimensionless? [free parameter]

% adjustment in gradient filtering
left_adj = 5; % pixel
right_adj = 5; % pixel

% total number of points in drawing the circle
bkt = 360;

% 1-threstrap in calculation
threstrap = 0.1; 

% naive parameter of water surfacea, difference is ~45-55 [outdated]
waterSurfaceTop = 180;
waterSurfaceBot = 230;

tInv = sf:interval:endFrame;

% cropped parameter
startHeight = 450;
heightLength = 700;
cropped_interval = [0 startHeight 720 heightLength];

% cell to storing metadata at different time frames
cellLength = floor((endFrame-sf)/interval)+1;
init_grey_cell = cell(1,cellLength);
xytrackcell = cell(1,cellLength);
rthetatrackcell = cell(1,cellLength);
edgetrackcell = cell(1,cellLength);
alphatrackcell = cell(1,cellLength);
interpolation_cell = cell(1,cellLength);

% initialize different data structure
LTINV = length(tInv); 
frameLab = zeros(1,LTINV);
timeLab = zeros(1,LTINV);
areaSet = zeros(1,LTINV);
Hboundset = []; % water height using boundary integral
Hbodyset = zeros(1,LTINV); % water height using body integral 
cross1set = []; % boundary integral cross pt (num of equilibrium)
cross2set = zeros(1,LTINV); % body integral cross pt (num of equilibrium)
newtorquezeroset = [];
% shape dynamics
adjustedcurvatureSet = zeros(1,LTINV);
curvature3Set = zeros(1,LTINV);

%% load and edit with calibration video
rr = fr.calibration_video(calibrationName,1);

close all

% water surface level
[head2top] = fr.calibration_watersurface(videoName,cropped_interval,sf);
wl_cm_s = head2top/rr;
disp(['water surface height is: ', num2str(wl_cm_s)]);

close all

%% additional parameters 
comSet = NaN(1,2);
% last stage of trap line
lasttrap = [0,0];
prevarea = 0;
savefig = 1;
wlineset = [];
closeoverind = 1; 
colorVec = parula(cellLength);
adj = 0;
plotset = 1;
manuelw = 1; 

%% global data cell
natpara = struct('wcs',wl_cm_s,'rr',rr,'numpt',numpt,'vrr',variationset);










