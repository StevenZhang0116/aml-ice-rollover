%--------------------------------------------------------------------------
% Rigid-body/theta-L simulation of floating ice
%
% Scott Weady, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------
clear
close all

addpath('functions');
target = 'save';
if isfolder(target); rmdir(target,'s'); end
mkdir(target)

%% System parameters and time stepping
global H beta rho_air rho_ice rho_water eps ds

% Parameters
R = 0.038; %radius of ice (m)
beta = 2.8e-5; %melt rate (m/s)
eps = 5e-8; %Gibbs-Thomson coefficient (m^2/s)
rho_air = 0; rho_water = 998; rho_ice = 917; %density (kg/m^3)
g = 9.8; %gravity (m/s^2)
dR = 0.01; %angular drag (1/s)
dT = 0.1; %velocity drag (1/s)

% Temporal discretization
tf = 60*60; %final time (s)
dt = 0.01; %time step (s)
tplt = 1; %plotting time (s)
Nt = floor(tf/dt); %number of time steps
nplt = floor(tplt/dt+0.1); %plotting times
nf = 0; %start frame counter
nt = 0; %start time step counter

%% Spatial discretization and initialization

% Setup grid (visualization only)
H = 3*R; Lx = 2*H; Ly = H + R; %box dimensions
dX = 1.5*R; dY = 2*R;

% Spectral discretization
L = 2*pi*R; %arclength 
Ns = 2048; %number of grid points
s = linspace(0,1,Ns+1)'; s = s(1:end-1); ds = s(2) - s(1); %arclength
iks = 2*pi*1i * [(0:Ns/2),(-Ns/2+1:-1)]'; %fourier modes
ik0 = iks; ik0(1) = 1;
filter = abs(iks) < 2*max(abs(iks))/3; %2/3 anti-aliasing filter
kinv_h = 1./iks; kinv_h(1) = 0;

% Initialize tangent angle and curvature
a = 1; %major axis of ellipse
arg = 2*pi*s - pi/2; %arclength parametrization
theta = pi - atan2(a*cos(arg),sin(arg)); %tangent angle
kappa = 2*pi*a*csc(arg).^2./(1 + a^2*cot(arg).^2); %curvature
kappa(isnan(kappa)) = 2*pi/a; %correct singularity

% Form body from tangent angle
rhsx_h = fft(cos(theta));
rhsy_h = fft(sin(theta));
x = real(ifft(kinv_h.*rhsx_h));
y = real(ifft(kinv_h.*rhsy_h));
X = L*[x y]; X = X - X(1,:);
X(:,2) = X(:,2) - mean(X(:,2)) + 0.756872*H;

% Adjust arclength
Xs = real(ifft(iks.*fft(X)));
L = sum(sqrt(sum(Xs.^2,2)))*ds; L0 = L;

% Store
body = struct('X',X,'kappa',kappa,'L',L,'theta',theta,'theta0',theta(1),'ohm',0,'L0',L);
f_m1 = struct('x0',0,'y0',0,'kappa',zeros(Ns,1),'L',0);
body_m1 = body;

% Euler step rigid body
A = calcarea(body,iks); A0 = A; Asave = [0 A]; %area
Vcm = [0,0]; %center of mass velocity
el = 0; %angular momentum
Fb = -g*buoyancy(body,iks); %buoyancy force
tau = -g*torque(body,iks); %torque
fvcm_m1 = [0,Fb]/(rho_ice*A) - (L/A)*dT*Vcm; %rhs for Vcm
fel_m1 = tau - dR*el; %rhs for angular momentum
Vcm_temp = Vcm + dt*fvcm_m1;
body.X = body.X + dt*Vcm_temp;

%% Visualization

sky = [189 236 240]/255;
water = [70 136 194]/255;
ice = [1 1 1];
interface = 0.33*[1 1 1];

fig = figure;
fig.Color = 'w';
fig.Units = 'inches';
fig.PaperSize = [8 7];
fig.PaperPosition = [0 0 fig.PaperSize];
fig.Position = [1 1 0 0] + fig.PaperPosition;

hold on
fill([0 Lx Lx 0 0],[H H Ly Ly H],sky);
fill([0 Lx Lx 0 0],[0 0 H H 0],water);
plot([0 Lx],[H H],'k','LineWidth',2)
axis equal tight, box on

axis([Lx/2-dX Lx/2+dX H-dY Ly])
xticklabels([]); yticklabels([]);

X = body.X; 
Xcm = centroid(body,iks);
dX = Xcm(1) - Lx/2;
XX = X; XX(:,1) = X(:,1) - dX; 
XX = [XX;XX(1,:)];
fh = fill(XX(:,1),XX(:,2),ice,'EdgeColor',interface,'LineWidth',2);
sh = scatter(Xcm(1)-dX,Xcm(2),20,'o','filled','MarkerEdgeColor','r',...
      'MarkerFaceColor','r');
ph = plot([XX(1,1) Lx/2],[XX(1,2) Xcm(2)],'r');

ax = gca;
ax.Units = 'inches';
ax.Position([1 2]) = [1 0.5];
ax.Position(3) = 0.99*(fig.PaperSize(1) - 2*ax.Position(1));
ax.Position(4) = diff(ax.YLim)/diff(ax.XLim)*ax.Position(3);

drawnow

timeset = [];
aa0set = [];

%% Main loop
while nt <= Nt && A/A0 > 0.01
  

  t = dt*nt;
  
  % Center of mass update
  L = body.L; %arclength
  A = calcarea(body,iks); %area
  Fb = -g*buoyancy(body,iks); %buoyancy  
  fvcm = [0,Fb]/(rho_ice*A) - (L/A)*dT*Vcm; %rhs for Vcm
  [Vcm,fvcm_m1] = ab2(Vcm,fvcm,fvcm_m1,dt); %update center of mass

  % Angular momentum update
  tau = -g*torque(body,iks); %torque
  fel = tau/(rho_ice*A) - (L/A)*dR*el; %rhs for angular momentum
  [el,fel_m1] = ab2(el,fel,fel_m1,dt); %update angular momentum
  body.ohm = rho_ice*A*el/inertia(body,iks); %angular velocity
  
  % Boundary update
  [body,body_m1,f_m1] = kappaL(body,body_m1,f_m1,Vcm,dt,s,iks,kinv_h,filter);
  A = calcarea(body,iks); Asave = [Asave;t,A];

  % Plot
  if mod(nt,nplt)==0
    
    clc
    mins = floor(t/60); secs = mod(t,60);
    fprintf('t = %2d min %2.0fs, A/A0 = %1.4f\n',mins,secs,A/A0)

    timeset(end+1) = mins*60+secs;
    aa0set(end+1) = A/A0;

    X = body.X; 
    Xcm = centroid(body,iks);
    dX = Xcm(1) - Lx/2;
    XX = X; XX(:,1) = X(:,1) - dX; XX = [XX;XX(1,:)];
    delete(fh), delete(ph), delete(sh)
    sub = 1:floor(length(X)/512):length(X);
    fh = fill(XX(sub,1),XX(sub,2),ice,'EdgeColor',interface,'LineWidth',2);
    sh = scatter(Lx/2,Xcm(2),20,'o','filled','MarkerEdgeColor','r','MarkerFaceColor','r');
    ph = plot([XX(1,1) Lx/2],[XX(1,2) Xcm(2)],'r');
    drawnow
        
    filename = sprintf('f-%d',nf);
%     save(strcat(target,'/',filename,'.mat'),'t','body');
%     print(strcat(target,'/',filename),'-dpng','-r300');
    nf = nf + 1;
    
  end

  nt = nt + 1;
  
end
