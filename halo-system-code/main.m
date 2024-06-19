%--------------------------------------------------------------------------
% Halo system for the ice flipping experiments
%
% Crafted by Scott Weady
% Revised by Zihan Zhang, Courant Institute
% Updated July 2023
%--------------------------------------------------------------------------

warning('off')
close all
set(0,'units','pixels');
res = get(0,'screensize'); %full screen resolution

func_result = allfunction();

global Nx Ny params im
Nx = res(4); Ny = res(3);
R = 0.5; %initial radius
X = 0; %initial center
Y = 0; %initial bottom
A = 0.1; %decay rate
params = struct('X',X,'Y',Y,'R',R,'A',A);

% Setup for controller
x = 120; y = 180;
dx = 120;
xlen = 200; ylen = 20;
pix = 1/Nx;
step = [pix 10*pix];
ui_bg = 'w';

ctrl = figure(1);
set(ctrl,'Position',[500 500 2*xlen 12*ylen],'Color',ui_bg,...
         'MenuBar','none','ToolBar','none')

fn = 'Cambria'; %font
fs = 12; %fontsize

% Radius slider
str = 'Sharpness'; 
loc = [x y xlen ylen];
uicontrol('style','slide','Position',loc,'Value',R,'SliderStep',step,...
                          'CallBack',{@func_result.changeDecay},'BackgroundColor',ui_bg,...
                          'min',0,'max',1);
loc = [x-dx y dx 20];
uicontrol('style','text','FontName',fn,'Position',loc,'String',str,...
                         'FontSize',fs,'BackgroundColor',ui_bg);

str = 'Radius';
loc = [x y-2*ylen xlen ylen];
uicontrol('style','slide','Position',loc,'Value',R,'SliderStep',step,...
                          'CallBack',{@func_result.changeRadius},'BackgroundColor',ui_bg);
loc = [x-dx y-2*ylen dx 20];
uicontrol('style','text','FontName',fn,'Position',loc,'String',str,...
                         'FontSize',fs,'BackgroundColor',ui_bg);

% Horizontal translator
str = 'X';
loc = [x y-4*ylen xlen ylen];
uicontrol('style','slide','Position',loc,'Value',X,'SliderStep',step,...
                          'CallBack',{@func_result.changeX},'BackgroundColor',ui_bg,...
                          'min',-1,'max',1);                    
loc = [x-dx y-4*ylen dx 20];
uicontrol('style','text','FontName',fn,'Position',loc,'String',str,...
                         'FontSize',fs,'BackgroundColor',ui_bg);

% Vertical translator
str = 'Y';
loc = [x y-6*ylen xlen ylen];
uicontrol('style','slide','Position',loc,'Value',Y,'SliderStep',step,...
                          'CallBack',{@func_result.changeY},'BackgroundColor',ui_bg,...
                          'min',-(Ny/Nx),'max',(Ny/Nx));
loc = [x-dx y-6*ylen dx 20];
uicontrol('style','text','FontName',fn,'Position',loc,'String',str,...
                         'FontSize',fs,'BackgroundColor',ui_bg);          
% Initialize figure
fig = figure(2);
set(fig,'Units','normalized','OuterPosition',[0.1 0.1 0.4 0.4],'Color','k',...
        'MenuBar', 'none','ToolBar', 'none')
      
I = func_result.calcI(X,Y,R,A);
im = imshow(I);
func_result.updateFig(I)










