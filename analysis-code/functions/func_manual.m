%--------------------------------------------------------------------------
% Manual operation on boundary data points
% Based on ginput2.m 
% Delete unwanted points -> Add Points -> Check
% Use mouse to operate: 
% Single left click: zoomin; 
% Double click: back to original view
% Right click: select datapoint and mark
% [The frame will be automatically centered based on the direction on the last two added
% points for better user experience]
% Left click (outside of frame) + Enter: stop and move on
% (Be careful otherwise may add unwanted point :()
% [Unfortunately no rollback option is provided internally in this function, so please be
% as possible in operations]
% Each finished frame will be required to be saved (or not) outside this function, so if
% error encountered, the user could change the starting frame and rerun the main function
%
% Zihan Zhang, Courant Institute
% Updated May 2023
%--------------------------------------------------------------------------

function [thisx,thisy] = func_manual(pic,datax,datay,fr)
% INPUT
% pic - background image
% datax, datay - extracted datapoints using alphashape-algo
% OUTPUT
% thisx, thisy - edited datapoints

    close all
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %% Delete Points
    imshow(pic); hold on; title('Deleted Points');
    scatter(datax,datay,'green','filled');
    [delx,dely] = ginput2('PlotOpt'); 
    % Delete the closest point in data
    % To guarantee accuracy, mark carefully
    for i = 1:length(delx)
        delpoint = [delx(i),dely(i)];
        distdic = sqrt((datax-delpoint(1)).^2+(datay-delpoint(2)).^2);
        [~,sind] = min(distdic);
        datax(sind) = []; datay(sind) = [];
    end
    disp(['Deleted ', num2str(length(delx)), ' Points!']);
    close all

    %% Add Points
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow(pic); hold on; title('Added Points');
    scatter(datax,datay,'green','filled');
    [addx,addy] = ginput2('PlotOpt'); 
    thisx = [datax;addx];
    thisy = [datay;addy];
    disp(['Added ', num2str(length(addx)), ' Points!']);
    close all

    %% CheckPoint
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    imshow(pic); hold on; title('Final Result')
    scatter(thisx,thisy,'blue','filled');
    pause
    close all

    [thisx,thisy] = fr.angle_sort(thisx,thisy);
end





