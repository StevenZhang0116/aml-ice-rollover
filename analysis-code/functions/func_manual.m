function [thisx,thisy] = func_manual(pic,datax,datay)
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