function [thisx,thisy] = func_manual(pic,datax,datay)
    close all
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %% Delete Points
    imshow(pic); hold on; title('Deleted Points');
    scatter(datax,datay,'green','filled')
    [xdel,ydel] = ginput2('PlotOpt'); 
    
    for i = 1:length(xdel)
        delpoint = [xdel(i),ydel(i)];
        distdic = sqrt((datax-delpoint(1)).^2+(datay-delpoint(2)).^2);
        [~,sind] = min(distdic);
        datax(sind) = []; datay(sind) = [];
    end
    disp(["Deleted ", num2str(length(xdel)), " Points!"])

    %% Add Points
    imshow(pic); hold on; title('Added Points');
    [addx,addy] = ginput2('PlotOpt'); 
    thisx = [datax;addx];
    thisy = [datay;addy];

    
end