%--------------------------------------------------------------------------
% Plot Boundary Evolution
%
% Steven Zhang, Courant Institute
% Updated Jan 2023
%--------------------------------------------------------------------------

ff = load('../../../process-data/sample_boundary.mat');

timeLab = ff.timeLab;
bdcell = ff.interpolation_cell;

close all
figure('units','normalized','outerposition',[0 0 0.8 0.8]);
breakp = 120; % s
totalframe = floor(timeLab(end)/breakp) - floor(timeLab(1)/breakp);

aa = 'alpha';

if totalframe > 12
    totalframe = 12;
end

wl = -0.02995; % water line
xx = 0:0.00001:0.1;
wlpt = [];
for i = 1:length(xx)
    wlpt(end+1,1) = xx(i);
    wlpt(end,2) = wl;
end
 
for i = 1:length(timeLab)
     if mod(timeLab(i),breakp) == 0
         time = timeLab(i);
         kk = timeLab(i)/breakp;
         if kk > 12
             break
         end
         subplot(2,totalframe/2,kk);
         dd = bdcell{i};
         xx = dd(:,1);
         yy = dd(:,2);
         ss = polyshape(xx,yy);
         hold on
         plot(xx,yy,'red','LineStyle','-','LineWidth',2)
         axis equal
         xlim([0 0.1])
         ylim([-0.1 0])
         title([num2str(time),'s'])
         [in,on] = inpolygon(wlpt(:,1),wlpt(:,2),xx,yy);
         plot(wlpt(~in,1),wlpt(~in,2),'o','MarkerSize',2,'Color','green')
         plot(ss)
         hold off
     end
end

sgtitle('Boundary Evolution Showcase','FontSize',14)
saveas(gcf,'boundary_trend.png')



