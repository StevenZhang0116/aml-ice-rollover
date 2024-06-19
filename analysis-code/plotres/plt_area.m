%--------------------------------------------------------------------------
% Measurement and evaluation of area
%
% Steven Zhang, Courant Institute
% Updated Jan 2023
%--------------------------------------------------------------------------
close all

path = '../datas/data-area/';
S = dir(fullfile(path,'**','*.mat'));
names = {S.name};
interval = 50;
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 1:length(names)
    matname = char(names(i));
    hold on
    if strcmp(matname,'simulation.mat') == 0
        times = load([path,matname]).timeLab;
        areas = load([path,matname]).areaSet;
        times = times(1:interval:1050); 
        areas = areas(1:interval:1050);
        areas = areas/areas(1);
        xx = polyfit(times,areas.^(5/8),1);
        tf = -xx(2)/xx(1);
        times = times/tf;
        plot(times,areas,'-o')
    else
        times = load([path,matname]).timeset;
        times = times(1:interval:1200);
        areas = load([path,matname]).aa0set;
        areas = areas(1:interval:1200);
        xx = polyfit(times,areas.^(5/8),1);
        tf = -xx(2)/xx(1);
        plot(times/tf,areas,'Marker','o','MarkerSize',10);
    end
end

newnames = cell(1,length(names));
for i = 1:length(names)
    nn = char(names(i));
    newnames{i} = nn(1:end-4);
end

tt = linspace(0,1,1000);
plot(tt,(1-tt).^(8/5),'LineStyle','-.','LineWidth',2);
hold off
xlabel('Time, t/t_f','FontSize',16)
ylabel('Area, A/A_0','FontSize',16)
legend(newnames,'Location','northeast','FontSize',12)
title('Calculated area from representative experiments and simulation','FontSize',16)
ax = gca;
ax.FontSize = 16; 


