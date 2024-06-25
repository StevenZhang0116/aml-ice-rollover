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
interval = 30;
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% for i = 1:length(names)
%     matname = char(names(i));
%     hold on
%     if strcmp(matname,'simulation.mat') == 0 %plot experiment
%         times = load([path,matname]).timeLab;
%         areas = load([path,matname]).areaSet;
%         times = times(1:interval:1050); 
%         areas = areas(1:interval:1050);
%         areas = areas/areas(1);
%         xx = polyfit(times,areas.^(5/8),1);
%         tf = -xx(2)/xx(1);
%         times = times/tf;
%         plot(times,areas,'-o')
%     else %plot simulation
%         times = load([path,matname]).timeset;
%         times = times(1:interval:1200);
%         areas = load([path,matname]).aa0set;
%         areas = areas(1:interval:1200);
%         xx = polyfit(times,areas.^(5/8),1);
%         tf = -xx(2)/xx(1);
%         plot(times/tf,areas,'Marker','o','MarkerSize',10);
%     end
% end
% 
% newnames = cell(1,length(names));
% for i = 1:length(names)
%     nn = char(names(i));
%     newnames{i} = nn(1:end-4);
% end
% 
% tt = linspace(0,1,1000);
% plot(tt,(1-tt).^(8/5),'LineStyle','-.','LineWidth',2);
% hold off
% xlabel('Time, t/t_f','FontSize',16)
% ylabel('Area, A/A_0','FontSize',16)
% legend(newnames,'Location','northeast','FontSize',12)
% title('Calculated area from representative experiments and simulation','FontSize',16)
% ax = gca;
% ax.FontSize = 16; 
% 

%% Bobae

%% Creating plot with gray experimental lines
gray = [128 128 128]/255;
deepyellow = [0.9290 0.6940 0.1250];
deepred = [0.6350 .0780 .1840];
deepgreen = [.466 .674 .188];

num_expts = 0;
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 1:length(names)
    matname = char(names(i));
    hold on
    if strcmp(matname,'simulation.mat') == 0 %plot experiment
    % if i == 1 %for bobae to plot just one experiment
        times = load([path,matname]).timeLab;
        areas = load([path,matname]).areaSet;
        % length(areas);
        % times = times(1:interval:1050); %original
        % areas = areas(1:interval:1050); %original
        times = times(1:interval:end); %bobae
        areas = areas(1:interval:end); %bobae
        areas = areas/areas(1); %normalize to be A/A0
        xx = polyfit(times,areas.^(5/8),1);
        tf = -xx(2)/xx(1);
        % max(times)
        times = times/tf;
        plot(times,areas,'-o','Color',gray)
        % Compute average:
        if i == 1
            avg_areas = areas;
        else
            avg_areas = avg_areas + areas;
        end

        % avg_areas = avg_areas + areas;
        expt_times = times;
        num_expts = num_expts+1;
    elseif strcmp(matname,'2022-06-22-b.mat') == 1 %Bobae, removing extra experiment 2022-06-22-b
        num_expts = num_expts;
    elseif strcmp(matname,'2022-06-23-a.mat') == 1 %Bobae, removing extra experiment 2022-06-23-a
        num_expts = num_expts;
    else %plot simulation
    % if strcmp(matname,'simulation.mat') == 1 %bobae just to plot
    % simulation alone
        times = load([path,matname]).timeset;
        % times = times(1:interval:1200);  %original
        times = times(1:interval:end); %bobae
        areas = load([path,matname]).aa0set;
        % areas = areas(1:interval:1200); %original
        areas = areas(1:interval:end); %bobae
        xx = polyfit(times,areas.^(5/8),1);%xx = polyfit(times,areas.^(8/5),1);
        tf = -xx(2)/xx(1); %tf = max(times)
        sim_times = times/tf;
        plot(sim_times,areas,'Marker','.','MarkerSize',25,'Color',deepred,'Linewidth',2);
        sim_areas = areas;
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
xlabel('time, $t/t_f$','Interpreter','Latex','FontSize',16)
ylabel('area, $A/A_0$','Interpreter','Latex','FontSize',16)
% legend(newnames,'Location','northeast','FontSize',12) %Original
legend({'experiment','','','','','','','','','','','simulation','$(1-t/t_f)^{8/5}$'},'Interpreter','Latex','Location','northeast','FontSize',12)
% title('Calculated area from representative experiments and simulation','FontSize',16)
title('Calculated area from representative experiments and simulation','FontSize',16)
ax = gca;
ax.FontSize = 16; 


%% Plot average areas
figure('Renderer', 'painters', 'Position', [10 10 900 600])
% num_expts = length(names) - 1; %number of experiments
avg_areas = avg_areas / num_expts;
plot(expt_times,avg_areas,'o','MarkerSize',8,'Color',deepgreen,'Linewidth',2);
hold on
plot(sim_times, sim_areas,'o','MarkerSize',8,'Color',deepred,'Linewidth',2)
plot(tt,(1-tt).^(8/5),'LineStyle','--','LineWidth',2,'Color','k');
legend({'experimental average','simulation','$(1-t/t_f)^{8/5}$'},'Interpreter','Latex','Location','northeast','FontSize',12)
xlabel('time, $t/t_f$','Interpreter','Latex','FontSize',16)
ylabel('area, $A/A_0$','Interpreter','Latex','FontSize',16)
title('Calculated area from simulation and averaging over experiments','FontSize',16)
ax = gca;
ax.FontSize = 16; 