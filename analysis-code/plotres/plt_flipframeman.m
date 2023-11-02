%--------------------------------------------------------------------------
% Plot frames with manuel operated ("twisted") datapoints 
% --> func_manual.m 
%
% Steven Zhang, Courant Institute
% Updated May 2023
%--------------------------------------------------------------------------
close all

foldername = '2022-06-30-c/';
actfolder = ['datas/flip-moments/',foldername];
S = dir(fullfile(actfolder,'**','*.mat'));
names = {S.name}; ll = length(names); 
nm = []; % register to numerics
for i = 1:ll
    nk = names(i); nr = nk{1}; 
    nm(end+1) = str2num(nr(1:end-4));
end
nm = sort(nm);
nm(end+1) = nm(1); nm(1) = []; % move the initial shape to last plot 
assert(mod(ll,2) == 1) % must have 2x+1 frames

figure('units','normalized','outerposition',[0 0 1 1]);
[ha,~] = tight_subplot(round(ll/2/2),2,[.01 .03],[.1 .01],[.01 .01]);

for ii = 1:ll/2
    axes(ha(ii));
    cnt = (ii-1)*2+1; 
    fframe = load([actfolder,num2str(nm(cnt)),'.mat']); pt1 = fframe.pt; 
    sframe = load([actfolder,num2str(nm(cnt+1)),'.mat']); pt2 = sframe.pt; 
    axis equal
    hold on
    scatter(pt1(:,1)-fframe.cx,pt1(:,2),'blue');
    xlim([-0.05 0.05])
    ylim([-0.1 -0.02])
    scatter(pt2(:,1)-sframe.cx,pt2(:,2),'red');
%     yline(wlineset(1),'LineStyle','--')
    legend([num2str(fframe.ts),'s'],[num2str(sframe.ts),'s'],'waterline');
end

% saveas(gcf,[foldername(1:end-1),'-fliprecord.jpg'])

% IC
figure()
iframe = load([actfolder,num2str(nm(end)),'.mat']); pt = iframe.pt; 
scatter(pt(:,1)-iframe.cx,pt(:,2),'green')
% yline(wlineset(1),'LineStyle','--')
title("Initial Shape after First Locking Down")
axis equal

% saveas(gcf,[foldername(1:end-1),'-IC.jpg'])








