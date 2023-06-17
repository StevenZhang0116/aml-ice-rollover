%--------------------------------------------------------------------------
% Plot frames with manuel operated ("twisted") datapoints 
% --> func_manual.m 
%
% Steven Zhang, Courant Institute
% Updated May 2023
%--------------------------------------------------------------------------

foldername = '2022-06-28-b/';
actfolder = ['../datas/flip-moments/',foldername];
S = dir(fullfile(actfolder,'**','*.mat'));
names = {S.name}; ll = length(names); 
nm = []; % register to numerics
for i = 1:ll
    nk = names(i); nr = nk{1}; 
    nm(end+1) = str2num(nr(1:end-4));
end
nm = sort(nm);
assert(mod(ll,2) == 0) % must have 2x frames

[ha,~] = tight_subplot(round(ll/2),2,[.01 .03],[.1 .01],[.01 .01]);

for ii = 1:ll/2
    axes(ha(ii));
    cnt = (ii-1)*2+1; 
    fframe = load([actfolder,num2str(nm(cnt)),'.mat']); pt1 = fframe.pt; 
    sframe = load([actfolder,num2str(nm(cnt+1)),'.mat']); pt2 = sframe.pt; 
    hold on
    scatter(pt1(:,1)-fframe.cx,pt1(:,2)-fframe.cy,'blue');
    scatter(pt2(:,1)-sframe.cx,pt2(:,2)-sframe.cy,'red');
    legend([num2str(fframe.ts),'s'],[num2str(sframe.ts),'s']);
end