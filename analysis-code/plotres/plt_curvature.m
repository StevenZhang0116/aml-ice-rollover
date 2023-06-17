%--------------------------------------------------------------------------
% Measurement and evaluation of curvature
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------

close all

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on

path = './data-curvature/';
subpath = char("normal/");
S = dir(fullfile([path,subpath],'**','*.mat'));
names = {S.name};
nn = length(names);
alllength = cell(1,nn);
alllength2 = cell(1,nn);
alltime = cell(1,nn);
lengthrecord = [];
lengthrecord2 = [];
flipset = cell(1,nn);

hold on
movage = 10;
curv_res = cell(1,nn+1);
% load the data
for i=1:nn
    matname = char(names(i));
    acc = load([path,subpath,matname]).adjustedcurvatureSet;
    acc2 = load([path,subpath,matname]).curvature3Set;
    timeLab = load([path,subpath,matname]).timeLab;
    % delete the timestamp before the first flip
    ftpath = './data-rolltime/';
    fltime = load([ftpath,matname]).result;

    % separate by flip [flip-wise average]
    newfltime = [0,fltime,0];
    tmp1 = []; tmp2 = [];
    prevind = 1;
    for j=2:length(newfltime)-1
        thisflip = newfltime(j);
        res = find(timeLab == thisflip);
        if length(res) > 0
            tmp1(end+1) = j-1;
            meanpoly = mean(acc(prevind:res));
            prevind = res;
            tmp2(end+1) = meanpoly;
        end
    end
    flipset{i} = [tmp1;tmp2];
    lengthrecord2(end+1) = length(tmp1);

    % first flip
    frflip = fltime(1);
    % truncate data/time before that 
    res = find(timeLab == frflip);
    acc = acc(res+1:end);
    acc2 = acc2(res+1:end);
    timeLab = timeLab(res+1:end);
    % standardize; set first flip time as 0
    timeLab = timeLab - timeLab(1);
    % moving average
    acc = movmean(acc,movage);
    acc2 = movmean(acc2,movage);
    timeLab = movmean(timeLab,movage);
    alllength{1,i} = acc;
    alllength2{1,i} = acc2;
    alltime{1,i} = timeLab;
    lengthrecord(end+1) = length(acc);
    plot(timeLab/60, acc,'LineStyle','--');
    % save trail result
    curv_res{i} = [timeLab/60;acc];
end

[longe,ii] = max(lengthrecord);
longtimelab = alltime{1,ii};
% calculate the mean curvature over all trails
longc = [];
for i=1:nn
    thecurve = alllength{1,i};
    fillnan = NaN(1,longe-length(thecurve));
    thecurve = [thecurve,fillnan];
    longc = [longc;thecurve];
end

plot(longtimelab/60,nanmean(longc),'-o','LineWidth',3)
% save mean result
curv_res{length(curv_res)} = [longtimelab/60;nanmean(longc)];
xlabel('Time, t(min)','FontSize',16)
ylabel('Polygonality','FontSize',16)
title(['Trailed-Average Polygonality using Curvature Method with Moving Average = ', ...
    num2str(movage),'s'],'FontSize',16)
xlim([0,1000/60])

save('curv_res.mat','curv_res');

% plot flip-wise average polygonality
[longe,ii] = max(lengthrecord2);
longtimelab = 1:1:longe;
longc = [];
for i=1:nn
    thecurve = flipset{i}(2,:);
    fillnan = NaN(1,longe-length(thecurve));
    thecurve = [thecurve,fillnan];
    longc = [longc;thecurve];
end

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
for i = 1:length(flipset)
    onedata = flipset{i};
    plot(onedata(1,:),onedata(2,:),'-o');
end
plot(longtimelab,nanmean(longc),'-o','LineWidth',3)
xlabel('Flip Time Count','FontSize',16)
ylabel('Average Polygonality','FontSize',16)
title('Trailed-Average Polygonality using Curvature Method categorized by Flipover Times', ...
    'FontSize',16)
xlim([1 8])




