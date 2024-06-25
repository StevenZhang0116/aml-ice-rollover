 %--------------------------------------------------------------------------
% Measurement and evaluation of flip time
% In measuring the flipping direction, 1 = clockwise, 0 = counterclockwise
% Count the flip time after the iceberg becomes steady in position
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------


close all

clear diff

path = '../datas/data-rolltime/';

delete([path,'/*.mat']);

gray = [128 128 128]/255;
deepyellow = [0.9290 0.6940 0.1250];
deepred = [0.6350 .0780 .1840];
deepgreen = [.466 .674 .188];

% name1 = '2022-06-07-a';
% result = [4*60+53,7*60+25,9*60+48,12*60+26,15*60+51,18*60+44,19*60+45,21*60+41,23*60+29,24*60+54,28*60+2];
% save([path,name1],"result");
% 
% name2 = '2022-06-07-b';
% result = [10*60+32,12*60+50,14*60+30,16*60+19,18*60+58,20*60+42,21*60+19,23*60+36,25*60+22,27*60+40];
% save([path,name2],"result");
% 
% name3 = '2022-06-14-a';
% result = [1*60+23,3*60+28,5*60+1,6*60+49,8*60+34,10*60+37,12*60+41,...
%     14*60+34,15*60+59,18*60+2,20*60+48,22*60+28];
% save([path,name3],"result");

% name4 = '2022-06-14-b';
% result = [2*60+15,4*60+11,6*60+11,8*60+8,10*60+51,12*60+52,14*60+15,15*60+40,17*60+47,20*60+45];
% direction = [0,0,1,1,0,1,0,1,0,0];
% save([path,name4],'result','direction');
% 
% name5 = '2022-06-15-a';
% result = [4*60+43,6*60+54,9*60+8,11*60+30,14*60+53,17*60+32,20*60+54,23*60+39];
% direction = [0,1,1,1,1,0,1,0];
% save([path,name5],'result','direction');
% 
% name6 = '2022-06-17-a';
% result = [3*60+21,5*60+47,7*60+55,10*60+10,13*60+11,15*60+51,19*60+11,21*60+36,24*60+46,27*60+23];
% direction = [1,0,0,0,0,0,0,1,0,1];
% save([path,name6],'result','direction');
% 
% name7 = '2022-06-17-b';
% result = [1*60+17,2*60+19,4*60+10,6*60+17,8*60+40,11*60+44,15*60+14,16*60+37,18*60+30,21*60+24,23*60+20];
% direction = [0,1,1,0,0,1,0,0,1,0,1];
% save([path,name7],'result','direction');

% name8 = '2022-06-17-c';
% result = [2*60+4,3*60+56,6*60+39,9*60+21,11*60+8,13*60+48,15*60+46,17*60+20,20*60+2,22*60+49,25*60+42];
% save([path,name8],"result");

% name9 = '2022-06-22-a';
% result = [3*60+37,6*60+16,8*60+48,11*60+23,15*60+16,17*60+24,21*60+13,24*60+17,26*60+3];
% save([path,name9],"result");

%% Bobae
%Bobae: no freetimes...
% name10 = '2022-06-22-b'; %DO NOT USE
% result = [2*60+27,4*60+41,6*60+46,9*60+4,11*60+4,13*60+22,14*60+33,15*60+59,17*60+48,19*60+51,22*60+18];
% direction = [0,1,0,1,1,0,0,1,0,1,0];
% save([path,name10],'result','direction');

% name11 = '2022-06-23-a'; %DO NOT USE
% result = [4*60+36,6*60+56,8*60+59,11*60+25,15*60+1,16*60+43,20*60+21,23*60+24,25*60+49];
% direction = [1,0,0,0,0,1,0,1,1];
% freetime = 0*60+52;
% save([path,name11],'result','direction','freetime');

name12 = '2022-06-28-a'; %use, double checked by Bobae
result = [6*60+25,8*60+53,11*60+16,13*60+36,16*60+6,18*60+54,20*60+25,21*60+44,23*60+25,25*60+33];
direction = [1,0,0,1,0,1,0,1,0,1]; %1 = right, 0 = left
freetime = 0*60+52;
save([path,name12],'result','direction','freetime');

name13 = '2022-06-28-b'; %use
result = [3*60+48,6*60+26,7*60+26,9*60+33,12*60+2,13*60+55,15*60+9,17*60+33,20*60+45,23*60+0,24*60+55,27*60+15];
direction = [0,1,0,1,1,0,1,1,0,1,0,1];
freetime = 0*60+52;
save([path,name13],'result','direction','freetime');

name14 = '2022-06-28-c'; %use
result = [9*60+31,11*60+48,14*60+12,16*60+16,18*60+51,21*60+35,23*60+40,25*60+39,27*60+33];
direction = [0,1,1,0,1,0,1,0,1];
freetime = 3*60+52;
save([path,name14],'result','direction','freetime');

name15 = '2022-06-29-a'; %use
result = [3*60+2,5*60+12,7*60+51,10*60+1,12*60+10,15*60+28,17*60+27,19*60+3,21*60+25,24*60+39,27*60+13];
direction = [0,1,0,1,1,0,1,0,1,0,0];
freetime = 1*60+22;
save([path,name15],'result','direction','freetime');

name16 = '2022-06-29-b'; %use
result = [2*60+17,5*60+2,7*60+0,9*60+50,12*60+19,15*60+30,17*60+36,19*60+11,21*60+30,23*60+54,26*60+50];
direction = [1,0,1,0,0,1,0,0,1,1,1];
freetime = 0*60+52;
save([path,name16],'result','direction','freetime');

name17 = '2022-06-29-c'; %use
result = [2*60+30,5*60+22,8*60+0,10*60+40,14*60+38,17*60+57,22*60+48,26*60+23,29*60+4];
freetime = 0*60+52;
direction = [0,1,1,1,1,0,1,0,1];
save([path,name17],'result','direction','freetime');

%INSERT 2022-06-30-a -- Bobae 06/21/24
%used very beginning of flip for timestamps.
name18 = '2022-06-30-a';
result = [2*60+6, 4*60+16, 6*60+4, 8*60+41, 11*60+26, 13*60+31, 16*60+16, 17*60+43, 18*60+47, 20*60+18, 22*60+26, 25*60+18];
freetime = 0*60+52; %not using this information anyway
save([path,name18],'result','direction','freetime')

% not that good one
name19 = '2022-06-30-b'; %use
result = [5*60+21,7*60+55,9*60+57,12*60+23,14*60+56,18*60+4,19*60+17,20*60+34,23*60+16,25*60+49];
direction = [1,0,0,1,0,1,0,0,0,0];
freetime = 3*60+32;
save([path,name19],'result','direction','freetime');

name20 = '2022-06-30-c'; %use
result = [1*60+51,4*60+27,6*60+34,9*60+29,12*60+26,14*60+52,16*60+58,19*60+2,22*60+13,24*60+38,26*60+35];
direction = [1,0,0,0,1,0,1,1,0,1,1];
freetime = 1*60+01;
save([path,name20],'result','direction','freetime');

name21 = '2022-07-02-b'; %use
result = [2*60+40,5*60+1,7*60+11,9*60+21,11*60+16,14*60+35,16*60+40,17*60+44,21*60+9,24*60+25,26*60+35];
direction = [0,1,0,1,1,0,1,1,0,1,1];
freetime = 0*60+52;
save([path,name21],'result','direction','freetime');

%INSERT 2022-07-05-a: Bobae 06/21/24 
%Notes: between flips 4-5, the ice hits the trap and doesn't flip
%successfully.
%Maybe delete last data point because quite small.
name22 = '2022-07-05-a';
% result = [2*60+6, 4*60+17, 6*60+3, 8*60+42, 11*60+27, 13*60+32, 16*60+16, 17*60+43, 18*60+47, 20*60+17, 22*60+25, 25*60+18];
result = [2*60+25, 4*60+22, 5*60+22, 7*60+13, 9*60+35, 11*60+37, 13*60+2, 14*60+53, 18*60+38, 20*60+51, 22*60+58, 25*60+4];
% direction = [0,1,0,1,1,0,1,1,0,1,1]; 
% freetime = 0*60+52;
save([path,name22],'result');%,'direction','freetime');

% IF USE THESE: MAKE SURE NAME NUMBERS ARE CHANGED FROM INSERTING OTHER
% INFORMATION.
% name21 = '2022-07-11-a-fixed';
% result = [2*60+3,4*60+37,6*60+39,8*60+58,12*60+14,14*60+27,16*60+59,19*60+50,22*60+2,23*60+52,25*60+12];
% direction = [0,1,1,1,1,0,1,0,1,0,1];
% save([path,name21],'result','direction');
% 
% name22 = '2022-07-08-c-fixed';
% result = [2*60+23,4*60+53,6*60+57,9*60+18,12*60+46,14*60+43,18*60+32,21*60+5,22*60+28,24*60+10,26*60+14];
% direction = [1,0,0,0,0,1,0,1,0,1,0];
% save([path,name22],'result','direction');
% 
% name23 = '2022-07-08-b-fixed';
% result = [1*60,3*60+13,5*60+32,8*60+10,10*60+4,12*60+38,14*60+9,15*60+31,17*60+41,19*60+39,21*60+40,24*60+52];
% direction = [0,1,0,1,1,0,0,1,0,1,1,0];
% save([path,name23],'result','direction');
% 
% name24 = '2022-07-08-a-fixed';
% result = [3*60+14,5*60+30,7*60+46,9*60+24,11*60+12,14*60+13,15*60+57,18*60+15,21*60+17,24*60+0];
% direction = [1,0,1,0,0,1,0,1,0,1];
% save([path,name24],'result','direction');

name25 = '2022-07-13-a'; %use
result = [1*60+38,3*60+58,6*60+5,8*60+43,10*60+29,12*60+46,14*60+7,15*60+50,17*60+42,19*60+13,19*60+52,21*60+2,22*60+55,...
    25*60+9,26*60+51, ];
direction = [1,0,1,0,0,1,0,1,0,1];
freetime = 0*60+52;
save([path,name25],'result','direction','freetime');

name26 = 'simulation';
% result = [5*60+47,7*60+0,8*60+50,10*60+20,12*60+41,14*60+28,16*60+40,18*60+22,19*60+55,...
%     21*60+20,22*60+34,23*60+55];
result = [3*60+35,5*60+10,7*60+30,10*60+5,12*60+40,15*60+15,17*60+45,21*60+40,24*60+45,26*60+5,28*60+25];
save([path,name26],'result','direction');

S = dir(fullfile(path,'**','*.mat'));
names = {S.name};

figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 1:length(names)
    fliptime = load([path,char(names(i))]).result;
    ind = 1:length(fliptime);
    hold on
    if i == length(names)
        plot(fliptime/60,ind,'-o','LineWidth',2,'Color',deepred)
    else
        plot(fliptime/60,ind,'-o','Color',gray)
    end
    % Checking that new inserted lines aren't too crazy:
    % if string(char(names(i))) == '2022-07-05-a.mat' 
    %     plot(fliptime/60,ind,'-o','LineWidth',5,'Color','r')
    %     i
    %     fliptime
    % end
    % if string(char(names(i))) == '2022-06-30-a.mat' 
    %     plot(fliptime/60,ind,'-o','LineWidth',5,'Color','k')
    %     i
    %     fliptime
    % end
end

newnames = cell(1,length(names));
for i = 1:length(names)
    nn = char(names(i));
    newnames{i} = nn(1:end-4);
end

hold off
legend(newnames,'Location','northwest','FontSize',12) %Original
legend({'experiment','','','','','','','','','','','','simulation'},'Interpreter','Latex','Location','northwest','FontSize',12)
title('Flip Time Plot')
ylabel('Cumulative Number of Flips','Interpreter','Latex','FontSize',16)
xlabel('Time, t(min)','Interpreter','Latex','FontSize',16)
title('Cumulative number of flips from experiments and simulation','FontSize',16,'Interpreter','Latex') 
ylim([0,16])















