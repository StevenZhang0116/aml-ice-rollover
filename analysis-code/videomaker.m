video = VideoWriter('bdtrackingsample.avi');
open(video); 

path = './bd-tracking/';
S = dir(fullfile(path,'**','*.jpeg'));
names = {S.name};
N = length(names);

for ii=1:N 
  I = imread([path,num2str(ii),'.jpeg']); 
  disp(ii)
  writeVideo(video,I); 
end
close(video); 