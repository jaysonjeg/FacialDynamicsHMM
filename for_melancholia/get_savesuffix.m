function [ii,savesuffix] = get_savesuffix(hmm,ID,maxFO,fehist, savefolder)

h1=hmm{ID};
files=dir(savefolder);
names={files(3:end).name};
ii=max(cellfun(@(x) get_number(x),names))+1; %get next figure number
savesuffix=sprintf('%i_maxF%.2f_fe%.3d.jpg',h1.K,mean(maxFO{ID}),fehist{ID}(end));
end

function y=get_number(x)
underscores=find(x=='_');
first_underscore=underscores(1);
indices=1:first_underscore-2;
y=eval(x(indices));
end
