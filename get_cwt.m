function [cfssraw,frq,coi] = get_cwt(cube,Fs,use_cwtfilterbank)
%{
Get continuous wavelet transform (cwt)

Inputs:
cube: nframes x nAUs x nsubjects
Fs: sampling frequency
use_cwtfilterbank: True (default) does with cwtfilterbank. False does with calculate_cwt. 

Outputs:
cfssraw: output array freqbins * nframes * nsubs * nAUs
frq: frequency bins
coi: cone of influence
frq and coi are in frequency space
%}

n=size(cube,1);

if use_cwtfilterbank
    fb = cwtfilterbank('SignalLength',n,... 
        'SamplingFrequency',Fs,...
        'VoicesPerOctave',12);
end

%this section is just to get number of freq bins so we can preallocate
%cfssraw (improved speed)
y=cube(:,1,1); 
if use_cwtfilterbank
    [wave,frq,coi,~] = wt(fb,y);
else
    [wave,frq,scale,coi] = calculate_cwt(y,n,Fs);
end
cfssraw=NaN(size(wave,1),size(cube,1),size(cube,3),size(cube,2));

tic;
for nsub=1:size(cube,3)
    temp={};
    for nAU=1:size(cube,2)
        y=cube(:,nAU,nsub);
        
        if use_cwtfilterbank
            [wave,frq,coi,~] = wt(fb,y);
        else
            [wave,frq,scale,coi] = calculate_cwt(y,n,Fs);
        end
        cfssraw(:,:,nsub,nAU)=wave;
    end
end
toc;

if use_cwtfilterbank
    frq=flip(frq)';
    coi=coi';
end

end

