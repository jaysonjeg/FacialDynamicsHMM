function [cfss7,frqvalues,n_ylabels,cfss6] = process_cfss(cfss,frq,usefreqbins,cutoff_upper,cutoff_lower)
%Prepare time-freq data for re-entry into HMM.

nframes=size(cfss,2);
nsubs=size(cfss,3);
nAUs=size(cfss,4);

nend=find(flipud(frq)>cutoff_upper,1); %default 1
nstart=find(flipud(frq)>cutoff_lower,1);

if usefreqbins
    n_ylabels=10; %how many frequency bins
    n_ylabels=n_ylabels-1;
    indices=round(linspace(nstart,nend,n_ylabels+1));
else
    indices=nstart:nend;  
end

%find mean frequencies of each frequency band
frq_flipped=flipud(frq);
frqvalues=frq_flipped(indices);
meanvalues=mean([frqvalues(1:end-1),frqvalues(2:end)],2); 
n_ylabels=length(frqvalues);

cfss2=cfss(fliplr(1:length(frq)),:,:,:);  %flip frequency to go from small to big
cfss3=cfss2(indices,:,:,:); clear cfss2; %include only relevant indices
cfss4=reshape(cfss3,[],n_ylabels,nframes,size(cfss3,3),nAUs); clear cfss3; %dim2 groups of size dim1
cfss5=squeeze(mean(cfss4,1)); clear cfss4; %average within each frequency band. Relevant if usefreqbins==true
cfss6=permute(cfss5,[2,3,1,4]); clear cfss5; %nframes * nsubs * nfreqbins * nAUs
cfss7=reshape(cfss6,nframes,[],n_ylabels*nAUs); %nframes * nsubs * datapoints(nfreqbins(fine)*nAUs(coarse))
%{
%Any subjects who have any features (freqbins*nAUs) with zero variance,
%will get first element replaced by 0.01. This is not good practice, and
%commented out by default.
ndatapoints=n_ylabels*nAUs;
for i=1:size(cfss7,2)
    for j=1:ndatapoints
        if var(cfss7(:,i,j))==0
            cfss7(1,i,j)=0.01;
        end
    end
end
%}
playtone();
end

