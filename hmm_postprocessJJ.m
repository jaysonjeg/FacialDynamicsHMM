function [nframesdown,ntotalframesdown,Tdown,vpathx,Gammax,Gammac,Masks_subject]=hmm_postprocessJJ(hmm,vpath,ID,nsubs,Gamma)

%Variables with 'down' are adjusted for downsampling
nframesdown=height(vpath{ID})/nsubs; %frames per person
ntotalframesdown=height(Gamma{ID}); %total n of frames with subject-concatenation
Tdown=ones(nsubs,1)*nframesdown; 

%Prepare Viterbi paths and Gamma
vpathx=reshape(vpath{ID},nframesdown,[]);
Gammax=reshape(Gamma{ID},nframesdown,[],hmm{ID}.K);
Gammac=squeeze(num2cell(permute(Gammax,[1,3,2]),[1,2]));

%Mask for each individual subject
temp=reshape(1:ntotalframesdown,nframesdown,[]);
Masks_subject={}; 
for i=1:size(temp,2)
    Masks_subject=[Masks_subject,temp(:,i)'];
end

end

