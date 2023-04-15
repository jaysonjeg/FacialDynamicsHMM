%{
Starting from array a_aus, get continuous wavelet transform using
get_cwt.m, then run HMM using HMM-MAR toolbox
This code requires variable a_aus in workspace, which is an array: ntimepoints x nAUs
(17) * nSubjects
%}

%Prepare data
format compact
aulist={'01','02','04','05','06','07','09','10','12','14','15','17','20','23','25','26','45'};
aunames=containers.Map({1,2,3,4,5,6,7,...
    8,9,10,11,12,13,...
    14,15,16,17},...
    {'InnerBrowRaiser','OuterBrowRaiser','BrowLowerer','UpperLidRaiser','CheekRaiser','LidTightener','NoseWrinkler'...
    'UpperLipRaiser','LipCornerPuller','Dimpler','LipCornerDepressor','ChinRaiser','LipStretcher'...
    'LipTightener','LipParts','JawDrop','Blink'});
these_AUs= [1,3,5:16];  %1:16; %don't include 17th AU which is AU45 blink, and AU02/05 which have constant 0 value in some subs in melancholia dataset
event_times=[40,70,96,125,199]; %Stimulus start times (sec) of each video in DISFA dataset
Fs = 10; %sampling freq of DISFA dataset

cube=a_aus(:,these_AUs,:); clear a_aus;
nframes=size(cube,1);
nAUs=size(cube,2);
nsubs=size(cube,3);

t=(0:size(cube,1)-1)/Fs; %list of timestamps
%%
%Optional: Supplementary Figure 1. Plot mean time series
cubem=squeeze(mean(cube,3));
figure;
for i=1:length(these_AUs)
    nplotsheight=ceil(sqrt(length(these_AUs)));
    subplot(nplotsheight,nplotsheight,i);
    nAU=these_AUs(i);
    plot(t,cubem(:,i),'b');
    title(sprintf('AU%s: %s',aulist{nAU},aunames(nAU))); ylabel('Intensity'); xlabel('Time (s)');
    draw_eventlinesInPlot(event_times,1)
    xlim([0,max(t)]);
end
clear cubem nplotsheight nAU i;

%%
%For each subject, perform continuous wavelet transform. 
maxframe=800; %sett to 800 to just include the comedy stimulus, so that we can make FIgure 1
[cfssraw,frq,coi] = get_cwt(cube(1:maxframe,:,:),Fs,true); %cfssraw is array:(81 frequency bins x 800 time points x 27 subjects x 14 action units)
[period,coi]=get_periods(frq,coi);
cfss=get_cwt_mag(cfssraw,cube(1:maxframe,:,:)); %Get just magnitude of CWT: nfreqs x ntimepoints x nsubs x nAUs
cfssm=squeeze(mean(cfss,3)); %Mean across subjects of CWT amplitude: 81 freq bins x 800 time points x 14 AUs
clear cfssraw

%%
%Optional: Diagnostic test to make sure pipeline reproducible
cubem=squeeze(mean(cube,3));
format long
cubem(1,1)
cfss(10,10,10,10)
format short
%{
Should be
0.376296296296296
0.008107829102573
%}

%%
%Optional but need this for Figure 1. Get continuous wavelet transform of group mean
cubem=squeeze(mean(cube(1:maxframe,:,:),3)); %Group mean time series
[mcfssraw,~,~]=get_cwt(cubem,Fs,true); %CWT of group mean time series
mcfss=get_cwt_mag(mcfssraw,cubem); %Amplitude of CWT of group mean time series: nfreqs x ntimepints x 1 x nAUs
%%
%Optional but need this for Figure 1. Mean of CWT minus CWT of mean - t-test for significant
%difference

nAU=find(strcmp(aulist(these_AUs),'12'));
diff=cfssm-mcfss;  %nfreqs x ntimepoints x nAUs
diff=diff(:,:,nAU); %nfreqs x ntimepoints

cfssn=cfss-mcfss; %Normalize, by subtracting cwt of group mean time series
cfsst=cfssm; %just to initialise size. Will store p-values of t-test for difference from group mean time series' value
tic
for i=1:size(cfsst,1)
    for j=1:size(cfsst,2)
        for k=nAU %1:size(cfsst,3)
            [h,p]=ttest(squeeze(cfssn(i,j,:,k))); %2-sided t-test for difference in group mean of cfss, from the value in mcfss
            cfsst(i,j,k)=p; %nfreqs x ntimepoints x nsubs
        end
    end
end
toc
pvalues=cfsst(:,1:maxframe,nAU); %nfreqs x ntimepoints
fdr=reshape(mafdr(pvalues(:),'BHFDR',true),size(pvalues)); %FDR correction
sig=(fdr<0.05); 
clear cfsst cfssn h p i j k mcfssraw

%%
%{
Make Figure 1
%}

n=5; %how many random sample participants
cmap_sig=turbo;
this_coi=coi(1:maxframe);
figure('Position',[100,100,1000,800]);
p=panel();

p.pack({0.3,0.7});
p(1).pack('h',{0.45,0.55});
p(2).pack('h',{0.98,[]});
p(2,1).pack({0.35,0.65});
p(2,1,1).pack('h',n);
p(2,1,2).pack('h',3);

ppics=p(1,1); %name sub-panels to make it easier to call them
p(1,2).pack({0.9,0.1});
plines=p(1,2,1);
psamples=p(2,1,1);
pmean1=p(2,1,2,1); 
pmean2=p(2,1,2,2);
pmean3=p(2,1,2,3);

%Panel A: Show pictures of participant SN003 face from 8 and 23 seconds
ppics.pack('h',2);
ppics(1).select();
temp=imread('C:\Users\Jayson\Google Drive\PhD\Project_Melancholia\MelancholiaHMM Paper submission\Nature Human Behaviour\Figures\SN003_0008_copy.png');
imshow(temp); ylim([0,maxframe]);
ppics(2).select();
temp=imread('C:\Users\Jayson\Google Drive\PhD\Project_Melancholia\MelancholiaHMM Paper submission\Nature Human Behaviour\Figures\SN003_0023_copy.png');
imshow(temp); ylim([0,maxframe]);

%Panel B: Plot raw AU time series, and group mean time series
plines.select();
xx=brewermap(2,'RdBu'); c2=xx(1,:); c1=xx(2,:);
rng(1); %set random seed
inds=randsample(size(cfss,3),n);
inds(1)=3; %always include SN003
hold on
for i=1:length(inds)
    values=cube(1:maxframe,nAU,inds(i));
    if i==1 %bold the plot of participant SN003
        plot(t(1:maxframe),smooth(values,10),'Color',c1,'LineWidth',2); 
    else
        plot(t(1:maxframe),smooth(values,10),'Color',c1);
    end
end
plot(t(1:maxframe),cubem(1:maxframe,nAU),'Color',c2,'LineWidth',2);
hold off
ylabel('Intensity'); xlabel('Time (s)'); axis tight

%Panel C: Plot n samples' individual CWTs
clims=[];
for i=1:n 
    psamples(i).select();
    plot_cwt(cmap_sig,squeeze(cfss(:,:,inds(i),nAU)),period,Fs,this_coi);
    temp=get(gca,'cLim'); clims(end+1)=temp(2); %get upper limit of each colorbar
end
%Set all clims to have the same upper limit
clims=max(clims); 
for i=1:n 
    psamples(i).select();
    set(gca,'cLim',[0,clims]);
end
colorbar('Position',[.93,.51,.015,.15]);

%Panel D: Plot mean of CWT
pmean1.select(); 
clims=[];
plot_cwt(cmap_sig,squeeze(cfssm(:,:,nAU)),period,Fs,this_coi);
clims(end+1)=max(get(gca,'cLim'));

%Panel E: Plot CWT of mean
pmean2.select(); 
plot_cwt(cmap_sig,squeeze(mcfss(:,:,:,nAU)),period,Fs,this_coi);
clims(end+1)=max(get(gca,'cLim'));
clims=max(clims);
pmean1.select(); set(gca,'cLim',[0,clims]); %Make clim common
pmean2.select(); set(gca,'cLim',[0,clims]);

%Panel F: Plot difference (D-E)
[ax1,ax2]=plot_cwt(cmap_sig,diff,period,Fs,this_coi,sig); 
pmean3.select([ax1 ax2]); %as in demopanelG
colorbar('Position',[.93,.17,.015,.15]);

%Place x/ylabels
p(2,1,1,1).select();
xlabel('Time (s)'); ylabel('Frequency (Hz)');
pmean1.select();
xlabel('Time (s)'); ylabel('Frequency (Hz)');

%Move margins
plines.marginleft=5;
p.de.marginright=10;
p.margintop=5; p.marginbottom=15; p.marginright=5; p.marginleft=15;
p(1,1).de.margin=0;
p(1,1).marginright=5;
p(1).marginbottom=0;
p(1,2).marginleft=0;
p(1,2,1).marginbottom=0;
p.fontname='TimesNewRoman';

%Label each panel
fontsize=14; 
annotation('textbox',[0.013,0.983,0,0],'String','a','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.45,0.983,0,0],'String','b','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.013,0.72,0,0],'String','c','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.013,0.465,0,0],'String','d','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.32,0.465,0,0],'String','e','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.63,0.465,0,0],'String','f','FitBoxToText','on','LineStyle','none','FontSize',fontsize);

%Arrows in line plot corresponding to first two face pictures
plines.select()
temp=[8,23];
for i=1:length(temp)
    [xaf,yaf]=ds2nfu(temp(i),0);
    annotation('arrow',[xaf,xaf],[yaf-0.04,yaf]);
end

%Red arrows for high freq activity in CWTmean plot
temp=[130,252,353,518,628,735]; 
temp=[130,252,518,735];
temp=temp/Fs;
for i=1:length(temp)
    plines.select();
    [xaf,yaf]=ds2nfu(temp(i),0);
    annotation('arrow',[xaf,xaf],[yaf-0.04,yaf],'Color',c2);
    pmean1.select()
    [xaf,yaf]=ds2nfu(temp(i),-3.1);
    annotation('arrow',[xaf,xaf],[yaf-0.04,yaf],'Color',c2);
    pmean2.select()
    [xaf,yaf]=ds2nfu(temp(i),-3.1);
    annotation('arrow',[xaf,xaf],[yaf-0.04,yaf],'Color',c2);
    pmean3.select()
    [xaf,yaf]=ds2nfu(temp(i),-3.1);
    annotation('arrow',[xaf,xaf],[yaf-0.04,yaf],'Color',c2);
end
clear pmean1 pmean2 pmean3 plines ppics

%%
%Get CWT of individual participants (for entire video). Not just comedy
[cfssraw,frq,coi] = get_cwt(cube,Fs,true); %Get CWT of individual pts
[period,coi]=get_periods(frq,coi);
cfss=get_cwt_mag(cfssraw,cube); %Get just magnitude of CWT: nfreqs x ntimepoints x nsubs x nAUs
cfssm=squeeze(mean(cfss,3)); %Mean across subjects of CWT amplitude: nfreqs x ntimepoints x nAUs
clear cfssraw
playtone();

%%
%Supplementary Figure 2. Show CWT plot.
figure; 
for nAU=1:length(these_AUs)
    subplot(4,4,nAU);  
    values=cfssm(:,:,nAU);
    plot_cwt(colormap('turbo'),values,period,Fs,coi);
    title(sprintf('AU%s: %s',aulist{these_AUs(nAU)},aunames(these_AUs(nAU)))); 
    xlabel('Time (s)'); ylabel('Frequency (Hz)'); colorbar;
    for i=1:length(event_times)
        xline(event_times(i));
    end
end

%%
%Optional: Check if any subjects/AUs have zero variance. These subjects/AUs will produce error with
%HMM
allvars=squeeze(var(cube,1))';
imagesc(allvars==0)
%%
%Prepare time-freq data for entry into HMM
clear cfssm cube
cfss=single(cfss);
usefreqbins=false; %true to use freqbins rather than whole cwt
cutoff_upper=4; %Upper cutoff frequency: default 5hz
cutoff_lower=0; %Lower cutoff frequency: default 0hz
[cfss7,frqvalues,n_ylabels,~] = process_cfss(cfss,frq,usefreqbins,cutoff_upper,cutoff_lower);
cfss8=reshape(cfss7,[],size(cfss7,3)); %ntimepoints(nframes(fine)*nsubs(coarse)) * datapoints
clear cfss7 cfss; 
playtone();
%%
%{
Prepare for HMM (using HMM-MAR toolbox)
%}
X=cfss8; clear cfss8;
T=repmat({[nframes]},nsubs,1); % length of data for each session
T=cell2mat(T);
Tsubject=T;

%%
%Optional: Supplementary Figure 3 - Plot relationship between nstates and
%free energy. Takes ~10mins
number_states=[2,4,8,16,32];
[windowsize_sec,smwin,configurations]=get_defaults_HMM(number_states);
for i=1:length(configurations)
    configurations{i}.Fs=Fs;
end
[hmm,Gamma,Gammaup,Xi,vpath,fehist,maxFO,newframerate] = runHMM_JJ1(configurations,X,T);
felast=cellfun(@(x) x(end),fehist);

figure; scatter(number_states,felast); xlabel('Number of states'); ylabel('Free energy');
set(gca,'XScale','log');

xticks([2,4,8,16,32,64]); ylim([1.7e6,2.2e6]);
%% Run HMM and perform statistical testing on the HMM results vs the conditions
[windowsize_sec,smwin,configurations]=get_defaults_HMM(8); %Get default configuration parameters
for i=1:length(configurations)
    configurations{i}.Fs=Fs;
end
[hmm,Gamma,Gammaup,Xi,vpath,fehist,maxFO,newframerate] = runHMM_JJ1(configurations,X,T);
%%
%{
Test for consistency
0.040378034226784   0.959540889968395   0.000001046779765   0.000000130682004   0.000077989429850   0.000000000046996   0.000000000000000
0.000001908866206
%}
format long
Gamma{1}(1000,:)
format short

%%
%Re-order states to put happy ones together, etc 
ord=[6,8, 2,4,1,7, 3,5]; %change old index 8 to second (new index 2), etc
ID=1; %We only had one HMM configuration, so ID=1
h1=hmm{ID};
vpath2=vpath;
for i=1:length(ord)
    vpath2{1}(vpath{1}==ord(i))=i;
end
vpath=vpath2; clear vpath2;
Gamma{ID}=Gamma{ID}(:,ord);
Gammaup{ID}=Gammaup{ID}(:,ord);
hmm{ID}.P=hmm{ID}.P(ord,ord);
hmm{ID}.Pi=hmm{ID}.Pi(ord);
hmm{ID}.Dir_alpha=hmm{ID}.Dir_alpha(ord);
hmm{ID}.Dir2d_alpha=hmm{ID}.Dir2d_alpha(ord,ord);
hmm{ID}.state=hmm{ID}.state(ord);

%%
%HMM post-processing including new variables to account for downsampling to 10Hz, reshape
%vpath/Gamma, and get individual subject masks
[nframesdown,~,~,vpathx,Gammax,~,~]=hmm_postprocessJJ(hmm,vpath,ID,nsubs,Gamma);
clear ntotalframesdown Tdown Gammac ii savesuffix Xi reorder
tic
windowsize=round(windowsize_sec*newframerate); %for most common Viterbi path and consistency plot, how many frames to average/smooth over
[vconsistency,vmostcommon,vconsistencies] = get_consistency_measures(windowsize,vpathx,hmm{ID}.K,true);
toc
playtone();
%[FO,ntrials] = getFractionalOccupancy (Gammaup{ID},T,configurations{ID}); %To get fractional occupancy of each state

%%
%Make 1000 surrogate Viterbi paths with each subject randomly circshifted.
%Need this for Figure 2, Panel D, grey shading
vpathxnull={};
n=1000;
for i=1:n
    i
    shifts=randi(size(vpathx,1),size(vpathx,2),1);
    temp=vpathx;
    for i=1:size(vpathx,2)
        temp(:,i)=circshift(vpathx(:,i),shifts(i));
    end
    vpathxnull{end+1}=temp;
end

%Find 5 and 95%ile of consistency
nullconsistencies=zeros(size(vpathx,1),n);
for i=1:n
    i
    [temp,~,~]=get_consistency_measures(windowsize,vpathxnull{i},hmm{ID}.K,false);
    nullconsistencies(:,i)=temp;
end
nullconsistencies=nullconsistencies(:);
nullconsistencies(isnan(nullconsistencies))=[];
nullconsistencies2=100*[prctile(nullconsistencies,5),prctile(nullconsistencies,95)]; 
sprintf('%.1f percent expected consistency from phase-randomised null',mean(nullconsistencies2))

%%
%{
Optional: Visualise HMM states with makehuman faces.
Need this for Figure 2, Panel A, pictures of faces
You can skip this step as I've provided the resulting faces in
folder RenderedFacesForFigure2.
Find each AU's contribution by summing across freq bands.
Contributions are rescaled to 0 to 1.
%}
n=5; %at the least, show the n largest AU contributions (rescaled to 0 to 1)
visualiseAUs=zeros(8,nAUs);

h1=hmm{ID};
means=getMeans(h1);
means2=reshape(means,[],n_ylabels,length(these_AUs));

for i=1:8
    temp=squeeze(means2(i,:,:));
    temp=sum(temp);
    temp_sorted=sort(temp,'descend');
    %{
    If more than n AUs have contribution > 0, show all positive contributions
    If less than n AUs have contributions > 0, show exactly the n largest
    contributions
    %}
    temp=rescale(temp,'InputMin',min(temp_sorted(n+1),0));  
    visualiseAUs(i,:)=temp;
end
%figure; imagesc(visualiseAUs); set(gca,'XTickLabel',aulist(these_AUs),'XTick',1:14); colorbar; xlabel('AU'); ylabels('HMM state');

for i=1:8 %Save data for makehuman in .json text format, but .facs file format
    aulist2=cellfun(@(x) ['AU',int2str(eval(x))],aulist,'UniformOutput',false);
    s=cell2struct(num2cell(visualiseAUs(i,:)'),aulist2(these_AUs));
    txt=jsonencode(s,'PrettyPrint',true);
    txt=strrep(txt,"AU","");
    DF=['C:\Users\Jayson\Documents\makehuman\v1py3\data\facs\','DISFA_HMM_n5_state',int2str(i),'.facs'];
    fid=fopen(DF,'w');
    raw=fwrite(fid,txt);
    fclose('all');
end
%Now go to makehuman, render these faces, and save, before returning to
%make Figure 2.

%%
%Figure 2

%Make panels
figure('Position',[100,100,1000,800]);
p=panel();
p.pack({.32,.68});
p(1).pack('h',{.98,.02});
p(2).pack('h',{.82,.18});
p(2,1).pack({.57,.08,.35});
pStateMeans=p(1,1); pV=p(2,1,1); pCommon=p(2,1,2); pConsistency=p(2,1,3); 
pTransition=p(2,2);
pStateMeans.pack({.65,.35});

%Panel A: Top row, avatar faces
pStateMeans(1).pack('h',8);
for i=1:8
    pStateMeans(1,i).select();
    temp=imread(['D:\FORSTORAGE\Data\Melancholia\MyCodeAndResults\makehuman_render\RenderedFacesForFigure2\n5_state',int2str(i),'.png']);
    %MODIFY THIS to where your makehuman renders
    temp=temp(70:end,210:end-150,:);
    imshow(temp);
end

%Panel A: Bottom row, State means
pStateMeans(2).pack('h',8);
h1=hmm{ID};
means=getMeans(h1);
means2=reshape(means,[],n_ylabels,length(these_AUs));
period=1./frqvalues;
for i=1:hmm{1}.K
    [col,row]=ind2sub([8,1],i);
    pStateMeans(row+1,col).select();  
    image=squeeze(means2(i,:,:));
    imagesc(1:length(these_AUs),flipud(log2(period)),image);
    set(gca, 'XTick', 1:length(these_AUs), 'XTickLabel', aulist(these_AUs));
    title(sprintf('State %i', i)); 
    set(gca,'clim',[-1,1]); %set(gca,'clim',makesymmetric(get(gca,'clim')));
    Yticks = (2.^(fix(log2(min(period))):fix(log2(max(period)))));
    Yticks=Yticks(1:2:end);
    YTickLabel=num2str(1./Yticks',1);
    set(gca,'YLim',log2([min(period),max(period)]), ...
        'YDir','reverse', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',{},...
        'layer','top',...
        'box','on',...
        'FontSize',6,...
        'FontSizeMode','manual')
    colormap(gca,brewermap(100,'*RdBu'))
end
colorbar('Position',[.96,.84,0.012,0.1]);
pStateMeans(2,1).select(); xlabel('Action unit'); ylabel('Frequency (Hz)'); set(gca,'YTickLabel',YTickLabel);

%Panel B: Viterbi paths
pV.select();
downsample=configurations{ID}.downsample;
colormap_Viterbi=brewermap(h1.K,'Accent');
colormap(gca,colormap_Viterbi);
imagesc(vpathx'); ylabel('Participant'); 
set(gca,'xtick',[])
set(gca,'xlim',[0.5,size(vpathx,1)+0.5]);
set(gca,'ylim',[0.5,size(vpathx,2)+0.5]);

%Panel C: Most common state
pCommon.select();
imagesc(vmostcommon); colormap(gca,brewermap(h1.K,'Accent')); 
set(gca,'xlim',[0.5,size(vpathx,1)+0.5]);
axis off;

%Panel D: Consistency
pConsistency.select(); 
tdown=(1:nframesdown)/downsample;
smwin=80;
plot(tdown,100*smoothdata(vconsistency,'Gaussian',smwin),'b'); 
ylow=100*smoothdata(prctile(vconsistencies,5,2),'Gaussian',smwin);
yhigh=100*smoothdata(prctile(vconsistencies,95,2),'Gaussian',smwin);
myshade(tdown,ylow',yhigh','b',0.3); %blue shading. Confidence intervals
myshade(tdown,nullconsistencies2(1),nullconsistencies2(2),'k',0.2); %grey shading. Null distribution
xlabel('Time (s)'); ylabel('Consistency (%)'); %ylim([0,100]);
xlim([0,max(tdown)]); 
draw_eventlinesInPlot(event_times,1);
axis tight;

%Panel E: Visualise transition matrix
pTransition.select();
layout='force'; %'layered' or 'force'
prc=80; %cutoff percentile. 80 is default. 0 just leave all connections in
adj=getTransProbs(h1);
pplot=plot_HMM_transitionmatrix(adj,layout,prc,colormap_Viterbi);
axis off;

%Some small edits
p.fontname='TimesNewRoman';
pV.marginbottom=10;
pCommon.marginbottom=5;
pConsistency.marginbottom=10;
pStateMeans.de.marginleft=0;
pStateMeans.de.marginright=2;
pStateMeans.de.marginbottom=0;
pStateMeans.de.margintop=0;
pTransition.margin=0;
pStateMeans.marginright=0;
p(1,2).margin=0;
p(1).marginbottom=10;

%Coloured rectangles 
for i=1:8
    pStateMeans(2,i).select();
    [xaf,yaf]=ds2nfu(0,log2(max(period))); %ds2nfu from https://au.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion
    annotation('rectangle',[xaf-0.006,yaf+0.03,0.018,0.05],'Color','none','FaceColor',colormap_Viterbi(i,:));
end

%Reduce font size 
for i=1:8
    pStateMeans(2,i).select();
    set(gca,'FontSize',8);
    set(gca,'TitleFontSizeMultiplier',10/8);
end

%Write video short descriptions 
pV.select();
video_descriptions={'Happy','Surprise','Disgust','Fear','Sad','Fear'};
temp=[0,event_times.*downsample,size(vpathx,1)]; %x coordinates demarcation between videos
temp=mean([temp(1:end-1);temp(2:end)]); %midpoint of times for each video
for i=1:length(temp)
    [xaf,yaf]=ds2nfu(temp(i),0);
    annotation('textbox',[xaf-0.03,yaf,0,0],'String',video_descriptions{i},'FitBoxToText','off','LineStyle','none');
end
clear video_descriptions;
[xaf,yaf]=ds2nfu(0,0);
annotation('textbox',[xaf-0.04,yaf,0,0],'String','Stimulus','FitBoxToText','off','LineStyle','none');

%Horizontal white lines between subjects
pV.select();
for i=1:size(vpathx,2)-1
    [xaf,~]=ds2nfu(0,i+0.5);
    [xaf2,yaf]=ds2nfu(size(vpathx,1),i+0.5);
    annotation('line',[xaf,xaf2],[yaf,yaf],'Color','w','LineWidth',2);
end

%Vertical black lines between stimuli
for i=1:length(event_times)
    xcoord=event_times(i)*downsample;
    pV.select();
    [xaf,yaf]=ds2nfu(xcoord,0.5);
    [xaf2,yaf2]=ds2nfu(xcoord,size(vpathx,2)+0.5);
    annotation('line',[xaf,xaf2],[yaf,yaf2],'Color','k','LineWidth',2);    
    pCommon.select();
    [xaf,yaf]=ds2nfu(xcoord,0.5);
    [xaf2,yaf2]=ds2nfu(xcoord,1.5);
    annotation('line',[xaf,xaf2],[yaf,yaf2],'Color','k','LineWidth',2);
end

%Label each panel
fontsize=14; 
annotation('textbox',[0.02,0.975,0,0],'String','a','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.02,0.68,0,0],'String','b','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.02,0.32,0,0],'String','c','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.02,0.265,0,0],'String','d','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.84,0.62,0,0],'String','e','FitBoxToText','on','LineStyle','none','FontSize',fontsize);

%%
%Optional: Supplementary videos: video examples of people in HMM states
nparticipant=6; %I did this for subjects 6, 27, and 7
downsample=configurations{ID}.downsample;
tdown=(1:nframesdown)/downsample;
figure; mypcolor(tdown,squeeze(Gammax(:,nparticipant,:))'); 

%{
Look at above plot to find time points when participants were in particular
HMM states.
Not all subjects exhibited all HMM states.
After getting representative time points in exemplar subjects, trim
videos with ffmpeg
Trim times given below -->

ExampleSubject01 is DISFA Subject 006
state 1: '00:23','00:27'
state 2: '00:36','00:40'
state 3: '1:32','1:36'
state 4: '1:22','1:26'
state 5: '1:56','2:00'
state 6: '2:28','2:32'

ExampleSubject02 is DISFA Subject 032
state 1: '00:30','00:34'
state 3: '1:11','1:15'
state 4: '1:23','1:27'
state 5: '2:45','2:49'
state 6: '3:09','3:13'
state 7: '3:34','3:38'
state 8: '2:33','2:37'

ExampleSubject03 is DISFA Subject 007
State 1: '00:24','00:28'
State 2: '00:35','00:39'
State 3: '1:20','1:24'
State 4: '2:01','2:05'
State 6: '3:32','3:36'
State 7: '3:04','3:08'
%}
