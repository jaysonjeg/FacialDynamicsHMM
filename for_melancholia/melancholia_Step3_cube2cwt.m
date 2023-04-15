%{
This script models melancholia AU dynamics. Continue after melancholia_analysis_mat2cube.m
My PC has 64GB RAM which was required to store the large CWT arrays
Also needs variables from melancholia_Step2_mat2cube: vY, vCN, vMEL
%}

%%
%Preparing for Figure 5
%Get CWT of individual participants. First run this code section with maxframe==3160 (corresponds to end of the first stimulu) (needed for Figure 5 panel B),
%Takes 12sec

maxframe=3160; %set to 8804 for all frames,or to 3160 for just stand-up comedy

cube=a.cube(1:maxframe,:,a_valid); %include valid subjects. ntimepoints x nAUs x nsubs
[cfssraw,frq,coi] = get_cwt(cube,Fs,true);
[period,coi]=get_periods(frq,coi);
cfss=get_cwt_mag(cfssraw,cube);  %Replace complex coeffs with just magnitudes. cfss is nfreqs x ntimepoints x nsubs x nAUs
clear cfssraw
cfssm=[]; %mean of magnitudes of CWT, for each group. nfreqs * ntimepoints * ngroups * nAUs
for nDiagnosis=1:3   
    inds=find(vY(:,nDiagnosis));
    cfssm=cat(3,cfssm,mean(cfss(:,:,inds,:),3)); %nfreqs x ntimepoints x ngroups x nAUs
end

%We are labelling those output variables that correspond to
%maxframe==3160 (or to stand-up comedy alone) with 'z' at the end
nAU=7; %AU12 is number 7. Focusing on this AU for Figure 5 panel B
valCNz=cfssm(:,:,invdxmap('CN'),nAU); 
valMELz=cfssm(:,:,invdxmap('MEL'),nAU);
periodz=period; coiz=coi;

playtone(1000,0.08,0.3);
clear inds cube

%%
%Preparing for Figure 5
%Get CWT of individual participants. Similar as above code section but with maxframe==8804 (corresponds to all stimuli together (equals nframes) (for Figure 5 panel A), then move to next step
%Takes 60 sec

maxframe=8804; %set to 8804 for all frames

cube=a.cube(1:maxframe,:,a_valid); %include valid subjects. ntimepoints x nAUs x nsubs
[cfssraw,frq,coi] = get_cwt(cube,Fs,true);
[period,coi]=get_periods(frq,coi);
cfss=get_cwt_mag(cfssraw,cube);  %Replace complex coeffs with just magnitudes. cfss is nfreqs x ntimepoints x nsubs x nAUs
clear cfssraw
cfssm=[]; %mean of magnitudes of CWT, for each group. nfreqs * ntimepoints * ngroups * nAUs
for nDiagnosis=1:3   
    inds=find(vY(:,nDiagnosis));
    cfssm=cat(3,cfssm,mean(cfss(:,:,inds,:),3)); %nfreqs x ntimepoints x ngroups x nAUs
end
playtone(1000,0.08,0.3);
clear inds cube

%%
%{
%Optional: Get CWT for mean ts of each group
means_CN=mean(a.cube(:,:,CN),3,'omitnan');
means_NONMEL=mean(a.cube(:,:,NONMEL),3,'omitnan');
means_MEL=mean(a.cube(:,:,MEL),3,'omitnan');
allmeans=cat(3,means_CN,means_MEL,means_NONMEL);
[mcfssraw,frq,~]=get_cwt(allmeans(1:maxframe,:,:),Fs,true);
[~,~,coi3160]=get_cwt(allmeans(1:3160,:,:),Fs,true);
mcfss=get_cwt_mag(mcfssraw,allmeans);
%}

%%
%Figure 5—figure supplement 1. Plot mean CWT maps for each diagnosis group
nstart=1; %put to 40, to cut off higher frequencies
t = (0:nframes-1)/Fs;
for nDiagnosis=0:2
    figure;  sgtitle(dxmap(nDiagnosis)); %colormap('jet'); 
    for nAU=1:length(these_AUs)
        subplot(4,4,nAU);
        values=cfssm(:,:,nDiagnosis+1,nAU);
        plot_cwt(inferno(),values,period,Fs,coi);
        title(sprintf('AU%s: %s',aulist{these_AUs(nAU)},aunames(these_AUs(nAU))));
        xlabel('Time (s)'); ylabel('Frequency (Hz)'); colorbar;
        xline(125.6); xline(125.6+168.84);
    end
end


%%
%Figure 5 preparation. Get arrays of CWT map differences for 2 groups
%Get null distribution of effect size (difference in mean between groups).
%Takes quite a long time (depending on Nsurrogates)

Nsurrogates=1000; %number of surrogates. Default 1000 (try with Nsurrogates=10 to test it out first)
plot_these_AUs=1:14; %show these AUs only, default 1:14
nstart=1; %set to 42, to cut off higher frequencies than 1.02 Hz
groups_of_interest={'MEL','CN'}; %first group minus second group

numbers=[sum(vCN),sum(vMEL)]; %group 1 CN, group 2 MEL. vCN has 1s for subjects that are controls, 0s for patients, but only indexes subjects who have a_valid==1

t = (0:nframes-1)/Fs;
ind=cellfun(@(x) invdxmap(x), groups_of_interest);

cfssm_d=zeros(size(cfssm,1),size(cfssm,2),size(cfssm,4)); %differences
%mcfss_d=cfssm_d; d_d=cfssm_d;
cfssm_d_perc=cfssm_d; cfssm_d_trinary=cfssm_d; cfssm_d_CI=cfssm_d;

for nAU=plot_these_AUs 
    string(nAU)+"/"+string(length(plot_these_AUs))
    tic
    cfssm_d(:,:,nAU)=cfssm(:,:,ind(1),nAU)-cfssm(:,:,ind(2),nAU);
    %mcfss_d=mcfss(:,:,ind(1),nAU)-mcfss(:,:,ind(2),nAU); %mcfss for 'CWTs of group means'
    %d_d=diff_of_cfssm - diff_of_mcfss;    
    %surrogates_AU=surrogates(:,:,:,nAU); %123 x 8804 x 1000
    surrogates_AU=NaN(size(cfss,1),maxframe,Nsurrogates); %each surrogate is a difference between surrogate groups' CWTs' means
    for n=1:Nsurrogates
        idx=randperm(size(cfss,3)); %random permutation of subjects in vCN and vMEL
        inds1=(idx(1:numbers(1))); %group 1
        inds2=(idx((numbers(1)+1):(numbers(1)+numbers(2)))); %group 2
        mean1=mean(cfss(:,:,inds1,nAU),3);
        mean2=mean(cfss(:,:,inds2,nAU),3);
        surrogates_AU(:,:,n)=mean2-mean1; %group 2 minus group 1
    end    
    cfssm_d_perc(:,:,nAU) = percentile_array(cfssm_d(:,:,nAU),surrogates_AU); %Get percentile relative to nulls
    %cfssm_d_trinary(:,:,nAU)=percentile_to_trinary(cfssm_d_perc(:,:,nAU),2.5);
    cfssm_d_CI(:,:,nAU)=exclude_inside_CI(cfssm_d(:,:,nAU),cfssm_d_perc(:,:,nAU),2.5);    
    sig=cfssm_d_CI;
    sig(sig~=0)=1;
    toc
end
playtone(1000,0.08,0.3);
%%
%Figure 5. Plot CWT group differences
cmap_a=inferno();
cmap_b=viridis();

figure('Position',[50,50,800,800]);
p=panel();
p.pack({.7,.3});

%Panel A
p(1).pack(4,4);
%plot_these_AUs=[1,4,5,9,12,13,14];
plot_these_AUs=1:14;
for nAU=plot_these_AUs 
    nAU
    [col,row]=ind2sub([4,4],nAU);
    if row==4
        col=col+1; %to centre the last row
    end   
    [ax1,ax2]=plot_cwt(cmap_a,cfssm_d(nstart:end,:,nAU),period,Fs,coi,sig(:,:,nAU));
    p(1,row,col).select([ax1,ax2]);
    [xaf,yaf]=ds2nfu(max(get(ax1,'XLim')),min(get(ax1,'YLim')));    
    set(ax1,'clim',[-0.18,0.06]); %manually set colour limits. Default [-0.18,0.05]
    set(ax2,'clim',[-0.18,0.06]); 
    title(sprintf('AU%s: %s',aulist{these_AUs(nAU)},aunames(these_AUs(nAU))));
    xline(ax1,125.6); xline(ax1,125.6+168.84);
    xline(ax2,125.6); xline(ax2,125.6+168.84);
end
p(1,3,1).select(); xlabel('Time (s)');ylabel('Frequency (Hz)'); 

p(2).pack('h',{.32,.32,.04,.32,});

%Panel B (left)
p(2,1).select();
plot_cwt(cmap_b,valCNz,periodz,Fs,coiz); set(gca,'clim',[0,0.25]); title('Controls');
xlabel('Time (s)'); ylabel('Frequency (Hz)'); 

%Panel B (right)
p(2,2).select();
axMEL=plot_cwt(cmap_b,valMELz,periodz,Fs,coiz); set(gca,'clim',[0,0.25]); title('Melancholia');

%Panel C
p(2,4).select();
plot_cwt(cmap_a,valMELz-valCNz,periodz,Fs,coiz); set(gca,'clim',[-0.18,0.06]); title('Difference');

%Put in lines at the laugh times
for j=[1,2,4]
    p(2,j).select();
    for i=1:length(laugh_times)
        xline(laugh_times(i));
    end
end

%Manually resize figure if needed then do this step
p.fontname='TimesNewRoman';
p.marginbottom=15;
p.margintop=10;
p(1).de.margin=10;
p(1).de.marginbottom=12;
p(1).marginbottom=12;
p(2).de.marginleft=10;
%p(1).de.marginleft=10;
%p.margin=15;
%p.de.margin=10;
%p.de.marginbottom=15;
clear nAU groups_of_interest

%Insert colour bars
csize=[0.01,0.1];
colorbar(axMEL,'Position',[.63,.15,csize(1),csize(2)]);

pos=[0.92,0.38];
colorbar(ax1,'Position',[pos(1),pos(2),csize(1),csize(2)],'TickLabels',{});
colorbar(ax2,'Position',[pos(1)+csize(1),pos(2),csize(1),csize(2)]);
fontsize=10; 
annotation('textbox',[pos(1)+0,pos(2)-0.01,0,0],'String','Melancholia < Control','FitBoxToText','on','LineStyle','none',...
    'FontSize',fontsize,'HorizontalAlignment','right','VerticalAlignment','Bottom');
annotation('textbox',[pos(1)+0,pos(2)+csize(2)+0.01,0,0],'String','Melancholia > Control','FitBoxToText','on','LineStyle','none',...
    'FontSize',fontsize,'HorizontalAlignment','right','VerticalAlignment','Top');

fontsize=14;
annotation('textbox',[.01,.99,0,0],'String','a','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[.01,.33,0,0],'String','b','FitBoxToText','on','LineStyle','none','FontSize',fontsize);

%%
%{
%Visualise pts/AUs with zero variance
allvars=squeeze(var(a.cube,1))';
imagesc(allvars==0)
%}
%%
%Prepare time-freq data for entry into HMM. 

usefreqbins=false; %default false. Need True for cwt2ML.m
cutoff_upper=5; %default 5hz. Set to a bit less than half of the Fs. So if Fs=10, then set this to 4.
cutoff_lower=0; %default 0hz
[cfss7,frqvalues,n_ylabels,cfss6] = process_cfss(cfss,frq,usefreqbins,cutoff_upper,cutoff_lower);
clear usefreqbins cutoff_upper cutoff_lower;
%Now go to melancholia_cwt2ML.m 

cfss7=cfss7(:,vCN|vMEL,:);
Y2=vY(vCN|vMEL,1:2);
v2CN=Y2(:,1); v2MEL=Y2(:,2);
clear cfss6;
ttrial=size(cfss7,1);
N=size(cfss7,2);
cfss8=reshape(cfss7,[],size(cfss7,3)); clear cfss7; %ntimepoints(nframes(fine)*nsubs(coarse)) * datapoints
playtone();

%{
clear cfss
save('all_variables_end_of_step_3.mat','-v7.3')
%}
%Now go to melancholia_cwt2hmm with hmmtype='timefreq' 



%%
function [array_trinary] = percentile_to_trinary(array_percentiles,cutoff)
%Given an array, convert values to -1, 0, or +1 depending on whether
%values are smaller than cutoff, or larger than 100-cutoff
array_trinary=zeros(size(array_percentiles));
array_trinary(array_percentiles<cutoff)=-1;
array_trinary(array_percentiles>(100-cutoff))=+1;
end 
 
function [out]=exclude_inside_CI(array,array_percentiles,cutoff)
%Given array, convert values inside the cutoff percentiles to zero
out=array; 
for i=1:size(array,1)
    for j=1:size(array,2)
        value=array_percentiles(i,j);
        if value > cutoff & value < (100-cutoff)
            out(i,j)=0;
        end
    end
end
end


