%{
Run HMM-MAR from https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide,
on melancholia data. Run melancholia_analysis.m first
This code generally uses CN=1, MEL=2, NONMEL=3 for indexing
%}

%In HMM-MAR/train/initGamma_random.m, set to rng(1) for reproducible results or rng('shuffle') for different randomization seed each time.
format compact
hmmtype='timefreq';%options are 'normal', 'timefreq'
savefolder='D:\\FORSTORAGE\\Melancholia_figures';

if strcmp(hmmtype,'normal')
    X=get_X(a);
    X=X(a_valid);
    for i=1:length(X) %replace 0 variance time series w a single 0.01
        array=X{i};
        for j=1:width(array)
            if var(array(:,j))==0
                X{i}(1,j)=0.1;
            end
        end  
    end
    N=length(X); %subjects
    ttrial=height(X{1}); %time points
    nregions=length(these_AUs); %number of AUs 
    %vY already defined in melancholia_analysis
    X=cell2mat(X);
elseif strcmp(hmmtype,'timefreq')
    X=cfss8;
    assert(size(X,2)>20);
    nregions=size(X,2);
end

T=repmat({[ttrial]},N,1); % length of data for each session
Tsubject=T;
T=cell2mat(T);
Tsubject=T;
nsubs=length(T);

%%
%Get default configuration parameters. Edit get_defaults_HMM.m if you want
%non-standardized HMM
[windowsize_sec,smwin,configurations]=get_defaults_HMM(8);
for i=1:length(configurations)
    configurations{i}.Fs=25;
end

%% Run HMM and perform statistical testing on the HMM results vs the conditions
[hmm,Gamma,Gammaup,Xi,vpath,fehist,maxFO,newframerate] = runHMM_JJ1(configurations,X,T);
windowsize=round(windowsize_sec*newframerate);
ID=1;

%%
%Reording states to put happy ones together. 
%ord=[6,8, 2,4,1,7, 3,5]; %change old index 8 to second (new index 2), etc

ord=[7,1,  6,2,5,   3,8,4]; %change old index 5 to third (new index 3) etc. Works for covtype='full' or 'diag'. Also for 'nostandard'.

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
Xi{ID}=Xi{ID}(:,ord,ord);

%%
test_group{ID}.p_fractional_occupancy=test_group{ID}.p_fractional_occupancy(ord);
%%..then rerun fractional occupancy
%%statistics. BH or Storey??

%%
rng(1); 
tt={}; 
test_group = cell(length(configurations),1);
options_test = struct();
options_test.subjectlevel = 0;
options_test.Nperm = 10000;
for i = 1:length(configurations)
    tt{i} = hmmtest(Gammaup{i},T,Tsubject,Y2,options_test,hmm{i});  
    test_group{i} = tt{i}.grouplevel; % only doing group-level testing
end

pvals=test_group{1}.p_fractional_occupancy;
fdr=mafdr(pvals);
[pvals,fdr]

%%
h1=hmm{ID};
[ii,savesuffix] = get_savesuffix(hmm,ID,maxFO,fehist,savefolder); %% Prepare save file name
%HMM post-processing including new variables for downsampling, reshape
%vpath/Gamma, and get individual subject masks
[nframesdown,ntotalframesdown,Tdown,vpathx,Gammax,Gammac,Masks_subject]=hmm_postprocessJJ(hmm,vpath,ID,nsubs,Gamma);
%[vconsistency,vmostcommon,vconsistencies] = get_consistency_measures(windowsize,vpathx,hmm{ID}.K,true);
means=getMeans(h1);
means2=reshape(means,[],n_ylabels,length(these_AUs));
playtone();

%%
%Get most common state at each time point, and consistency. Takes a while.
[vconsistencyCN,vmostcommonCN,vconsistenciesCN] = get_consistency_measures(windowsize,vpathx(:,v2CN),h1.K,true);
[vconsistencyMEL,vmostcommonMEL,vconsistenciesMEL] = get_consistency_measures(windowsize,vpathx(:,v2MEL),h1.K,true);

%Lists of which time points in super-array belong to CN vs MEL subjects
Masks_groups=zeros(ntotalframesdown,1);
index=1;
for nsub=1:nsubs
    Masks_groups(index:index+Tdown(nsub)-1)=find(vY(nsub,:));
    index=index+Tdown(nsub);
end
Masks_groups={find(Masks_groups==1)',find(Masks_groups==2)',find(Masks_groups==1|Masks_groups==2)'}; 
playtone();

%Get transition probabilities for each group
[P,Pi] = getMaskedTransProbMats (X,Tdown,h1,Masks_groups,Gamma{ID},Xi{ID}); %using masks for each group
Pn=cellfun(@(x) myTransProbs(x),P,'UniformOutput',false); %remove self-connnection
Ppersist=cell2mat(cellfun(@(x) diag(x),P,'UniformOutput',false)); %persistence probabilities


%%
%{
Visualisation of Melancholia->HMM states with makehuman faces
Find each AU's contribution by summing across freq bands
Contributions shown are rescaled to 0 to 1
%}
n=5; 
visualiseAUs=zeros(8,nAUs);
for i=1:8
    temp=sq(means2(i,:,:));
    temp=sum(temp);
    temp_sorted=sort(temp,'descend');
    %{
    If more than n AUs have contribution > 0, show all positive contributions
    If less than n AUs have contributions > 0, show exactly the n largest
    contributions
    e.g. DISFA_HMM_state1
    %}
    %
    temp=rescale(temp,'InputMin',min(temp_sorted(n+1),0));  
    visualiseAUs(i,:)=temp;
end
figure; imagesc(visualiseAUs); set(gca,'XTickLabel',aulist(these_AUs),'XTick',1:14); colorbar;

%%
%Save data for FACSHuman in .json text format, but .facs file format
for i=1:8
    aulist2=cellfun(@(x) ['AU',int2str(eval(x))],aulist,'UniformOutput',false);
    s=cell2struct(num2cell(visualiseAUs(i,:)'),aulist2(these_AUs));
    txt=jsonencode(s,'PrettyPrint',true);
    txt=strrep(txt,"AU","");
    DF=['C:\Users\Jayson\Documents\makehuman\v1py3\data\facs\','Melancholia_HMMdiag_n5_state',int2str(i),'.facs'];
    fid=fopen(DF,'w');
    raw=fwrite(fid,txt);
    fclose('all');
end
%Now go to makehuman and render these faces

%%
%Figure 6

%Prepare the layout
colormap_Viterbi=brewermap(h1.K,'Accent');
figure('Position',[50,50,1000,700])
p=panel();

%p.pack({.8,.2}); pV=p(1); pConsistency=p(2); pV.pack(2); %When HMM is
%nonstandradised

p.pack({.3,.7});
p(1).pack('h',{.99,[]})
p(1,1).pack({.65,.45});
pImages=p(1,1,1); pStateMeans=p(1,1,2); %Panel A 

p(2).pack('h',{.5,.5});
pV=p(2,1); pV.pack(2); %Panel B
p(2,2).pack({.15,.30,[]});
p(2,2,3).pack('h',{.55,.37,[]});
pCommon=p(2,2,1); %Panel C
pCommon.pack(2);
pConsistency=p(2,2,2); %Panel D
pGraph=p(2,2,3,1); %Panel E
pMatrix=p(2,2,3,2); %Panel F
downsample=configurations{ID}.downsample;
pImages.pack('h',8);
pStateMeans.pack('h',8);

%Panel A
for i=1:8
    pImages(i).select();
    temp=imread(['C:\Users\Jayson\Documents\makehuman\v1py3\render\melancholia_HMM_processed\n5_state',int2str(i),'.png']);
    temp=temp(60:end-50,210:end-160,:);
    imshow(temp);
end

period=1./frqvalues;
for i=1:hmm{1}.K
    %[col,row]=ind2sub([8,1],i);
    pStateMeans(i).select();  
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
        'FontSizeMode','manual',...
        'xlim',[0.5,size(image,2)+0.5])
    colormap(gca,brewermap(100,'*RdBu'))
end
colorbar('Position',[.945,.79,0.012,0.1]);
pStateMeans(1).select();xlabel('Action unit'); ylabel('Frequency (Hz)'); set(gca,'YTickLabel',YTickLabel); 

%Panel B
v2={v2CN,v2MEL};
for j=1:2
    pV(j).select();
    colormap(gca,colormap_Viterbi);
    imagesc(vpathx(:,v2{j})'); 
    if j==1
        ylabel('Participant'); 
    end
    set(gca,'xtick',[])
    set(gca,'xlim',[0.5,size(vpathx,1)+0.5]);
    set(gca,'ylim',[0.5,sum(v2{j})+0.5]);
end
clear v2;

%Panel C
pCommon(1).select();
imagesc(vmostcommonCN); colormap(gca,colormap_Viterbi); 
set(gca,'xlim',[0.5,size(vpathx,1)+0.5]);
axis off;
pCommon(2).select();
imagesc(vmostcommonMEL); colormap(gca,colormap_Viterbi); 
set(gca,'xlim',[0.5,size(vpathx,1)+0.5]);
axis off;

tdown=(1:nframesdown)/downsample;

%Panel D
smwin=80;
pConsistency.select(); 
plot(tdown,100*smoothdata(vconsistencyCN,'Gaussian',smwin),'b'); 
ylow=100*smoothdata(prctile(vconsistenciesCN,5,2),'Gaussian',smwin);
yhigh=100*smoothdata(prctile(vconsistenciesCN,95,2),'Gaussian',smwin);
myshade(tdown,ylow',yhigh','b',0.3)
%myshade(1:length(ylow),nullconsistencies2(1),nullconsistencies2(2),'k',0.2);

hold on
plot(tdown,100*smoothdata(vconsistencyMEL,'Gaussian',smwin),'k'); 
ylow=100*smoothdata(prctile(vconsistenciesMEL,5,2),'Gaussian',smwin);
yhigh=100*smoothdata(prctile(vconsistenciesMEL,95,2),'Gaussian',smwin);
myshade(tdown,ylow',yhigh','k',0.3)
%myshade(1:length(ylow),nullconsistencies2(1),nullconsistencies2(2),'k',0.2);
hold off

xlabel('Time (s)'); ylabel('Consistency (%)'); %ylim([0,100]);
xlim([0,max(tdown)]); temp=get(gca,'ylim'); ylim([temp(1),99]);
draw_eventlinesInPlot(event_times,1);

%Panel E. Transition array: Digraph representation, for all subjects
pGraph.select();
layout='force'; %'layered' or 'force'
prc=80; %cutoff percentile. 80 is default. 0 just leave all connections in
pplot=plot_HMM_transitionmatrix(Pn{3},layout,prc,colormap_Viterbi);
axis off

%Panel F: Transition array: MEL-CN
pMatrix.select();
imagesc(Pn{2}-Pn{1}); 
xlabel('To'); ylabel('From');
set(gca,'ylim',[0.5,8.5]);
set(gca,'xlim',[0.5,8.5]);
set(gca,'XTick',1:8);
set(gca,'YTick',1:8);
caxis(makesymmetric(get(gca,'clim')));
colormap(gca,brewermap(100,'*BrBG'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5 - continued
p.fontname='TimesNewRoman';
pCommon.margin=5;
pStateMeans.marginleft=5;
p(1,1).de.marginleft=0;
p(1,1).de.marginright=2;
p(1,1).de.marginbottom=0;
p(1,1).de.margintop=5;
pStateMeans.margintop=0;

%Coloured rectangles 
for i=1:8
    pStateMeans(i).select();
    [xaf,yaf]=ds2nfu(0,log2(max(period)));
    annotation('rectangle',[xaf-0.006,yaf+0.023,0.016,0.05],'Color','none','FaceColor',colormap_Viterbi(i,:));
end
%Reduce font size 
for i=1:8
    pStateMeans(i).select();
    set(gca,'FontSize',8);
    set(gca,'TitleFontSizeMultiplier',10/8);
end

%Write video short descriptions 
video_descriptions={'Comedy','Sad','Funny'};
temp=[0,event_times.*downsample,size(vpathx,1)]; %x coordinates demarcation between videos
temp=mean([temp(1:end-1);temp(2:end)]); %midpoint of times for each video
pV(1).select();
[xaf,yaf]=ds2nfu(0,0);
annotation('textbox',[xaf-0.03,yaf,0,0],'String','Stimulus','FitBoxToText','off','LineStyle','none');
for j=1
    pV(j).select();
    for i=1:length(temp)
        [xaf,yaf]=ds2nfu(temp(i),0);
        annotation('textbox',[xaf-0.03,yaf,0,0],'String',video_descriptions{i},'FitBoxToText','off','LineStyle','none');
    end
end
clear video_descriptions;

%MEL and CN labels
pV(1).select()
fontsize=10;
[xaf,yaf]=ds2nfu(size(vpathx,1)/2,sum(v2CN));
annotation('textbox',[xaf,yaf+0.035,0,0],'String','Control','FitBoxToText','off','LineStyle','none','FontSize',fontsize,'HorizontalAlignment','center');
pV(2).select()
[xaf,yaf]=ds2nfu(size(vpathx,1)/2,sum(v2MEL));
annotation('textbox',[xaf,yaf+0.035,0,0],'String','Melancholia','FitBoxToText','off','LineStyle','none','FontSize',fontsize,'HorizontalAlignment','center');
pCommon(1).select();
[xaf,yaf]=ds2nfu(size(vpathx,1)/2,1);
annotation('textbox',[xaf,yaf+0.04,0,0],'String','Control','FitBoxToText','off','LineStyle','none','FontSize',fontsize,'HorizontalAlignment','center');
pCommon(2).select();
[xaf,yaf]=ds2nfu(size(vpathx,1)/2,1);
annotation('textbox',[xaf,yaf+0.04,0,0],'String','Melancholia','FitBoxToText','off','LineStyle','none','FontSize',fontsize,'HorizontalAlignment','center');

%Colour label for Transition matrix 
pMatrix.select();
pos=[.95,.155];
csize=[0.012,0.1];
fontsize=11; 
colorbar('Position',[pos(1),pos(2),csize(1),csize(2)],'TickLabels',{});
annotation('textbox',[pos(1)+0.0,pos(2)+0.1,0,0],'String','Melancholia > Control','LineStyle','none',...
    'FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','Bottom');
annotation('textbox',[pos(1)+0.0,pos(2)+csize(2)-0.1,0,0],'String','Melancholia < Control','LineStyle','none',...
    'FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','Top');
clear csize pos

v2={v2CN,v2MEL};
%Horizontal white lines between subjects
for j=1:2 
    pV(j).select();
    for i=0:sum(v2{j})
        [xaf,~]=ds2nfu(0,i+0.5);
        [xaf2,yaf]=ds2nfu(size(vpathx,1),i+0.5);
        annotation('line',[xaf,xaf2],[yaf,yaf],'Color','w','LineWidth',1);
    end
end
%Vertical black lines between stimuli
for i=1:length(event_times)
    xcoord=event_times(i)*downsample;
    for j=1:2
        pV(j).select();
        [xaf,yaf]=ds2nfu(xcoord,0.5);
        [xaf2,yaf2]=ds2nfu(xcoord,sum(v2{j})+0.5);
        annotation('line',[xaf,xaf2],[yaf,yaf2],'Color','k','LineWidth',2);    
        pCommon(j).select();
        [xaf,yaf]=ds2nfu(xcoord,0.5);
        [xaf2,yaf2]=ds2nfu(xcoord,1.5);
        annotation('line',[xaf,xaf2],[yaf,yaf2],'Color','k','LineWidth',2);
    end
end
clear v2

%Label each panel
fontsize=13; 
annotation('textbox',[0.03,0.975,0,0],'String','a','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.03,0.68,0,0],'String','b','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.52,0.69,0,0],'String','c','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.52,0.59,0,0],'String','d','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.52,0.35,0,0],'String','e','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[0.73,0.35,0,0],'String','f','FitBoxToText','on','LineStyle','none','FontSize',fontsize);


%%
%Histogram and print p-values for fractional-occupancies across groups

[FO,ntrials] = getFractionalOccupancy (Gammaup{ID},T,configurations{ID});

pvals=test_group{1}.p_fractional_occupancy;
fdr=mafdr(pvals);
[pvals,fdr]

groupcolours={'b','k'}; %CN, MEL, NONMEL
figure;
for i=1:h1.K
    nedge=ceil(sqrt(h1.K));
    subplot(nedge,nedge,i); hold on;
    for j=1:width(Y2)
        histogram(FO(Y2(:,j),i),'FaceColor',groupcolours{j});
    end
    hold off;
    
    ylabel('Fractional occupancy');
    pval=test_group{ID}.p_fractional_occupancy(i);
    title(sprintf('State %i, p=%.2f',i,pval));  
end
sgtitle('Fractional Occupancy of each state');
legend('CN','MEL');

%saveas(gcf,fullfile(savefolder,sprintf('%id_%s',ii,savesuffix)));

%figure; plot(maxFO);title('maxFO'); %too high indicates individual states are mapping onto individual subjects

%%
%Plot transition probabilities for each group separately
%{
figure; 
subplot(2,2,1); 
bargraph=bar(Ppersist(:,[1,2]));
bargraph(1).FaceColor='b';
bargraph(2).FaceColor='k';
minimum=min(min(Ppersist)); maximum=max(max(Ppersist)); range=maximum-minimum;
padding=0.2;
ylim([minimum-padding*range,1]); legend({'CN','MEL'}); title('Persistence probabilities'); xlabel('HMM state');

subplot(2,2,3); imagesc(Pn{1}); title('Transition probabilities: CN'); colorbar;
subplot(2,2,4); imagesc(Pn{2}); title('Transition probabilities: MEL'); colorbar;
subplot(2,2,2); imagesc(Pn{2}-Pn{1}); title('Transition probabilities: MEL-CN'); colorbar;
%Each row sums to 1, so it is a right stochastic matrix
%columns are 'to' states. rows are 'from' states. Q(t+1) = Q(t) * MATRIX

%saveas(gcf,fullfile(savefolder,sprintf('%ie_%s',ii,savesuffix)));

%}
%%
%Which transition probabilities and their differences across groups can account for difference in fractionlal occupancies? - compared across groups (w pvals). Quite
%confusing code. 

[P2,P2i] = getMaskedTransProbMats (X,Tdown,h1,Masks_subject,Gamma{ID},Xi{ID}); %Using masks for each participant separately
P2n=cellfun(@(x) myTransProbs(x),P2,'UniformOutput',false); %remove self-connnection

meandiffs=zeros(h1.K,h1.K); pvals=meandiffs; prights=meandiffs; plefts=meandiffs;
for i=1:h1.K
    for j=1:h1.K
        temp=cellfun(@(x) x(i,j),P2);
        x=temp(v2CN); y=temp(v2MEL);
        [h,p]=ttest2(x,y,'Tail','both');
        [h,pleft]=ttest2(x,y,'Tail','left');
        [h,pright]=ttest2(x,y,'Tail','right');
        
        meandiffs(i,j)=mean(y)-mean(x);
        pvals(i,j)=p; prights(i,j)=pright; plefts(i,j)=pleft;
    end
end

figure;
subplot(2,2,1); imagesc(meandiffs>0); title('y(MEL) > x(CN)??');
subplot(2,2,2); imagesc(pvals<0.05); title('pval<0.05?');
subplot(2,2,3); imagesc(plefts<0.05); title('plefts<0.05? for CN<MEL');
subplot(2,2,4); imagesc(prights<0.05); title('prights<0.05? for CN>MEL');

%old 1,7,4 is new 2,1,8
states=1:h1.K; %list of states
T1a=plefts(2,setdiff(states,[]))'; %Transitions from state 2 and 1: expect CN < MEL (left)
T1b=plefts(1,setdiff(states,[]))';
T2a=prights(setdiff(states,[]),2); %Transitions to 2 and 1: expect CN > MEL (right)
T2b=prights(setdiff(states,[]),1);
T3=prights(8,setdiff(states,[]))'; %Transition from 8: expect CN > MEL (right)
T4=plefts(setdiff(states,[]),8); %Transition to 8: expect CN < MEL (left)
T5=[T1a,T1b,T2a,T2b,T3,T4];
clear T1a T1b T2a T2b T3 T4;
fdr=reshape(mafdr(T5(:),'BHFDR',false),size(T5,1),size(T5,2));
figure; imagesc(T5<0.05); %columns are [From 2, From 1, To 2, To 1, From 8, To 8]

Tpersist=[prights(1,1),prights(7,7),plefts(4,4)]; %Persistence probabilities. CN>MEL for 1 and 7, CN<MEL for 4

%fdr=reshape(mafdr(pvals(:)),size(pvals,1),size(pvals,2));
%figure; imagesc(fdr<0.05); title('pfdr<0.05?');
%%
%Visualise transition matrix separately for each group (Optional)

layout='force'; %'layered' or 'force'
prc=20; %cutoff percentile. 80 is default. 0 just leave all connections in

figure;
groupnames={'CN','MEL'}; %has to be in this order because of Masks
for i=1:2
    subplot(2,1,i);
    adj=Pn{i};
    pplot=plot_HMM_transitionmatrix(adj,layout,prc,colormap_Viterbi);
    title(groupnames{i});
end
sgtitle('Transition matrix');
fig=gcf;
pos=fig.Position;
fig.Position(4)=pos(4)*2;
fig.Position(2)=pos(2)-500;
%saveas(gcf,fullfile(savefolder,sprintf('%if_%s',ii,savesuffix)));

%%
function [M]=get_X(m)
%Get the variable X in preparation for HMM-MAR
M=squeeze(num2cell(m.cube,[1,2])); %CHANGE CUBE TO CUBEN TO GET NORMALIZED VERSION
end

function y=myTransProbs(x)
x(eye(height(x))==1)=0;
multiplier=1./sum(x,2);
y=x.*multiplier;
end

function y=mycell2mat(x)
%{
x is 1-dimensional cell array of 2-dimensional arrays
Combine into a 3D array
%}
y=zeros(size(x{1},1),size(x{1},2),length(x));
for i=1:size(y,3)
    y(:,:,i)=x{i};
end
end

