%{
Notes for Phil to read
The first few code sections labelled 'Preparation of data' (approx lines 0-310 has some info
specific to this study which is not as relevant for you. These sections
prepare many variables. You could just skip this section and start at line
310 with load('all_variables.mat')

Of the variables, the important ones are:
Variable a.cube here is the same as variable a_aus(:,these_AUs,:) in the
DISFA code. So probably best to rename your a_aus(:,these_AUs,:) to a.cube,
then to try some of these code sections.
Variable a.cubes contains same info as a.cube, but split by each stimulus
viewed. Needed for Figure 4
Variable a_valid is a list of 1 or 0, saying whether each participant has
valid data or not. Valid subjects (1s in a_valid) have all videos, don't have any action units with all zeros (causes HMM error), and are either CN or MEL

Info about the diagnostic group of each subject is stored in variables 'Y',
'CN', 'MEL' and 'NONMEL'
Varibles 'vY','vCN','vMEL','vNONMEL' are similar except only for subjects
in a_valid
Other important variables needed for most of the code are aulist,
these_AUs, these_diagnoses, these_stims, dxmap, invdxmap, Fs, event_times,
laugh_times

Unnecessary code sections are commented out.
%}

%%
%Preparation of data

Fs = 25; %sampling freq
these_stims= {'BillCosby','TheChamp','weather'}; %which stimuli to use: All available are {'BillCosby','TheChamp','weather','CryFreedom'};
these_AUs= [1,3,5:16];  %1:16; %don't include 17th AU which is AU45 blink, and AU02/05 which have 0 value in some subs
these_diagnoses={'CN','MEL','NONMEL'}; %controls, melancholia, nonmelancholic
nAUs=length(these_AUs);

dxmap=containers.Map({0,1,2},{'CN','MEL','NONMEL'}); %index 0 maps to 'controls', etc
invdxmap=containers.Map({'CN','MEL','NONMEL'},{1,2,3});

event_times=[125.6,125.6+168.84]; %the cutoff times (in sec) between the stimuli
laugh_times=[15,20,27,31,37,43,55,72,78,86,90,94,98,104,108,114,120]; %from my annotation, within BillCosby stimulus

format compact
%Following 3 files have the openface output .csvs, for all subjects, combined into a .mat file
%'data' field is a cell array of time series data. size n subjects * m stimuli
%'success' field is 'success' column of OpenFace output .csv
mat_bc=load('brisbane_controls.mat'); 
mat_sc=load('sydney_controls.mat');
mat_sp=load('sydney_patients.mat');

aulist=cellfun(@(x) x(3:end-2),mat_sc.aulist,'UniformOutput',false);
aunames=containers.Map({1,2,3,4,5,6,7,...
    8,9,10,11,12,13,...
    14,15,16,17},...
    {'Inner Brow Raiser','Outer Brow Raiser','Brow Lowerer','Upper Lid Raiser','Cheek Raiser','Lid Tightener','Nose Wrinkler'...
    'Upper Lip Raiser','Lip Corner Puller','Dimpler','Lip Corner Depressor','Chin Raiser','Lip Stretcher'...
    'Lip Tightener','Lip Parts','Jaw Drop','Blink'});
%%
%Preparation of data. Maps indices of aulist to AU names. Combines the
%controls and patients in the 'sydney patients' 'sydney control's and
%'brisbane controls' groups into 2 groups

stims_bc=[1,2,3,4]; %bc' represents 'brisbane controls'. indices of BillCosby, TheChamp, weather, CryFreedom
stims_s=[1,2,3,5]; %'s' represents Sydney. This code section is because the sydney and brisbane control subjects watched different videos (with some overlap)
stim_frames=containers.Map({'BillCosby','TheChamp','weather','CryFreedom'},[3140,4221,1443,3999]); %no of frames for each stimulus (excluding Interview)
stims=containers.Map({'BillCosby','TheChamp','weather','CryFreedom'},[1,2,3,4]); %code 1-4 for each stimulus
n_stims=length(stims_bc);

%Get single array for controls and patients. Remove data for unused
%stimuli

%Controls
c=struct();
c.data=[mat_sc.data(:,stims_s);mat_bc.data(:,stims_bc)]; 
c.dataPresent=[mat_sc.dataPresent(:,stims_s);mat_bc.dataPresent(:,stims_bc)]; %documents whether this participant watched this stimulus or not
c.nframes=[mat_sc.nframes(:,stims_s);mat_bc.nframes(:,stims_bc)]; %some subjects have 1 less or 1 too many frames
c.participants=[mat_sc.participants;mat_bc.participants];
c.stimuli=mat_sc.stimuli(stims_s);
c.success=[mat_sc.success(:,stims_s);mat_bc.success(:,stims_bc)]; %the 'success' column from the OpenFace .csv

%Patients
p=struct();
p.data=mat_sp.data(:,stims_s);
p.dataPresent=mat_sp.dataPresent(:,stims_s);
p.nframes=mat_sp.nframes(:,stims_s);
p.participants=mat_sp.participants;
p.stimuli=mat_sp.stimuli(stims_s);
p.success=mat_sp.success(:,stims_s);

%%
%Preparation of data
c.cubes=make_datacubes(c,these_AUs,aulist,stim_frames,these_stims); %each cube is nFrames x nAUs x nSubs Note that 'these_AUs' is applied in here so we're left with 14 AUs (excluding blink, etc)
p.cubes=make_datacubes(p,these_AUs,aulist,stim_frames,these_stims);
c.cube=join_cubes(c); 
%cube is shaped 8804 (nframes) * 14 (nAUs) * 54 (nsubjects), same as a_aus
%!!!
%join_cubes: if any participant has NaN in any video/AU, make that entire participant NaN for the cube
p.cube=join_cubes(p);
[c.cuben,c.medians,c.baselines]=normalize_entire_median(c);
[p.cuben,p.medians,p.baselines]=normalize_entire_median(p);

C_has_all_videos=(~isnan(c.medians(1,:)))';
P_has_all_videos=(~isnan(p.medians(1,:)))';

%{
%%
allclipsfile=fullfile('G:\My Drive\PhD\Project_Melancholia\CopiedData','RachelScott','FromRachelForSNG','Cleaned DATA stimuli and diagnosis - JAYSONSVERSION.xlsx');

%Read Cleaned DATA stimuli and diagnosis from RachelScott. Missing this
%data for Brisbane controls
T1=readtable(allclipsfile,'Sheet','Sydney_Patients');
IDs=cellfun(@(x) x(1:end-1),T1.DataPresent_1_yes,'UniformOutput',false);
inds=find_matching_inds(p.participants,IDs); %make sure no nans
T2=T1(inds,{'Gender','Diagnosis','Age','MIHistory','FamilyMIHistory','TakesAnyMeds','SSRI','Dual_actionADT','TCAOrMAOI','MoodStabilizer','Antipsychotic'});
T2.MIHistory(T2.MIHistory==99)=NaN; %make 99s into NaNs
T2.FamilyMIHistory(T2.FamilyMIHistory==99)=NaN;
p.T=T2;
mels=p.T.Diagnosis==1;
nonmels=p.T.Diagnosis==2;
%%
T1syd=readtable(allclipsfile,'Sheet','Sydney_controls','ReadVariableNames',true);
T1bris=readtable(allclipsfile,'Sheet','Brisbane_controls','ReadVariableNames',true);
variable_list={'DataPresent_1_yes','Age','Gender','Diagnosis','Other','MIHistory','FamilyMIHistory','TakesAnyMeds','SSRI','Dual_actionADT','TCAOrMAOI','MoodStabilizer','Antipsychotic'};
T1both=[T1syd(:,variable_list);T1bris(:,variable_list)];
IDs=cellfun(@(x) x(1:end-1),T1both.DataPresent_1_yes,'UniformOutput',false);
inds=find_matching_inds(c.participants,IDs); %make sure no nans
T2=make_genderDiagnosisTable(T1both,inds,{'Gender','Diagnosis','Age','MIHistory','FamilyMIHistory','TakesAnyMeds','SSRI','Dual_actionADT','TCAOrMAOI','MoodStabilizer','Antipsychotic'});
T2.MIHistory(T2.MIHistory==99)=NaN; %make 99s into NaNs
T2.FamilyMIHistory(T2.FamilyMIHistory==99)=NaN;
c.T=T2;
%}
%%
%Combine controls and participants 
a=struct(); %struct containing both
a.participants=[c.participants;p.participants];
a.cubes=cellfun(@(x,y) cat(3,x,y),c.cubes,p.cubes,'UniformOutput',false);
a.cube=cat(3,c.cube,p.cube); %8804 frames * 14 AUs * 120 participants
%a.cuben=cat(3,c.cuben,p.cuben);
%a.cubesn=cellfun(@(x) x-median(x,1),a.cubes,'UniformOutput',false);
a.medians=cat(2,c.medians,p.medians);
%a.baselines=cat(2,c.baselines,p.baselines);
a.data=[c.data;p.data];
a.dataPresent=[c.dataPresent;p.dataPresent];

load('diagnoses.mat'); %adds a variable called 'diagnoses' to the workspace
CN=diagnoses==0;
MEL=diagnoses==1;
NONMEL=diagnoses==2;

%%
%{
%Make a list of valid participants (a_valid).
%valid subjects (1s in a_valid) have all videos, don't have any action units with all zeros (causes HMM error), and are either CN or MEL
%}
Y=[CN,MEL,NONMEL];

A_has_all_videos=(~isnan(a.medians(1,:)))'; %whether subjects have all videos or not
allvars=sq(var(a.cube,1))==0;
A_nozerovars=all(~allvars,1)'; %pts with no AUs with all 0s

a_valid=A_has_all_videos&A_nozerovars&(CN|MEL); 

vY=Y(a_valid,:); %vY can be used to find indices of each group in the 'a_valid' subgroup
vCN=CN(a_valid); 
vMEL=MEL(a_valid);
vNONMEL=NONMEL(a_valid);

%{
diagnoses_valid=[];
for i=1:size(Y,1)
    diagnoses_valid(i)=find(Y(i,:))-1;
end
%}

%{
%only relevant for analyses involving non-melancholic participants
a_bigvalid=A_has_all_videos&(CN|MEL|NONMEL);
Ybig=Yold(a_bigvalid,:);
vbigCN=Ybig(:,1); vbigMEL=Ybig(:,2); vbigNONMEL=Ybig(:,3);
%}
%%

%Final preparation (actually important)
nframes=size(a.cube,1);
nsubs=length(diagnoses);
t = (0:nframes-1)/Fs; %timestamps

%clear unnecessary variables to reduce memory usage
clear acubet acubetblock allchunks allchunks2 cube cube2 cube_truncated;
clear IDs inds mat_bc mat_sc mat_sp means seq sexplained slatent smu sscore stsquared texplained tlatent tmu tscore ttsquared;
clear c p subdata T1 T1both T1bris T1syd T2 vars z;
clear acubes acubesblock mels
clear type support stims_s stims_bc s;
clear A_has_all_videos A_nozerovars allvars 
clear ans C_has_all_videos P_has_all_videos 

%{
%%
%Optional: Get demographics for Table 1

sum(CN&A_has_all_videos)
sum(MEL&A_has_all_videos)

age_CN=a.T{CN&A_has_all_videos,'Age'};
age_MEL=a.T{MEL&A_has_all_videos,'Age'};
gender_CN=a.T{CN&A_has_all_videos,'Gender'};
gender_MEL=a.T{MEL&A_has_all_videos,'Gender'};

sprintf('CN age: mean %.1f, std %.1f',nanmean(age_CN),nanstd(age_CN))
sprintf('MEL age: mean %.1f, std %.1f',nanmean(age_MEL),nanstd(age_MEL))
[h,p]=ttest2(age_CN,age_MEL); sprintf('Group difference in age mean: p= %.3f',p)

sprintf('CN: %i male, %i female',sum(gender_CN==1),sum(gender_CN==2))
sprintf('MEL: %i male, %i female',sum(gender_MEL==1),sum(gender_MEL==2))
eligible_inds=((a.T.Diagnosis==1 | a.T.Diagnosis==0) & A_has_all_videos);
[table,chi2,p,labels] = crosstab(a.T{eligible_inds,'Gender'},a.T{eligible_inds,'Diagnosis'});
sprintf('Group difference in gender: chi-squared p=%.3f',p)

clear age_CN age_MEL gender_CN gender_MEL h p eligible_inds table chi2 labels

sum(MEL & A_has_all_videos & a.T.TakesAnyMeds==1)
sum(MEL & A_has_all_videos & a.T.TakesAnyMeds==0)
sprintf('Total participants in melancholia without missing data: %i', sum(MEL & A_has_all_videos & (a.T.TakesAnyMeds==1 | a.T.TakesAnyMeds==0)))
sum(MEL & A_has_all_videos & a.T.SSRI==1)
sum(MEL & A_has_all_videos & a.T.Dual_actionADT==1)
sum(MEL & A_has_all_videos & a.T.TCAOrMAOI==1)
sum(MEL & A_has_all_videos & a.T.MoodStabilizer==1)
sum(MEL & A_has_all_videos & a.T.Antipsychotic==1)

sum(CN & A_has_all_videos & a.T.TakesAnyMeds==1)
sum(CN & A_has_all_videos & a.T.TakesAnyMeds==0)
sprintf('Total participants in controls without missing data: %i', sum(CN & A_has_all_videos & (a.T.TakesAnyMeds==1 | a.T.TakesAnyMeds==0)))
sum(CN & A_has_all_videos & a.T.SSRI==1)
sum(CN & A_has_all_videos & a.T.Dual_actionADT==1)
sum(CN & A_has_all_videos & a.T.TCAOrMAOI==1)
sum(CN & A_has_all_videos & a.T.MoodStabilizer==1)
sum(CN & A_has_all_videos & a.T.Antipsychotic==1)

%%
%Optional: Get neuropsych data for patient group. In control group only 8 subjects
%have neuropsych
neuropsych=fullfile('G:\My Drive\PhD\Project_Melancholia\CopiedData','RachelScott','FromRachelForSNG','neuropsych_questionnaires_feb2014 - JAYSONSVERSION.csv');
T1=readtable(neuropsych);
IDs=T1.ID;
L=length(a.participants);
has_neuropsych=logical(zeros(L,1)); %stores whether each pt has neuropsych or not
QIDS_total=NaN(L,1); %stores QIDS_total
CORE=NaN(L,4); 
%stores Ni (non-interactiveness), RT (retardation), AG (agitation), and final column is CORE_total
Qsort=NaN(L,32); %stores Qsort

Qsortlist={};
for i=1:32
    Qsortlist{i}=['Qsort_',int2str(i)];
end
for i=1:L
    participant=a.participants{i}(1:end-4);
    z= strcmp(IDs,participant) | strcmp(IDs,participant(4:end));
    if any(z) %If the participant has a row in neuropsych .csv file
        has_neuropsych(i)=true;
        index=find(z); %find the correct row in neuropsych .csv file
        assert(length(index)==1);
        
        QIDS_total_value=T1.QIDS_total(index); %get QIDS_total
        if ~isnan(QIDS_total_value) & QIDS_total_value<50 %value>50 is usually a NaN
            QIDS_total(i)=QIDS_total_value;
        end
        CORE_values=table2array(T1(index,{'CORE_NI','CORE_RT','CORE_AG','CORE_total'}));
        if ~any(isnan(CORE_values)) & all(CORE_values<50)
            CORE(i,:)=CORE_values;
        end
        Qsort_values=table2array(T1(index,Qsortlist));
        if ~strcmp(Qsort_values{1},'NA')
            Qsort_values=str2double(Qsort_values);
            if ~any(isnan(Qsort_values)) & all(Qsort_values<50)
                Qsort(i,:)=Qsort_values;
            end
        end    
    end
    
end
sprintf('%i participants have neuropsych out of %i in entire group',sum(has_neuropsych) , length(a.participants))
a.participants(~(has_neuropsych)); %print this to see which participants don't seem to have neuropsych
a.QIDS_total=QIDS_total;
a.CORE=CORE;
a.Qsort=Qsort;

%%
%Optional: Plot all time series for a particular AU and video
this_AU=3; 
adata=a.data(CN,stims('BillCosby')); %choose group and stimulus
alldata=cat(3,adata{:}); %remove dataNotPresent rows
alldata=squeeze(alldata(:,this_AU,:));


times=[3141:3141+4220,12803-3999:12803];
%[3140,4221,1443,3999]

alldata=squeeze(a.cube(times,this_AU,CN)); %uncomment this to use all videos concatenated

figure;
nsubs=width(alldata);
nplotsheight=ceil(sqrt(nsubs));
for i=1:nsubs
    subplot(nplotsheight,nplotsheight,i);
    plot(alldata(:,i));
    title(sprintf('%i',i));
end
%}
%%
%ANOVA part 1: group x video x AU
%%%%%You can also start here with load('all_variables.mat') %%%%

these_stims_ANOVA={'BillCosby','TheChamp'};
ANOVA_stims_for_plot={'Comedy','Sad Movie'}; %stimulus labels for the Figure
cubes=a.cubes;
temp={{'01','04','15'},{'06','12'}}; %Which AUs get summed for each emotion valence group
ANOVA_AUgroups={'Sadness','Happiness'}; %emotion labels for the Figure
these_AUs_ANOVA=cellfun(@(y) cell2mat(cellfun(@(x) find(strcmp(aulist(these_AUs),x)),y,'UniformOutput',false)), temp,'UniformOutput',false);
aulist_ANOVA=aulist(these_AUs(these_AUs_ANOVA{1})); %confirm correct AUs

%%
%ANOVA part 2: Prepare data for ANOVA of median activation for AUs
a_valid_ANOVA=(CN|MEL);
y=[]; G={};
for i=1:length(these_stims_ANOVA)
    stimname=these_stims_ANOVA{i};
    nstim=find(strcmp(these_stims,stimname));
    stimtype=ANOVA_stims_for_plot{i};
    cube=cubes{i};
    for nsub=1:nsubs
        if a_valid_ANOVA(nsub)
            diagnosis=diagnoses(nsub);
            diagnosis_name=these_diagnoses{diagnosis+1};
            matrix=cube(:,:,nsub);
            if a.dataPresent(nsub,nstim)  
                for nAUset=1:length(these_AUs_ANOVA) %for each AU;
                    temp=these_AUs_ANOVA{nAUset}; %list of AU id numbers, to add
                    ts=sum(matrix(:,temp),2);                   
                    y=[y;median(ts)]; %add median
                    AU=ANOVA_AUgroups{nAUset};
                    G=[G;{diagnosis_name,stimtype,AU}];
                end
                
            end
        end
    end
end

%%
%ANOVA part 3: Do ANOVA (only use terms with AU in it). 


these_terms=[0,0,1;1,0,1;0,1,1;1,1,1];
[pval,anovatbl,stats,terms]=anovan(y,{G(:,1),G(:,2),G(:,3)},'model',these_terms,'varnames',{'Diagnosis','Stim','AU'}); % This line makes result: A significant three-way interaction was found between clinical group...
%stimulus, and facial valence (p=0.003)"

figure;
results = multcompare(stats,'Dimension',[1,2,3]) %'results' array has p-values for post-hoc comparisons using Tukey's honestly significant difference criterion. %Also makes Figure 4 - supplement 1
%see documentation at https://www.mathworks.com/help/stats/multcompare.html

%%
%Prepare for Figure 4. get bootstrap versions of the mean at each
%timepoint. Read the corresponding description in the paper first, to
%understand this
CNboot=mybootmean(a.cube(:,:,CN),100);
MELboot=mybootmean(a.cube(:,:,MEL),100);

%%
%Prepare for Figure 4.
means_CN=mean(a.cube(:,:,CN),3,'omitnan');
means_MEL=mean(a.cube(:,:,MEL),3,'omitnan');

smwin=100; %smoothing factor for the confidence intervals
CNbootstd=std(CNboot,[],3);
MELbootstd=std(CNboot,[],3);

% 5% and 95% confidence bands
CNlowhigh={smoothdata(prctile(CNboot,5,3),'Gaussian',smwin),smoothdata(prctile(CNboot,95,3),'Gaussian',smwin);};
MELlowhigh={smoothdata(prctile(MELboot,5,3),'Gaussian',smwin),smoothdata(prctile(MELboot,95,3),'Gaussian',smwin);};

%Mean +- 1 stdev
%CNlowhigh={means_CN-CNbootstd, means_CN+CNbootstd}; %show mean +- one SD
%MELlowhigh={means_MEL-MELbootstd, means_MEL+MELbootstd};

%%
%Figure 4—figure supplement 3: Plot mean time series and shade confidence bands, for all AUs
clf;
for i=1:length(these_AUs)
    subplot(4,4,i);
    temp=aulist{these_AUs(i)};
    i=find(strcmp(aulist(these_AUs),temp)); 
    mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);
    xlabel('Time (s)'); ylabel('Intensity'); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
end
legend('Control','Melancholia')

%%
%Prepare for Figure 4. Principal component analysis
acubes=permute(a.cube(:,:,a_valid),[1,3,2]); %nframes * nsubs * nfeatures
%acubes=cfss7;%UNCOMMENT this line (data from melancholia_dynamics.m, to use time-freq bins as inputs
acubesblock=reshape(acubes,size(acubes,1)*size(acubes,2),size(acubes,3));%facial AU arrays of all subject's concatenated along time axis
addpath('D:\FORSTORAGE\ToolboxesMatlab\pca_ica\'); %From https://www.mathworks.com/matlabcentral/fileexchange/38300-pca-and-ica-package
[scoeff,~,slatent,~,~,~] = pca(acubesblock); playtone(); %spatial PCA 
means_CN_spca=squeeze(mean(acubes(:,vCN,:),2))*scoeff; %mean of controls' time series, in new basis (linear combination of AUs)
means_MEL_spca=squeeze(mean(acubes(:,vMEL,:),2))*scoeff; 

%%
%{
Prepare for Figure 4. Visualisation of spatial PCA component 1 with
makehuman face for lower right corner of Figure 4
%This code section saves AUs for the face at 'C:\Users\Jayson\Documents\makehuman\v1py3\data\facs\','melancholia_spca1.facs'.
%After this code section go to makehuman program and render the face. Save the result as 'C:\Users\Jayson\Documents\makehuman\v1py3\render\melancholia_spca_processed\melancholia_spca1.png'

Details: -->
If more than n AUs have contribution > 0, show all positive contributions
If less than n AUs have contributions > 0, show exactly the n largest
contributions
Contributions shown are rescaled to 0 to 1
%}
n=5; %at the least, show the n largest AU contributions (rescaled to 0 to 1)
temp=scoeff(:,1); %1st PCA component
temp_sorted=sort(temp,'descend');
temp=rescale(temp,'InputMin',min(temp_sorted(n+1),0));

%Now save data for FACSHuman in .json text format, but .facs file format
aulist2=cellfun(@(x) ['AU',int2str(eval(x))],aulist,'UniformOutput',false);
s=cell2struct(num2cell(temp),aulist2(these_AUs));
txt=jsonencode(s,'PrettyPrint',true);
txt=strrep(txt,"AU","");
DF=['C:\Users\Jayson\Documents\makehuman\v1py3\data\facs\','melancholia_spca1.facs'];
fid=fopen(DF,'w');
raw=fwrite(fid,txt);
fclose('all');

%%
%Generates Figure 4
means_these_AUs={'01','04','15','06','12'};
means_these_AUs_inds=cell2mat(cellfun(@(x) find(strcmp(aulist(these_AUs),x)),means_these_AUs,'UniformOutput',false));

clf;
p=panel();
p.pack({.3,.3,.4});
p(1).pack('h',3); %default {1/6, 1/3, 1/3, 1/6}
p(2).pack('h',3);
p(3).pack('h',{.6,.2,.2});

p(1,1).select(); temp='06'; i=find(strcmp(aulist(these_AUs),temp)); 
mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);
xlabel('Time (s)'); ylabel('Intensity'); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
p(1,2).select(); temp='12'; i=find(strcmp(aulist(these_AUs),temp)); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);

p(2,1).select(); temp='01';i=find(strcmp(aulist(these_AUs),temp)); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);

p(2,2).select(); temp='04';i=find(strcmp(aulist(these_AUs),temp)); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);
p(2,3).select(); temp='15';i=find(strcmp(aulist(these_AUs),temp)); xline(125.6); xline(125.6+168.84); set(gca,'xlim',[0,max(t)])
mytitle=(sprintf('AU%s: %s',aulist{these_AUs(i)},aunames(these_AUs(i)))); myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle);
%legend('Healthy control','Melancholic depression','Position',[0.8,0.8,0.1,0.1])



%For BillCosby, Display mean time series one group, and 1 (AU or SPCA component) together
%with laugh time annotations as vert lines
p(3,1).pack({.9,.1});
p(3,1,1).select();
cmap=brewermap(2,'RdBu');
thisblue=cmap(2,:);
thisred=cmap(1,:);
indices=1:3140; %indices corresponding to the stimulus
t_stim=t(indices);
ncomponent=1; %or nAU
values1=means_CN_spca(indices,ncomponent); %or means_CN_spca
values2=means_MEL_spca(indices,ncomponent);

%CNlowhigh needs to be 1 x 2 cell array, each of which is 3140 x 1 double

CNlowhigh_spca=cellfun(@(x) x(1:length(t_stim),:)*scoeff,CNlowhigh,'UniformOutput',false);
MELlowhigh_spca=cellfun(@(x) x(1:length(t_stim),:)*scoeff,MELlowhigh,'UniformOutput',false);

myplotmeans(t_stim,CNlowhigh_spca,MELlowhigh_spca,values1,values2,1,'');
xlabel('Time (s)'); ylabel('Intensity'); xlim([0,max(t_stim)]);
for i=1:length(laugh_times)
    xline(laugh_times(i));
end

p(3,2).select();
temp=imread(['C:\Users\Jayson\Documents\makehuman\v1py3\render\melancholia_spca_processed\melancholia_spca1.png']); %This is a face image generated by makehuman, from file 'melancholia_spca1.facs'
temp=temp(20:end,180:end-160,:);
imshow(temp);
xlim([0,size(temp,2)]); ylim([0,size(temp,1)]);

fontsize=10; 
annotation('rectangle',[.79,.2,0.025,0.07],'Color','none','FaceColor',thisblue);
annotation('rectangle',[.79,.1,0.025,0.07],'Color','none','FaceColor',thisred);
annotation('textbox',[.79+.02,.2+.01,0.02,0.05],'String','Control','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[.79+.02,.1+.01,0.02,0.05],'String','Melancholia','FitBoxToText','on','LineStyle','none','FontSize',fontsize);

fontsize=14;
annotation('textbox',[.02,1,0,0],'String','a','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[.02,.36,0,0],'String','b','FitBoxToText','on','LineStyle','none','FontSize',fontsize);

p.marginleft=10;
p.marginbottom=10;
p(1).marginbottom=15;
%p(1,2).marginleft=0; p(1,3).marginright=0;
p(3).marginbottom=0;
p.marginbottom=5;
p(2).marginbottom=4;
p(3).marginleft=5; p(3).marginright=5;

p.fontname='TimesNewRoman';
fontsize=8;

set(p(1,1).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);
set(p(1,2).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);
set(p(2,1).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);
set(p(2,2).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);
set(p(2,3).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);
set(p(3,1,1).axis,'FontSize',fontsize,'FontSizeMode','manual','YTick',[]);

%clear CNlowhigh CNbootstd MELlowhigh MELbootstd nplotsheight

colours=brewermap(3,'Accent');
thickness=.015;

%Dark2, Accent, *Dark2, *Set1, Set2, Set3, *Set 3
p(1,1).select(); mystimulusannotate(colours,thickness);
p(1,2).select(); mystimulusannotate(colours,thickness);
p(2,1).select(); mystimulusannotate(colours,thickness);
p(2,2).select(); mystimulusannotate(colours,thickness);
p(2,3).select(); mystimulusannotate(colours,thickness);
p(3,1,1).select(); mystimulusannotate(colours,thickness,1); %only comedy

fontsize=10; 
pos=[.69,.62]; boxsize=[.05,thickness]; textpos=[.05,-.01];
annotation('rectangle',[pos(1),pos(2)+.3,boxsize(1),boxsize(2)],'Color','none','FaceColor',colours(1,:));
annotation('rectangle',[pos(1),pos(2)+.2,boxsize(1),boxsize(2)],'Color','none','FaceColor',colours(2,:));
annotation('rectangle',[pos(1),pos(2)+.1,boxsize(1),boxsize(2)],'Color','none','FaceColor',colours(3,:));
annotation('textbox',[pos(1)+textpos(1),pos(2)+.3+textpos(2),0.02,0.05],'String','Stand-up comedy','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[pos(1)+textpos(1),pos(2)+.2+textpos(2),0.02,0.05],'String','Sad movie','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
annotation('textbox',[pos(1)+textpos(1),pos(2)+.1+textpos(2),0.02,0.05],'String','Weather report','FitBoxToText','on','LineStyle','none','FontSize',fontsize);
clear pos boxsize;


%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%

function [mcuben,medians,baseline]=normalize_entire_median(m)
medians=median(m.cube,'omitnan');
mcuben=m.cube-medians;
medians=squeeze(medians);
baseline=squeeze(median(m.cube(1:100,:,:),'omitnan')); %first 4 sec (assumes 25fps
end

function cubes=make_datacubes(m,these_AUs,aulist,stim_frames,these_stims)
%{
Make m.data into a 4D data-cube
Use stims in these_stims, AUs in these_AUs
Fill empty data with NaNs
%}
N=length(these_stims);
cube1=repmat({[]},size(m.data,1),N);
for sub=1:size(m.data,1)
    for n=1:N
        stimnum=find(strcmp(m.stimuli,[these_stims{n},'_wav'])); %get numeric code for the stimulus
        assert(logical(stimnum));
        stimframes=stim_frames(these_stims{n});
        data=m.data{sub,stimnum};
        if ~isempty(data)
            cube1{sub,n}=data(:,these_AUs);
        else
            cube1{sub,n}=NaN(stimframes,length(these_AUs));
        end
    end
end
cubes={};
for i=1:N
    stimframes=stim_frames(these_stims{i});
    cube1stimi=cell2mat(cube1(:,i)');
    cube_stimi=reshape(cube1stimi,stimframes,length(these_AUs),[]);
    cubes=[cubes,cube_stimi];
end
end

function mcube=join_cubes(m)
%join m.cubes into one cube
mcube=cat(1,m.cubes{:});
anyNANs=squeeze(any(isnan(mcube(:,1,:)),1)); %label time series with any NaNs
mcube(:,:,anyNANs)=NaN;
end

function Means=getmeans(m)
%Find mean activation
Means=[];
for i=1:length(m.cubes)
    Means=cat(1,Means,mean(m.cubes{i},3,'omitnan'));
end
end

function inds=find_matching_inds(X,Y)
%{
For each item in X, output index of matching item in cell array Y
X and Y are cell arrays of strings
If not found, use NaN
%}
inds=NaN(length(X),1);
for i=1:length(X)
    item=X{i};
    result=strcmp(Y,item);
    if sum(result)==1
        inds(i)=find(result);
    end 
end
end


function tbl=make_genderDiagnosisTable(T1,inds,params)
%{
Given table T1, get rows in index list inds, and columns in cell array
params
NaN items of inds become NaN rows of T1
%}
tbl=table();
for i=1:length(inds)
    ind=inds(i);
    if isnan(ind)
        tbl=[tbl;repmat({NaN},1,length(params))];
    else
        tbl=[tbl;T1(ind,params)];
    end
end
end

function out=mybootmean(array3,n_bootstrap)
    %Get bootstrap repeats of values in array3
    func=@(x) mean(x,'omitnan');
    out=NaN(size(array3,1),size(array3,2),n_bootstrap);
    for i=1:size(array3,2)
        i
        matrix=sq(array3(:,i,:));
        out(:,i,:)=bootstrp(n_bootstrap,func,matrix')';
    end
end

function []=myplotmeans(t,CNlowhigh,MELlowhigh,means_CN,means_MEL,i,mytitle)
cmap=brewermap(2,'RdBu');
thisblue=cmap(2,:);
thisred=cmap(1,:);
hold on;
plot(t,means_CN(:,i),'Color',thisblue); 
plot(t,means_MEL(:,i),'Color',thisred);
myshade(t,CNlowhigh{1}(:,i)',CNlowhigh{2}(:,i)',thisblue,0.3);
myshade(t,MELlowhigh{1}(:,i)',MELlowhigh{2}(:,i)',thisred,0.3);
hold off;
title(mytitle); %xlabel('Time (s)'); ylabel('Intensity');
max_xlim=floor(max(t)/10)*10;
set(gca,'XTick',[0,max_xlim]);
end

