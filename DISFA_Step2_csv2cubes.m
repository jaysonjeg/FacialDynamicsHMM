%{
Given .csv OpenFace results, get the AU time series into nice arrays
%}

%From OpenFace outputs, extract table of facial AUs over time
in_folder='D:\FORSTORAGE\OpenFace-master\OpenFace-master\MyResults\DISFA_au_static\'; %MODIFY THIS
dir2=dir(in_folder);
in_filenames={dir2.name}; in_filenames=in_filenames(3:end); 
t_aus={};successes={};  %to store data

%{
t_aus will store AUs for each subject as tables
succeses will have the sucess column for each subject
%}

for i=1:length(in_filenames) %for each subject's file
    in_filename=in_filenames{i};
    T1=readtable(strcat(in_folder,in_filename,'/',in_filename(1:end-4),'.csv')); %read csv file
    t_aus{i}=T1(:,{'AU01_r','AU02_r','AU04_r','AU05_r','AU06_r','AU07_r','AU09_r','AU10_r','AU12_r','AU14_r','AU15_r','AU17_r','AU20_r','AU23_r','AU25_r','AU26_r','AU45_r'});
    successes(:,end+1)={table2array(T1(:,'success'))};
end
clear T1 in_filename in_folder dir2 in_filenames; 

%%
%{
DISFA subjects have either 4854 or 4855 frames in their
videos. We pad the end of the shorter videos with last frame repeated.
%}
maxlen = max(cellfun(@(x) length(x), successes)); %find maximum no of frames across all videos
successes = cell2mat(cellfun(@(a) cat(1,a,zeros(maxlen-length(a),1)), successes, 'UniformOutput',0)); %pad 'successes' with zeros
t_aus=cellfun(@(t)[t;repmat(tail(t,1),maxlen-size(t,1),1)] , t_aus, 'UniformOutput',0); %pad 't_aus' with last frame repeated
clear maxlen;

%%
%Deal with missing frames
missing_successes=sum(successes,2) < size(successes,2); %list of logicals for each frame. 1 if that frame is missing in ANY subject, 0 otherwise
sprintf('percentMissingFrames=%.3f percent', 100*sum(missing_successes)/size(successes,1)) %Display what % of frames are missing in ANY subject?
%imagesc(successes); xlabel('subject'); ylabel('frame'); %shows missing frames in more detail, for each subject

%Remove missing frames (only ~0.5% of frames in DISFA)
for subject=1:length(t_aus)
    temp=t_aus{subject};
    temp(missing_successes,:)=[];
    t_aus{subject}=temp;
end

%%
%Convert (cell array of tables) to (cell array of arrays)
a_aus={}; 
for subject=1:length(t_aus)
    a_aus{subject}=table2array(t_aus{subject});
end

%Convert to 3D array
nsubs=length(a_aus); %number of subjects
nAUs=size(a_aus{1},2);
nframes=size(a_aus{1},1);
a_aus=reshape(cell2mat(a_aus),nframes,nAUs,[]); %Frames x AUs x Subjects array, goes into next piece of code

clear missing_successes subject successes t_aus temp;
