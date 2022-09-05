%{
This code calls OpenFace version 2.2.0 (https://github.com/TadasBaltrusaitis/OpenFace) to extract
 action unit time series from DISFA videos (http://mohammadmahoor.com/disfa/)
Before running this code, start Matlab at current folder:
OpenFace-master/OpenFace-master/MyResults

If your data is already in a 3D array (ntimepoints x nAUs x nSubjects), you
can go straight to DISFA_Step3_cubes2cwt2hmm.m
%}
format compact
in_folder='D:\FORSTORAGE\EmotionDatabases\DISFA\Videos_LeftCamera'; %MODIFY THIS: DISFA video files are here
out_folder = 'D:\FORSTORAGE\Data\Melancholia\temp'; %MODIFY THIS. Output will be saved here
executable = '"D:/FORSTORAGE/OpenFace-master/OpenFace-master/x64/Release/FeatureExtraction.exe"'; %points to the executable in OpenFace-master folder

%Get list of folders and files
%in_folder and in_files will store the folders and .avi files respectively
dir2=dir(in_folder);
in_files={dir2.name}; in_files=in_files(3:end); 

%run feature extraction using the produced folder and file lists, and save
for i=1:length(in_files) %i indexes folders in in_folders
    i,
    command=sprintf('%s -f "%s" -au_static -out_dir "%s" -aus', executable, fullfile(in_folder,in_files{i}), fullfile(out_folder, in_files{i}));
    %This command is based on https://github.com/TadasBaltrusaitis/OpenFace/wiki/Command-line-arguments
    dos(command);
end
