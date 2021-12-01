function [windowsize_sec,smwin,configurations]=get_defaults_HMM(number_states)
%Get default parameters for HMM and subsequent analysis

windowsize_sec=8;  %window size for averaging, for plotting common Viterbi path, and for consistency measures
smwin=1; %Gaussian smoothing window for line plots. Default 1 = no smoothing

template_configuration = struct();
template_configuration.order = 0; 
template_configuration.dropstates = 1; %Drop unnecessary states
template_configuration.verbose = 1;
template_configuration.cyc = 500;
template_configuration.initcyc = 10;
template_configuration.covtype="full";
template_configuration.pca=10;
template_configuration.downsample=10; %Downsample time series, speeds up inference. Default 10Hz
template_configuration.DirichletDiag=1;

%Alternative options which can be uncommented
%template_configuration.filter=[0,3];
%template_configuration.zeromean = 1; 
%template_configuration.standardise=0;
%template_configuration.BIGNbatch=5;

i = 1; 
configurations = {}; 
for k = number_states
    configurations{i} = template_configuration;
    configurations{i}.K = k;
    i = i + 1;
end

end

