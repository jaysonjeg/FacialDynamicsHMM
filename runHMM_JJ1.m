function [hmm,Gamma,Gammaup,Xi,vpath,fehist,maxFO,newframerate] = runHMM_JJ1(configurations,X,T)

%Make sure initGamma_random is set to rng(1) for reproducibility

t = (0:T(1)-1)/configurations{1}.Fs;
L = length(configurations);
Gamma = cell(length(configurations),1);
hmm={}; Xi={}; vpath={};fehist={}; Gammaup={}; maxFO={};

tic
for i = 1:length(configurations)
    tic
    disp([num2str(i) ' of ' num2str(L)])
    [hmm{i},Gamma{i},Xi{i},vpath{i},~,~,fehist{i}] = hmmmar(X,T,configurations{i});
    if isfield(configurations{i},'downsample')
        Gammaup{i}=padGamma(Gamma{i},T,configurations{i}); %upsampled Gamma
        tdown=downsample(t,ceil(configurations{i}.Fs/configurations{i}.downsample)); %common to all configurations
    else
        Gammaup{i}=Gamma{i};
    end
    maxFO{i} = getMaxFractionalOccupancy(Gammaup{i},T,configurations{i});
    toc
    
end
toc

if isfield(configurations{1},'downsample')
    newframerate=configurations{1}.downsample;
else
    newframerate=Fs;
end

playtone();
end

