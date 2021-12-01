function [vconsistency, vmostcommon,vconsistencies] = get_consistency_measures(windowsize,vpathx,numstates,getnull)
%{
Get measures of between-subject consistency
INPUTS:
    window size (number of frames) in each direction, Default=0 look at single time point
    vpathx contains Viterbi paths
    numstates is number of HMM states = h1.K
    getnull is optional. If it's true, calculate and return
    vconsistencies. vconsistencies is used in Figure 2, Panel D, blue
    shading
%}
    
    num_of_values=size(vpathx,2); 
    nframes=size(vpathx,1);
    vmostcommon=NaN(nframes,1); %vector of most common state
    vconsistency=NaN(nframes,1); %vector of % of subjects having this state
    
    func= @(matrix) max(arrayfun(@(n) sum(any(matrix'==n,1)), 1:numstates))/num_of_values;
    n_bootstrap=100;
    vconsistencies=NaN(nframes,n_bootstrap);%stores bootstrap estimates of vconsistency
     
    for i=1:nframes
        if i-windowsize>=1 & i+windowsize<=nframes
            xs=vpathx(i-windowsize:i+windowsize,:);
            counts=arrayfun(@(n) sum(any(xs==n,1)), 1:numstates);
            %for each HMM state, how many subjects expressed it at least once in this
            %timewindow  
            [maxval,maxind]=max(counts);    
            vmostcommon(i)=maxind;
            vconsistency(i)=maxval/num_of_values;

            if nargin==4 & getnull %calculate bootstrap vconsistencies 
                bootstat=bootstrp(n_bootstrap,func, xs');
                vconsistencies(i,:)=bootstat';
            end
        end
    end
    vmostcommon=vmostcommon';
    vmostcommon(1:numstates)=1:numstates; %to ensure all states are present in time series
end