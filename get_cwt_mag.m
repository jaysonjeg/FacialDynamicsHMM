function [out] = get_cwt_mag(values,cube)
%{
Given CWT image 'values' (complex), return either amplitude or power
Optionally normalise by std or variance of input array

Inputs: 
values: nfrqs x ntimes x nsubs x nAUs
%}

square_amplitude=false; %True: Use magnitude (abs(value)).^2 instead of amplitude (abs(value))
normalize_by_variance=false; %normalize by variability or not. Setting to true makes no difference, because HMM normalizes anyway

out=abs(values);
if square_amplitude
    out=out.^2;
end

if normalize_by_variance & nargin==2
    if square_amplitude   
        variability=var(cube,1); 
    else
        variability=std(cube,1);
    end
        
    variability=reshape(variability,[size(variability,3),size(variability,2)]);
    out=out./permute(repmat(variability,1,1,size(values,1),size(values,2)),[3,4,1,2]); 
end

end

