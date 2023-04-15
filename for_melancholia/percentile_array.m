function [array_perc] = percentile_array(array,array_surrogates)
%{
INPUTS: array is M x N, array_surrogates is M x N x Nsurrogates
OUTPUTS: array_perc is M x N. Has the percentile of each value in array,
relative to the distribution in array_surrogates
%}
array_perc=array;  
Nsurrogates=size(array_surrogates,3);
diff2=repmat(array,1,1,Nsurrogates); %make it same size as surrogatesAU
s2=sort(array_surrogates,3);
x1=diff2<s2; 
for i=1:size(array,1)
    for j=1:size(array,2)
        index=find(x1(i,j,:),1);
        if isempty(index)
            value=Nsurrogates;
        else
            value=find(x1(i,j,:),1);
        end
        value=value*100/Nsurrogates; %convert to percentile 0-100
        array_perc(i,j)= value;
    end
end
end

