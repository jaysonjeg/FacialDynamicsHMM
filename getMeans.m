function meansn=getMeans(hmm)
%Get means for all states of the hmm
means=[];
for k=1:hmm.K
    means=[means;getMean(hmm,k)];
end
multipliers=1./max(abs(means)')'; %rescale so that largest absolute value becomes + or - 1
meansn=means.*multipliers;
end

