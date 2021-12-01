function [adj] = network_cutoffpercentile(adj,prc)
%Given adjacency matrix adj, convert edges less than percentile 'prc' to
%zeros
    cutoff=prctile(adj,prc);
    adj(adj<cutoff)=0; %don't visualise

end

