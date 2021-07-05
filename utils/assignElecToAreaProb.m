function elec_area = assignElecToAreaProb(group_prob)
    % sample a random number between 0 and 1 
    %rng('shuffle');
    p = rand(size(group_prob,1),1); % * ones(1,size(gp,2));
    % compute cumulative sum
    a = cumsum(group_prob,2) > p;
    sa = sum(a,2) < 1;
    [~,elec_area] = max(a,[],2);
    elec_area(sa) = nan;
end