function [xx] = leverage_score_sample_with_replacement(Q,k)
% Leverage score sampling

lev = vecnorm(Q').^2;

pp = cumsum(lev);
pp = pp/pp(end); % pp = pp/k


xx = rand(k,1);
for ii = 1:numel(xx)
    xx(ii) = min(find(xx(ii)<=pp,1));
end
%xx = unique(xx);
end