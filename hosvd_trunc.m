function [core, F] = hosvd_trunc(core, F, tol)
% Recompress core and factors to tolerance TOL using deterministic STHOSVD

d = ndims(core);
F_QR = cell(1,d); R_QR = cell(1,d);
for j = 1:d
    [F_QR{j}, R_QR{j}] = qr(F{j},0);
end
core = pagetmprod(core, R_QR, [1:d]);
% order does not matter a lot here as factors are already not too big

sz_new = size(core);
[FF2,core, sv] = mlsvd(core,sz_new,[1:d]); % no truncation here
for j = 1:d
    sv{j} = sv{j} / sv{j}(1); % make them relative preparing for precise rank-estimation
end     
kk = zeros(d,1);
for j = 1:d
    kk(j) = find(sv{j} > tol, 1, 'last');   % find where to truncate
end
% HOSVD truncation of core takes place here:
if d == 3
    core = core(1:kk(1), 1:kk(2), 1:kk(3));
elseif d == 4
    core = core(1:kk(1), 1:kk(2), 1:kk(3), 1:kk(4));
elseif d == 5
    core = core(1:kk(1), 1:kk(2), 1:kk(3), 1:kk(4), 1:kk(5));
end
% and here:
for j = 1:d
        FF2{j} = FF2{j}(:,1:kk(j));
end

% HOSVD truncation of factors takes place here:
for j = 1:d
    F{j} = F_QR{j}*FF2{j};
end

end