% Fixed-rank Example 6 in 
% B. Hashemi and Y. Nakatsukasa, RTSMS: Randomized Tucker with single-mode
% sketching, submitted, 2023.

clc, clear, close all

%% Generate 4D Hilbert tensor
d = 4;
n = 150*ones(1,d);
[i1,i2,i3,i4] = ndgrid(1:n(1));
A = 1./(i1+i2+i3+i4-(d-1));
clear i1 i2 i3 i4
Frob_A = frob(A);

%% Set parameters
order = 1:d;  % processing order
num_ex = 6;   % number of input ranks to try
T = 5;        % repeat T times for each rank to take average
%T = 2; 

%% Call RTSMS and RHOSVDSMS
clc
k = 4;
for j = 1:T
    j
    for i = 1:num_ex
        i
        r = 5*i*ones(1,d);
    
        tic
        [core, F] = rtsms(A, k, order, r, 'rank');
        t_rtsms(i,j) = toc;
                
        [C_rhosvdsms, F_rhosvdsms] = hosvd_trunc_rank(core, F, r);
        t_rhosvdsms(i,j) = toc;

    res_rhosvdsms(i,j) = frob(A - pagetmprod(C_rhosvdsms, F_rhosvdsms, [1:d]))/Frob_A
    res_rtsms(i,j) = frob(A - pagetmprod(core, F, [1:d]))/Frob_A
    end

end

%% Take average
avg_res_rtsms = geomean(res_rtsms')
avg_res_rhosvdsms = geomean(res_rhosvdsms')
avg_t_rtsms = mean(t_rtsms')
avg_t_rhosvdsms = mean(t_rhosvdsms')
avg_res_rtsms
avg_t_rtsms

%%
MS = 'markersize'; ms = 14;
FS = 'fontsize'; fs = 15;
LW = 'linewidth'; lw = 2;

figure
subplot(121)
R = 5*(1:num_ex);

semilogy(R, avg_res_rtsms,'o-b',MS, ms, LW, lw)
hold on
semilogy(R, avg_res_rhosvdsms,'x-.', 'Color',[0.4 0.4 0.4], MS, ms, LW, lw);
legend('RTSMS','RHOSVDSMS')
grid on
axis padded
xlabel('rank')
xticks(R)
set(gca,'MinorGridAlpha', 0) % too many grid lines; Let some of them disappear!
title('computed residual')
set(gca, FS, fs)

subplot(122)
plot(R, avg_t_rtsms,'o-b',MS, ms, LW, lw), hold on
plot(R, avg_t_rhosvdsms,'x-.', 'Color',[0.4 0.4 0.4], MS, ms, LW, lw);
grid on
axis padded
xlabel('rank')
xticks(R)
title('time (sec)')
set(gca, FS, fs)
shg