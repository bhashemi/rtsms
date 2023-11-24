% Adaptive rank Examples 1-3 in 
% B. Hashemi and Y. Nakatsukasa, RTSMS: Randomized Tucker with single-mode
% sketching, submitted, 2023.

clc, clear, close all

%% Generate low-rank tensor A of size n(1) x n(2) x n(3)

% Example 1
%%{
n = [600 600 600];
f = @(x,y,z) 1./(5 + x.^2 + y.^2 + z.^2); 
%%}

%%
% Example 2
%{
n = [800 1200 300]; 
% Wagon's function from Page 99 of F. Bornemann, D. Laurie, S. Wagon and J. 
% Waldvogel, The SIAM 100-Digit Challenge, SIAM, 2004.
f = @(x,y,z) exp(sin(50*x)) + sin(60*exp(y)).*sin(60*z) + ...
             sin(70*sin(x)).*cos(10*z) + sin(sin(80*y)) - ...
             sin(10*(x+z)) + (x.^2 + y.^2 + z.^2)/4;
%}

%%
% Example 3
%{
n = [1000 1000 1000];
% cooked up by Nick Trefethen 16 May 2016
f = @(x,y,z) sqrt(x.^2+y.^2+z.^2);
%}

%%
% sample points in 1-dimension (Chebyshev pts in [-1, 1])
m = n - 1;
x = sin(pi*(-m(1):2:m(1))/(2*m(1))).';  % (Use of sine enforces symmetry.)
y = sin(pi*(-m(2):2:m(2))/(2*m(2))).';  
z = sin(pi*(-m(3):2:m(3))/(2*m(3))).';  
% sample points in d-dimensions:
[xx,yy,zz] = ndgrid(x,y,z);

A = f(xx,yy,zz); 
clear xx yy zz f n
d = ndims(A);
Frob_A = frob(A);  % used several times in residual computations

%% Set parameters
order = [3 2 1];               % processing order
k = 4;                         % sketch parameter
all_tols = logspace(-2,-14,7); % input tolerances to try 1e-02   1e-04 ... 1e-14
T = 5;                         % repeat T times for each tolerance to take average
%T = 2; 

%%
fprintf('\n--- RTSMS ---\n')

%HOSVD_truncation = false;
HOSVD_truncation = true;      % try RHOSVDSMS as well
rank_rtsms = cell(size(all_tols)); rank_rhosvdsms = cell(size(all_tols));

for i = 1:numel(all_tols)
    i
    tol = all_tols(i)

    for j = 1:T        
        rtsms(A, k, order, tol); % 1st time slower. Ignore it
        j
        tic        
        [core_rtsms, F_rtsms] = rtsms(A, k, order, tol);
        t_rtsms(i,j) = toc;
        res_rtsms(i,j) = frob(A - pagetmprod(core_rtsms, F_rtsms, [1:d]))/Frob_A;
        rank_rtsms{i}(j,:) = size(core_rtsms);

        if HOSVD_truncation
            hosvd_trunc(core_rtsms, F_rtsms, tol);            
            tic
            [core_rhosvdsms, F_rhosvdsms] = hosvd_trunc(core_rtsms, F_rtsms, tol); % renamed best code
            t_rhosvdsms(i,j) = t_rtsms(i,j) + toc;
            res_rhosvdsms(i,j) = frob(A - pagetmprod(core_rhosvdsms, F_rhosvdsms, [1:d]))/Frob_A;
            rank_rhosvdsms{i}(j,:) = size(core_rhosvdsms);
        end
    end
    
end
avg_res_rtsms = geomean(res_rtsms')
avg_t_rtsms = mean(t_rtsms')
std_rtsms1 = std(res_rtsms');
std_rtsms2 = std(-log10(res_rtsms'));
for jj = 1:numel(rank_rtsms)
    avg_rank_rtsms(jj) = round(mean(mean(rank_rtsms{jj})));
end

if HOSVD_truncation
    avg_res_rhosvdsms = geomean(res_rhosvdsms')
    avg_t_rhosvdsms = mean(t_rhosvdsms')
    std_rhosvdsms1 = std(res_rhosvdsms');
    std_rhosvdsms2 = std(-log10(res_rhosvdsms'));
    for jj = 1:numel(rank_rtsms)
        avg_rank_rhosvdsms(jj) = round(mean(mean(rank_rhosvdsms{jj})));
    end
end

%% Plot the results
MS = 'markersize'; ms = 14;
FS = 'fontsize'; fs = 15;
LW = 'linewidth'; lw = 2;

subplot(131)
p1 = loglog(all_tols, avg_res_rtsms,'o-b', MS, ms, LW, lw);
hold on
loglog(all_tols, avg_res_rhosvdsms,'x:', MS, ms, 'Color',[0.4 0.4 0.4], LW, lw)

set(gca, 'XDir','reverse')    % from large tol to small tol is better
grid on 
set(gca, 'MinorGridAlpha', 0) % too many grid lines; Let some of them disappear!
%xlabel('input tolerance')
title('computed residual')
xt = 10.^(-14:2:-2); xticks(xt), yticks(xt), xlim([5e-15 3e-2])
ylim([min(avg_res_rtsms)/5 max(avg_res_rtsms)*5])
axis padded
set(gca, FS, fs)

subplot(132)
p1 = semilogx(all_tols, avg_rank_rtsms,'o-b', MS, ms, LW, lw);
hold on
semilogx(all_tols, avg_rank_rhosvdsms,'x:', MS, ms, 'Color',[0.4 0.4 0.4], LW, lw)
set(gca, 'XDir','reverse') 
grid on
set(gca, 'MinorGridAlpha', 0) 
title('computed ranks')
xt = 10.^(-14:2:-2); xticks(xt), xlim([5e-15 3e-2]), 
ylim([1 max(avg_rank_rtsms)+1])
set(gca, FS, fs)

subplot(133)
semilogx(all_tols, avg_t_rtsms,'o-b', MS, ms, LW, lw)
hold on
semilogx(all_tols, avg_t_rhosvdsms, 'x:', MS, ms, 'Color',[0.4 0.4 0.4], LW, lw)

leg = legend({'RTSMS','RHOSVDSMS'},'Location','nw');
set(leg, FS, fs-6)
set(gca, 'XDir','reverse')
grid on, set(gca, 'MinorGridAlpha', 0)
title('time (sec)')
xticks(xt), xlim([5e-15 3e-2])
set(gca, FS, fs)
shg