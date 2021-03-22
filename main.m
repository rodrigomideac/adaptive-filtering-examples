close all
%% Simulation properties
% The simulation will take the expected value of MSE, MSD, and EMSE
% by performing many realizations.
N = 600;                                  % Total samples
M = 10;                                     % Filter length
R = 300;                                    % Realizations Qty

%% Optimum filter definition
% A FIR filter needs to be used as optimum plant.
% To use MATLAB filter function, we need to pass the poles vector. Passing 
% [1] just means the optimum filter doesn't have poles.
A = [1];                                    
B = randn(M,1);                             % Create a random FIR filter
M = 4
B = [1 0.5 -1 2]
sig2v = 1e-2;                               % Noise level

%% Other definitons
mu = 0.01;                                  % Adaptive filter step size
mu_nlms = 0.2 %2 * mu*M / (2+ mu*M)

%% Global Accumulators
% These accumulators are responsible for storing desired metrics between
% realizations
% MSE: Mean Square Error
% EMSE: Excess Mean Square Error
% MSD: Mean Square Deviation

% FIR LMS Adaptive Filter accumulators
mse_lms = zeros(N,1);
emse_lms = zeros(N,1);
msd_lms = zeros(N,1);

% FIR NLMS Adaptive Filter accumulators
mse_nlms = zeros(N,1);
emse_nlms = zeros(N,1);
msd_nlms = zeros(N,1);

%% Run realizations
alpha_ = 0.0;
beta_ = 1- alpha_;

parfor r=1:R
    %% Signal definitions
    % Define input signal
    u = randn(N,1);    
    u = get_u_autoregressive_process_given_alpha_beta(N,1,alpha_,beta_)';
    
    % Define noise signal as white, with variance sig2v
    v = sqrt(sig2v)*randn(N,1);
    
    % Run input signal through optimum filter. This will be used in the
    % filters to measure EMSE
    yo = filter(B,A,u);
    
    % Define desired signal by adding noise to the optimum filter output
    d = yo + v;
    
    %% Estimate Covariance Matrix
    % The matrix covariance matrix Ru will later on be used to calculate
    % the theoretical EMSE. It can be estimated by taking the mean of
    % outer product of a sliding window with size M of the input signal.
    ru = zeros(M,M);
    u_reg = zeros(1,M);
    for n = 1:N
        u_reg =[u(n) u_reg(1:M-1)];
        ru = ru + 1/N .* (u_reg' * u_reg);
    end
    
    mu_nlms = 2 * mu * trace(ru) / (2+ mu * trace(ru));

    %% Run the LMS Filter
    [r_msd,r_mse,r_emse,w] = fir_lms(u,d,yo,B,M,N,mu,false);
    mse_lms = mse_lms + (1/R) .* r_mse;
    emse_lms = emse_lms + (1/R) .* r_emse;
    msd_lms = msd_lms + (1/R) .* r_msd;
    
    %% Run the NLMS Filter
    [r_msd,r_mse,r_emse,w] = fir_lms(u,d,yo,B,M,N,mu_nlms,true);
    mse_nlms = mse_nlms + (1/R) .* r_mse;
    emse_nlms = emse_nlms + (1/R) .* r_emse;
    msd_nlms = msd_nlms + (1/R) .* r_msd;
end

%% Asserts
% It is very easy to make a mistake with the accumulator vector dimensions.
% If that was the case, the accumulators would, at this point, be matrices
% instead of vectors. The assertions just make sure that all accumulators
% are actually column vectors with size N.
assertDimensions(mse_lms,zeros(N,1));
assertDimensions(emse_lms,zeros(N,1));
assertDimensions(msd_lms,zeros(N,1));
assertDimensions(mse_nlms,zeros(N,1));
assertDimensions(emse_nlms,zeros(N,1));
assertDimensions(msd_nlms,zeros(N,1));

%% Plot MSE
fig = figure;
hold all;
plot(10*log10(mse_lms),'DisplayName','LMS');
plot(10*log10(mse_nlms),'DisplayName','NLMS');
title('MSE');
xlabel('n')
ylabel('MSE (dB)')
legend('show')
y = 10*log10(sig2v);
line([0,N],[y,y],'Color','r','DisplayName','\sigma_v^2 (noise level)');
saveas(fig,'./img/mse.png')

%% Plot EMSE
fig = figure;
hold all;
plot(10*log10(emse_lms),'DisplayName','LMS');
plot(10*log10(emse_nlms),'DisplayName','NLMS');
title('EMSE');
xlabel('n')
ylabel('EMSE (dB)')
legend('show')
lms_theoretical_emse = 10*log10(mu*sig2v*trace(ru)/2);
line([0,N],[lms_theoretical_emse ,lms_theoretical_emse ],'Color','r','DisplayName','Theoretical EMSE');
saveas(fig,'./img/emse.png')

%% Plot MSD
fig = figure;
hold all;
plot(10*log10(msd_lms),'DisplayName','LMS');
plot(10*log10(msd_nlms),'DisplayName','NLMS');
title('MSD');
xlabel('n')
ylabel('MSD (dB)')
legend('show')
lms_theoretical_msd = 10*log10(mu*sig2v*M/2);
line([0,N],[lms_theoretical_msd ,lms_theoretical_msd ],'Color','r','DisplayName','Theoretical MSD');
saveas(fig,'./img/msd.png')

%%
% This function calculates a matrix U with being composed as white noise
% where can be used as input u at an adaptive filter.
% Number of lines of U and d are determined by k

% Arguments:
%   N: How many samples will U and d hold (scalar)
%   k: number of nodes (scalar)
%   M: coeffient vector size (scalar)
%   w_o: M x 1 optimum coefficients vector

% Returns:
%   U: k x N input signal matrix
%   d: k x N desired signal

%%
function [U] = get_u_autoregressive_process_given_alpha_beta(N,k,alpha,beta)
    assertDimensions(alpha,ones(k,1));
    z = randn(k,N);
    U = zeros(k,N);
    
    for node=1:k
        for i=2:N
            U(node,i) = alpha(node)*U(node,i-1) + beta(node)*z(node,i);
        end
    end
  
end
