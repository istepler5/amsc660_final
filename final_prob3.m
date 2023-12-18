% AMSC 660 Final Problem 3

% Use the mean_mag function code to calculate mean magnetization for the
% set of values beta = 0.2:0.01:1. Initially, set all spins up for each
% value of beta. Desirably, use kmax = 10^8 iterations for each value of
% beta. Plot the graph of the compute mean magnetization as a function of
% beta and superimpose it with the graph of the analytic expression for
% mean magnetization. Also, plot mean +/- sqrt(Var(m)). 

beta = 0.2:0.01:1;
kmax = 1e8;
n = length(beta);
mu = zeros(n,1);
variance = zeros(n,1);
analytic_mu = zeros(n,1);

% all spins are set up in the initialization for mean_mag function
for i = 1:n
    [mu(i),variance(i)] = mean_mag(beta(i),kmax);
    if beta(i) > 0.4408
        analytic_mu(i) = (1 - (sinh(2*beta(i)))^(-4) )^(1/8);
    else
        analytic_mu(i) = 0;
    end
end

figure;
grid on;
plot(beta,analytic_mu,'b-');
hold on;
plot(beta,mu,'r.');
plot(beta,mu+sqrt(variance),'r--');
plot(beta,mu-sqrt(variance),'r--');
hold off;
xlabel('beta'); ylabel('mean magnetization');


