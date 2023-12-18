function [mu,variance] = mean_mag(beta,kmax)
% AMSC 660 Final Problem 2
% Write a code to compute the mean magnetization mu(beta) for the 2D Ising
% model by the Metropolis algorithm. Use a 30x30 lattice with periodic
% boundary conditions. Your program should compute the running mean and
% running variance for the magnetization. 

%% Initialization
N = 30;                 % dimension of the square lattice
S = ones(30,30);        % initialize all spins up, so s_ij = 1 for all i,j

% find the initial magnetization of the initial state
m = magnetiz(S,N);
mu = m;
variance = 0;

%% Calculating mean magnetization loop
for k = 1:kmax
    % first, pick a random site (i,j) and propose to flip the spin at it,
    % and calculate the energy difference delta_H. 
    rand_site = randi(N,2,1);
    i = rand_site(1);
    j = rand_site(2);
    delta_H = energy_diff(S,i,j,N);

    if delta_H < 0
        accept = true;
    else
        u = rand;
        if u < exp(-beta*delta_H)
            accept = true;
        else
            accept = false;
        end
    end

    if accept == true
        % flip the spin we proposed to flip, then calculate the
        % magnetization of the new state
        S(i,j) = -S(i,j);
        m = magnetiz(S,N);
    end

    % update the running mean magnetization and running variance
    mu = (k*mu + m) / (k+1);
    variance = ((k-1)*variance + (m - mu)^2) / k;

    if mod(k,1e7) == 0
        fprintf('beta = %d, iter = %d\n',beta,k);
    end
end
end


%% functions
function m = magnetiz(S,N)
% calculates the magnetization of the state S
m = (1/N^2) * sum(S,'all');
end

function delta_H = energy_diff(S,i,j,N)
% this function will calculate delta_H (the difference in energy) between
% the current state and the proposed state with spin flip at random_site. 

% find nearest neighbor indices
i_minus = mod(i-1,N);
i_plus = mod(i+1,N);
j_minus = mod(j-1,N);
j_plus = mod(j+1,N);

% taking modulus will screw up the last column and top row (indexed at N)
% of S and make the index 0 (can't be done in Matlab), so if the modulus is
% 0, we add N (=30) to it
    if i_minus == 0
        i_minus = i_minus + N;
    end
    if i_plus == 0
        i_plus = i_plus + N;
    end
    if j_minus == 0
        j_minus = j_minus + N;
    end
    if j_plus == 0
        j_plus = j_plus + N;
    end

% calculate delta_H
delta_H = 2*S(i,j) * (S(i_minus,j)+S(i_plus,j)+S(i,j_minus)+S(i,j_plus));
end
