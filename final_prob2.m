% AMSC 660 Final Project #2

% Compute positions of vertices of a graph with adjacency matrix given by
% the file Adjacency_matrix.csv (N=91 vertices). Choose initial positions
% for vertices for example as Gaussian random variables with mean zero and
% standard deviation N. Then use an appropriate deterministic optimizer to
% optimize the positions of the vertices, then draw the graph. 

adjacency_table = readtable("Adjacency_matrix.csv");
A = table2array(adjacency_table);
N = length(A);

%% initialize vertex positions
% choose x and y as Gaussian random variables with mean 0 and standard
% deviation N
r = normrnd(0,sqrt(N),N,2);
x = r(:,1); 
y = r(:,2);

tol = 1e-4;
kmax = 1e5;

[x,y,gnorm] = adam_energy(x,y,A,N,kmax,tol);

plot_graph(x,y,A);


figure(2);
semilogy(1:length(gnorm),gnorm);
xlabel('number of iterations');
ylabel('log of norm of the forces function');
title('Norm of the Force vs. iteration')

function [x,y,gnorm] = adam_energy(x,y,A,N,kmax,tol)
% computes the minimum of the energy function for the adjacency matrix
% using the forces function and the Adam optimizer
eta = 0.001;
beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;

w = [x;y];      % initialize the variables we want to optimize (length 2N)

m = zeros(length(w),1);
v = zeros(length(w),1);
gnorm = zeros(kmax,1);

for k = 1:kmax
    g = -forces(x,y,A);
    m = beta1*m + (1 - beta1).*g;
    v = beta2*v + (1 - beta2).*(g.^2);
    mhat = m ./ (1 - beta1^k);
    vhat = v ./ (1 - beta2^k);
    w = w - (eta * mhat) ./ (sqrt(vhat) + epsilon);
    x = w(1:N);
    y = w((N+1):end);

    gnorm(k) = norm(-g);

    if gnorm(k) < tol
        break
    end

    if mod(k,2500) == 0
        fprintf('iteration = %d, norm of the residual is %d\n',k,gnorm(k))
    end
end

end





