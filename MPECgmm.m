function [ob] = MPECgmm(alpha,beta,delta, sigma, lambda)
% This function calculate the crterion function given alpha beta, delta,
% and sigma
global x share G N Mkt price z

% given alpha beta delta sigma
% the moment condition is 

g_delta = zeros(N,1);

for g = 1:G
    g_delta(Mkt==g) = delta(g);
end

xi = g_delta - (x * beta - alpha * price);
g =  z' * xi ./G; %?
W = inv(z'*z/G);

% the objected function will be
ob = g'* W * g - lambda * sum((findMktShare(g_delta,sigma,x) - share));


end