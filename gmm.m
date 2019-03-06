function [Q] = gmm(theta)
% This is the GMM function of BLP 2015 
% - theta = [beta,alpha];
% - Q is the criterion function
% - z is the instrument

global z x delta price 
% parameter to estimate 
beta = theta(1:size(x,2),:);
alpha = theta(size(x,2)+1,:);

xi = delta - (x * beta - price * alpha);  
% the moment condition
m =  z' * xi;
Omega = (z' * xi) * ( z' * xi)' ;
% generate the size of the sample 
n = size(m,1);
W = inv(Omega);
% W = eye(size(theta,1));
% the GMM criterion function 
Q = n * m' * W * m;


end
