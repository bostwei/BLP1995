function [s_j] = findMktShare(delta,sigma,x)
% This function provide a method that given and vik, we compute a market
% share
% set sample times
I = 1000; % set sample time
N = size(x,1);
% generate preference shock vik ~ N(0,1) for each perference
v = normrnd(0,1,size(x,2),I);
s_ij = zeros(N,I);
for i = 1:I 
   v_i = v(:,i);  
% calculate the utility over different goods and charateristics
enu = exp(delta + sum((x .* (ones(N,1) * v_i') .* (ones(N,1) * sigma')),2) );  % numerator
den = 1+ sum(exp(delta + sum((x .* (ones(N,1) * v_i') .* (ones(N,1) * sigma')),2) ),1); % denorminator
s_ij(:,i) = enu./den;
end

s_j = sum(s_ij,2)/I;

end