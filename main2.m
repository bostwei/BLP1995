%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 1 for Econ 762
% Shiyan Wei
% 03/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This main script is doing the IV Logit Regression for replicating the BLP 1995
clear
clc

%% Data Use and Generation 
load BLP_1999.mat

global x sigma N


N = size(car_id,1); % number of obs.
cons = ones(N,1); % constant
x = [cons, hpwt, air, mpd, space,trend]; % grouping the x variable

[model_group,ID] = findgroups(model_id); % initiate assign group number to each group
G = size(ID,1);% number of groups

%% Step 1: Given delta, compute the mkt share 
% generate EVT1 for each group, eij
e = evrnd(0,1,N,1);



% Initialized guess of sigma
sigma = ones(size(x,2),1);

% initial guess of delta
delta = 0.5 * ones(N,1);
%-------------- step 2 find delta(sigma) -----------------
% initialized loop parameter
Tol = 0.1;
MaxIter = 1000;
diff = 100;
iter = 1;

while diff>Tol && iter < MaxIter
%-------------- step 1 find Mkt Share --------------
s_j = findMktShare(delta);

% construct contract mapping
delta1 = delta + log(share)- log(s_j);

% calculate the difference
diff = max(abs(delta1 - delta)./abs(delta));

fprintf('This is %d itertion. The difference is %.3f \n', iter, diff);
% update guess of delat
delta = delta1;

iter = iter +1;
end

%-------------- step 3 find alpha and beta given Moment Condition -----------------


