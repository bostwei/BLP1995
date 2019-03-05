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

global x z price W delta


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
s_j = findMktShare(delta,sigma,x);

% construct contract mapping
delta1 = delta + log(share)- log(s_j);

% calculate the difference
diff = max(abs(delta1 - delta)./abs(delta));

fprintf('This is %d itertion. The difference is %.3f \n', iter, diff);
% update guess of delat
delta = delta1;

iter = iter +1;
end

%-------------- step 3 find alpha and beta using GMM estimator -----------------
% Here we use the price of other goods as the instrument( Hausmann Instrument) 
price_til = zeros(N,1);
for j = 1: N
    price_til(j,:) = sum(price,1) - price(j,:);
end
z = [x,price_til];

% Initial guess of theta
beta = ones(size(x,2),1);
alpha = 1;

theta = [beta;alpha];

%-------- using 2step GMM estimator ---------------
% The weighted matrix 
W = eye(size(theta,1));
% using the 1 step GMM estimator
[theta_hat,fval]=fminsearch(@gmm,theta);

% beta_hat = theta_hat(1:size(x,2),:);
% alpha_hat = theta_hat(size(x,2)+1,:);
% 
% xi_hat = delta - (x * beta_hat - price * alpha_hat);  
% % the moment condition
% g_hat =  z'* xi_hat;
% Omega = g_hat * g_hat'./N;
% W = inv(Omega);
% [theta_hat,fval]=fminsearch(@gmm,theta_hat);

%% --------- step 4 search for the optimal theta ------------------
% Initialized guess of sigma
s = linspace(0.1,10,10)';

sigma = cartesian(s,s,s,s,s,s); 
SN = size(sigma,1);


for i = 1:SN
    
    sigma_i = sigma(i,:)';

    % initial guess of delta
    delta = 40 * ones(N,1);
    % find delta(sigma)
    % initialized loop parameter
    Tol = 0.1;
    MaxIter = 1000;
    diff = 100;
    iter = 1;

    while diff>Tol && iter < MaxIter
    %-------------- step 1 find Mkt Share --------------
    s_j = findMktShare(delta,sigma_i,x);

    % construct contract mapping
    delta1 = delta + log(share)- log(s_j);

    % calculate the difference
    diff = max(abs(delta1 - delta)./abs(delta));

%     fprintf('This is %d itertion. The difference is %.3f \n', iter, diff);
    % update guess of delat
    delta = delta1;

    iter = iter +1;
    end


    % Initial guess of theta
    beta = ones(size(x,2),1);
    alpha = 1;

    theta = [beta;alpha];

    %-------- using 2step GMM estimator ---------------
    % The weighted matrix 
    W = eye(size(theta,1));
    % using the 1 step GMM estimator
    [theta_hat_i,fval_i]=fminsearch(@gmm,theta);
    theta_hat(:,i) = theta_hat_i;
    fval(:,i) = fval_i;
    
    fprintf('This is %d itertion. \n', i);
end