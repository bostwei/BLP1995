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

global x_j price_j xi_j


N = size(car_id,1); % number of obs.
cons = ones(N,1); % constant
x = [cons, hpwt, air, mpd, space,trend]; % grouping the x variable

[model_group,ID] = findgroups(model_id); % initiate assign group number to each group
G = size(ID,1);% number of groups

% Initial guess of the parameter
beta = ones(size(x,2),1); 
alpha = 1;

xi_j = ones(G,1);
%% Traditional Logit Model for Mkt Share
% Assume the outside goods share be 0

% generate EVT1 for each group
e = evrnd(0,1,G,1);

% calculate the share of each group
share_j = splitapply(@sum,share,model_group);
share_j = share_j/sum(share); % standardized of share_j 

% calculate the x of group by mean 
cons_j = splitapply(@mean,cons,model_group);
hpwt_j = splitapply(@mean,hpwt,model_group);
air_j = splitapply(@mean,air,model_group);
mpd_j = splitapply(@mean,mpd,model_group);
mpg_j = splitapply(@mean,mpg,model_group);
space_j = splitapply(@mean,space,model_group);
trend_j = splitapply(@mean,trend,model_group);

% x_j = [cons_j, hpwt_j, air_j, mpd_j, mpg_j, space_j,trend_j];
x_j = [cons_j, hpwt_j, air_j, mpd_j, space_j,trend_j];
% calculate the price of group by mean
price_j = splitapply(@mean,price,model_group);
%------------------- Step 1 --------------------------------------------
% given initial guess of delta, compute the mkt share

% calcuate initial guess of the delta 
delta0_j = x_j * beta - alpha * price_j  +  xi_j;


% calculate the implied mkt share
 s_hat = exp(delta0_j)/(1 + sum(exp(delta0_j)));  

% ----------------- Step 2 ----------------------------------------------
% find the vector delta that equate the observed mkt share and predicted
% mkt share by contract mapping

% initialized loop parameter
Tol = 0.1;
MaxIter = 1000;
diff = 100;
iter = 1;

while diff>Tol && iter < MaxIter
% given inital guess of delta calulate implicited mkt share    
s_hat = exp(delta0_j)/(1 + sum(exp(delta0_j)));  

% construct contract mapping
delta1_j = delta0_j + log(share_j)- log(s_hat);

% calculate the difference
diff = max(abs(delta1_j - delta0_j));

fprintf('This is %d itertion. The difference is %.3f \n', iter, diff);
% update guess of delat
delta0_j = delta1_j;

iter = iter +1;

end
%----------------- Step 3 computing xi ---------------------------
% Here we use BLP instrument that Zj = sum xi, where i~=j ( May have the
% Over ID problem)

% Here we use the price of other goods as the instrument( Hausmann Instrument) 
for j = 1: G
    z(j,:) = sum(price_j,1) - price_j(j,:);
end
z = [ones(G,1),z];
% then the moment condition for the estimation is E[xi*z] = 0
% since this is linear model, we do this by 2sls
[par_z,~,~,Cov_z] = mvregress(z,price_j); % 1st stage regression

p_j_hat = z * par_z; % predict the p_j_hat

X = [x_j, -p_j_hat];

[par,~,E,CovB] = mvregress(X, delta1_j); 
% - par is paramter of beta and alpha
% - E is residual of the linear regression, xi here
% - CovB is the varaiabce-covariance matrix of the coefficient

% calcualating xi
xi_j = E;








