%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 1 for Econ 762
% Shiyan Wei
% 03/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This main script is to use MPEC to speed up the computation
clear all
clc
warning('off')
%%% Data Use and Generation 
load BLP_1999.mat

global x z price share G N Mkt price lambda



N = size(car_id,1); % number of obs.
cons = ones(N,1); % constant

x = [cons, hpwt, air, mpd, space,trend]; % grouping the x variable

[Mkt,Mkt_ID] = findgroups(year1); % initiate assign group number to each group
Mkt_ID = Mkt_ID-70;

G = size(Mkt_ID,1);% number of groups

% Here we use the price from other mkt as the instrument( Hausmann Instrument) 

g_price = splitapply(@mean,price,Mkt); % the mean price in a mkt
gg_price = zeros(N,1); % average group price

for g=1:G
    gg_price(Mkt==g) = g_price(g);
end

gg_num = zeros(N,1); % member number in a group

for g=1:G
    gg_num(Mkt==g) = size(Mkt(Mkt==g),1);
end


price_til = (sum(price,1) - gg_price .* gg_num)./(N- gg_num); % mean price in other mkt 

z = [x,price_til];

%% Step 1: Given delta, compute the mkt share 

% Initialized guess of sigma
sigma = ones(size(x,2),1);

% initial guess of delta
delta = ones(G,1); % given initial guess of delta for each group

% initial guess of beta
beta = ones(size(x,2),1);
alpha = 1;

lambda =2;
% Initialized paramter
theta = [alpha;beta;delta;sigma];
% Initialized guess
alpha_ub = 60;
alpha_lb = -60;
beta_ub = [-5;9;4;2;6;1];
beta_lb = [-9;-9;-4;-2;-6;-1];
delta_ub = 80*ones(G,1);
delta_lb = zeros(G,1);
sigma_ub = 5 * ones(size(x,2),1);
sigma_lb = 0 * ones(size(x,2),1);
ub = [alpha_ub;beta_ub;delta_ub;sigma_ub ];
lb = [alpha_lb;beta_lb;delta_lb;sigma_lb ];
func = @MPECgmm;

opts = optimoptions('fmincon','UseParallel',true );

% [theta_hat,Qval]=particleswarm(func,1);%,[],[],[],[],lb,ub);%,[],opts);

lambda_loop =linspace(2,10,10);
for l = 1:size(lambda_loop)  
lambda = lambda_loop(l);
tic
[theta_hat(:,l),Qval(:,l)]=fmincon(func,theta,[],[],[],[],lb,ub);%,[],opts);
toc
end
