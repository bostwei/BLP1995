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

global x z price share G N Mkt price



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
delta = Mkt_ID; % given initial guess of delta for each group

% initial guess of beta
beta = ones(size(x,2),1);
alpha = 1;

lambda =1;

[Q] = MPECgmm(alpha,beta,delta, sigma,lambda);





