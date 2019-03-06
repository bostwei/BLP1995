%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 1 for Econ 762
% Shiyan Wei
% 03/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This main script is try to do gradian search for the parameter space
clear all
clc
warning('off')
%% Data Use and Generation 
load BLP_1999.mat

global x 


N = size(car_id,1); % number of obs.
cons = ones(N,1); % constant
x = [cons, hpwt, air, mpd, space,trend]; % grouping the x variable

[model_group,ID] = findgroups(model_id); % initiate assign group number to each group
G = size(ID,1);% number of groups

%% Step 1: Given delta, compute the mkt share 
% generate EVT1 for each group, eij
e = evrnd(0,1,N,1);

% Here we use the price of other goods as the instrument( Hausmann Instrument) 
price_til = zeros(N,1);
for j = 1: N
    price_til(j,:) = sum(price,1) - price(j,:);
end
z = [x,price_til];

%% --------- step 4 search for the optimal theta ------------------
% Initialized guess of sigma
s = linspace(0.1,5,10)';

sigma = cartesian(s,s,s,s,s,s); 

SN = size(sigma,1);

% segment of partition
seg = 100;
% length of the partition
l = ceil(SN/seg);

for index = 1:seg
    % begin index of seg
    begin = (index-1) * l+1;
    endd =  (index) * l;
    % partition sigma into different portion
    genvarname('sigma_seg');
    eval(['sigma_seg','= sigma(begin:endd,:)']);
    filename=strcat('input/sigma',int2str(index));
    save( filename, 'sigma_seg');
end

save('input/base.mat','z','x', 'price','share');
%% create a reading calabrating file

d = size(sigma1,1);

for i =1: d     
    [theta1,theta2,theta3,theta4,theta5,theta6,theta7,fval_i] = SearchTheta(sigma(d));
end
save('result.mat');
return