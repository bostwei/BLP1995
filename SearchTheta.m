function [theta_hat_i,fval_i] = SearchTheta(s1,s2,s3,s4,s5,s6)
global x share W

sigma = [s1,s2,s3,s4,s5,s6]; 

N  = size(x,1);
    
    sigma_i = sigma';

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
%     theta1 = theta_hat_i(1);
%     theta2 = theta_hat_i(2);
%     theta3 = theta_hat_i(3);
%     theta4 = theta_hat_i(4);
%     theta5 = theta_hat_i(5);
%     theta6 = theta_hat_i(6);
%     theta7 = theta_hat_i(7);
    
end

