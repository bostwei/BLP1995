function [theta_hat_i,fval_i] = SearchTheta(sigma,x,share)

global delta
% sigma = [s1,s2,s3,s4,s5,s6]; 

N  = size(x,1);
    
    sigma_i = sigma';
%     sigma_i = ones(6,1);
    
    % initial guess of delta
    delta = 2 * ones(N,1);
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

    % fprintf('This is %d itertion. The difference is %.3f \n', iter, diff);
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
    
    % using the 1 step GMM estimator
    [theta_hat_i,fval_i]=fminsearch(@gmm,theta);
    
end

