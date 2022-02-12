function u = AM_tracelasso_subproblem(a,X,XtX,mu,rho1,rho2,err,u,maxiter,sigma)

%%  min_u a'*u + \mu*\|Xdiag(u)\|_* £¬s.t. \|Xu\| = 1

%% min_u a'*u + \mu*\|P\|_* ,  Xu = y, Xdiag(u) = P, \|y\| = 1

%%  \min_u a'*u + \mu*\|P\|_* + rho1/2\|Xu - y + lambda_1/rho_1\|^2  + rho2/2\|Xdiag(u) - P + Lambda_2/rho2\|^2 

lambda1 = zeros(size(X,1),1); Lambda2 = zeros(size(X));
iter = 1;
while(iter<maxiter)
    temp = X*u+ lambda1/rho1;
    y = temp/norm(temp);
    
    temp = X*diag(u) + Lambda2/rho2;
    [P,~] = prox_nuclear(temp,mu/rho2);
    
    
    A = rho1*XtX + rho2*diag(diag(XtX)) + sigma*eye(size(X,2)) ;  b = rho1*X'*(y - lambda1/rho1) + rho2*diag(X'*(P-Lambda2/rho2)) - a + sigma*u;
    u = A\b;
    

    
    lambda1 = lambda1 + rho1*(X*u - y);
    Lambda2 = Lambda2 + rho2*(X*diag(u) - P);
    iter = iter + 1;
    
    if(norm(X*u-y)<err || norm(X*diag(u)-P)<err)
        break;
    end

    
end

end