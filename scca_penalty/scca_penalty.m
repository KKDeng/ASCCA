function [A,B] = scca_penalty(sigmaX,sigmaY,sigmaXY,X,Y,u0,v0,lambda1,lambda2)


r = size(u0,2);  [n,p] = size(X); q = size(Y,2);

Omega = eye(n);  A = []; B = []; R = [];


for i = 1:r
    [u,v] = scca_ipls(X,Y,sigmaX,sigmaY,Omega,u0(:,i),v0(:,i),lambda1,lambda2);
    
    A = [A,u]; B = [B,v]; R = [R,u'*sigmaXY*v];
    
    Omega = eye(n) - (X*A*diag(R)*B'*Y')/n; 
    
end

