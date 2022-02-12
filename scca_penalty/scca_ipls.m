function [u,v] = scca_ipls(X,Y,sigmaX,sigmaY,Omega,u0,v0,lambda1,lambda2)
n = size(X,1);


opts.maxIts = 1000; opts.tol = 1e-4; opts.printEvery = 0;
maxiter = 100;

for i = 1:maxiter
    
    hatX = Omega'*X*u0;
    [ v, ~, ~ ] = solver_L1RLS( Y, hatX , lambda2, v0, opts );
    if norm(v) ~=0
    v = v / sqrt(v'*sigmaY*v);
    end
    hatY = Omega'*Y*v;
    [ u, ~, ~ ] = solver_L1RLS( X, hatY , lambda1, u0 , opts);
    if norm(u) ~= 0
    u = u / sqrt(u'*sigmaX*u);
    end
    
    if max(norm(u-u0),norm(v-v0)) < 1e-4
        break
    end
    
    u0 = u; v0 = v;
end



end

