function rgrad = egrad2rgrad(X, egrad,invM)
        
%% generalized stiefel manifold
        % First, scale egrad according the to the scaled metric in the
        % Euclidean space. Ideally, B should be preprocessed to ease
        % solving linear systems, e.g., via Cholesky factorization.
        egrad_scaled = invM*egrad;
        
        % Second, project onto the tangent space.
        % rgrad = egrad_scaled - X*symm((B*X)'*egrad_scaled);
        %
        % Verify that symm(BX'*egrad_scaled) = symm(X'*egrad).
        symm = @(D) (D + D')/2;
        rgrad = egrad_scaled - X*symm(X'*egrad);
        
    end