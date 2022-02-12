function dobj = dualobj(b, y, Z, L, U)
    dobj = b' * y;
    for k = 1:length(Z)
        n = length(Z{k});
        mZ = -Z{k};
        G = get_content(L, k, n, -inf(n, n));
        G(mZ == 0) = 0;
        Uk = get_content(U, k, n, inf(n, n));
        G(mZ > 0) = Uk(mZ > 0);
        dobj = dobj - sum(sum(mZ.*G));
    end
end

function Xk = get_content(X, k, n, default)
    if iscell(X)
        Xk = X{k};
    elseif length(X) == n
        Xk = X
    elseif length(X) == 1
        Xk = X * ones(n, n);
    elseif isempty(X)
        if nargin < 4
            Xk = [];
        else
            Xk = default;
        end
    else
        error('fail to get content');
    end
end
