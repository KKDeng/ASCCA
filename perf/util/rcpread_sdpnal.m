function [blk, At, C, b] = rcpread(fname, K, n)
    W0 = importdata(fname);
    
    if isstruct(W0)
        W0 = W0.data; 
    end

    if nargin < 3
        n = size(W0, 1);
    end
    W0 = W0(1:n, :);
    C{1} = sparse(-(W0*W0'));

    blk{1,1} = 's';
    blk{1,2} = n;
    
    Acell = cell(1, n + 1);

    % tr(X) = K
    Acell{1} = speye(n);
    
    % Xe = e
    for i = 1:n
        Ai = zeros(n, n);
        Ai(i, :) = 1;
        Ai = sparse(Ai + Ai') / 2;
        Acell{i+1} = Ai;
    end
    
    At = svec_sdpnal(blk, Acell, 1);
    
    b = [K; ones(n, 1)];
end
