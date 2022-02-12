function [pinf,dinf,relgap,etaC1,etaC2,etaK1,etaK2] = rcp_compute_crit(X,S,Z,y,blk,At,C,b,L,U)


[V,D] = eig(S{1});
D(D<0) = 0;
S{1} = V*D*V';  %% res



pinf = norm(AXfun_sdpnal(blk,At,X) - b)/(1 + norm(b));
dinf = ops(ops(ops(ops(Atyfun_sdpnal(blk,At,y),'+',S),'+',Z), ...
    '-',C),'norm')/(1 + ops(C,'norm'));
pobj = full(ops(ops(C,'.*',X),'sum'));
dobj = b'*y;
relgap = abs(pobj - dobj)/(1+abs(pobj)+abs(dobj));
etaC1 = abs(full(ops(ops(S,'.*',X),'sum')))/(1+ops(X,'norm')+ops(S,'norm'));
etaC2 = abs(full(ops(ops(Z,'.*',X),'sum')))/(1+ops(X,'norm')+ops(Z,'norm'));

X_eig = eig(X{1});
X_eig_neg = X_eig(X_eig < 0);
S_eig = eig(S{1});
S_eig_neg = S_eig(S_eig < 0);


etaK1 = norm(X_eig_neg) / (1 + norm(X_eig));
%etaK2 = norm(S_eig_neg) / (1 + norm(C{1},'fro'));%  re-scale_res
etaK2 = norm(S_eig_neg) / (1 + norm(S_eig));     % re-scale_res_1

% rescale_res_2
% etaK2 = norm(S_eig_neg) / (1 + norm(C{1},'fro'));
% temp = etaK2; etaK2 = dinf;dinf = temp;  
end
