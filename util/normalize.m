function [X,Y,XtY,M1,M2] = normalize(X,Y)
[m,p] = size(X);
[~,q] = size(Y);

X = X-  repmat(mean(X,1),m,1);
Y = Y - repmat(mean(Y,1),m,1);
%m = 2;
X = X/sqrt(m-1); Y = Y/sqrt(m-1);
XtY = X'*Y; % covariance
XtX = X'*X;
YtY = Y'*Y;
[VX, DX] = eig(XtX);DX = diag(DX);
gamma = 1e-2;
if svds(X,1,'smallest') < 1e-6*max(abs(DX))
    M1 = (1-gamma)*XtX + gamma*eye(p);
else
    M1 = XtX;
end

[VX, DX] = eig(YtY);DX = diag(DX);
if svds(Y,1,'smallest') < 1e-6*max(abs(DX))
    M2 = (1-gamma)*YtY + gamma*eye(q);
else
    M2 = YtY;
end


end