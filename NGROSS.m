function [betaE u]=NGROSS(X, y, w, betaROS, r)  % NG-ROSS regression

% r=1, has intercept; other r, no intercept

[n p]=size(X);
sigEs=(y-X*betaROS)'*diag(w)*(y-X*betaROS)/sum(w);
lam=log(n)*sigEs/2;
Xn=diag(sqrt(w))*X; yn=diag(sqrt(w))*y;
if r==1
    Xn=Xn;
    [betaE u]=NG(Xn, yn, lam, betaROS(2:p), r); 
else
    [betaE u]=NG(Xn, yn, lam, betaROS, r);
end