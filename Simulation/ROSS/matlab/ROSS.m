function betaROSS=ROSS(X, y, w, betaROS, betaNGROSS, r)  % ROSS regression

% r=1, has intercept; other r, no intercept

[n p]=size(X);
sigEs=(y-X*betaROS)'*diag(w)*(y-X*betaROS)/sum(w);
lam=log(n)*sigEs/2;
Xn=diag(sqrt(w))*X; yn=diag(sqrt(w))*y;
if r==1
    [betaROSS K]=AdaLassoAOEM(Xn, yn, lam, betaROS(2:p), betaNGROSS, r); 
else
    betaROSS=AdaLassoAOEM(Xn, yn, lam, betaROS, betaNGROSS, r);
end