function [beta K]=AdaLassoAOEM(X, y, lam, b0, betaE0, R)
% Accelerating OEM algorithm for computing adaptive lasso
% model is y=Xb+e, X,y should not be centered, so in general regression set-up, the first column of X is (1,...,1)'
% R=1, has intercept, adaptive lasso estimator of b is min |y-Xb|^2+2*lam(|b_2|/|b01|+...+|b_{p+1}|/|b0p|),
% other R, no intercept, adaptive lasso estimator of b is min |y-Xb|^2+2*lam(|b_1|/|b01|+...+|b_p|/|b0p|),
% b0 is the initial estimate for adaptive
% betaE0 is initial point for iteration
% beta is the estimator, K is the times of iteration

[A c]=Northar(X);  % add some rows to make X orthogonal
[n p]=size(X); u0=X'*y;
lamp=[0; lam*ones(p-1,1)]; b0p=[1; b0]; % only for the case with intercept
k=1; u=u0+A*betaE0;
if R==1
    Mb=sign(u).*max((abs(u)-lamp./abs(b0p))/c,0);
else
    Mb=sign(u).*max((abs(u)-lam./abs(b0))/c,0);
end
if rd(betaE0,Mb)<1e-4
    betaE(:,k)=Mb;
else
    u=u0+A*Mb;
    if R==1
        MMb=sign(u).*max((abs(u)-lamp./abs(b0p))/c,0);
    else
        MMb=sign(u).*max((abs(u)-lam./abs(b0))/c,0);
    end
    r=Mb-betaE0; v=MMb-Mb-r; gamma=-sqrt((r'*r)/(v'*v));
    betaE(:,k)=betaE0-2*gamma*r+gamma^2*v;
end
sr=rd(betaE0,betaE(:,k));    % stop rule
while sr>=1e-6
    k=k+1; u=u0+A*betaE(:,k-1);
    if R==1
        Mb=sign(u).*max((abs(u)-lamp./abs(b0p))/c,0);
    else
        Mb=sign(u).*max((abs(u)-lam./abs(b0))/c,0);
    end
    if rd(betaE(:,k-1),Mb)<1e-4
        betaE(:,k)=Mb;
    else
        u=u0+A*Mb;
        if R==1
            MMb=sign(u).*max((abs(u)-lamp./abs(b0p))/c,0);
        else
            MMb=sign(u).*max((abs(u)-lam./abs(b0))/c,0);
        end
        r=Mb-betaE(:,k-1); v=MMb-Mb-r; gamma=-sqrt((r'*r)/(v'*v));
        betaE(:,k)=betaE(:,k-1)-2*gamma*r+gamma^2*v;
    end
    sr=rd(betaE(:,k-1), betaE(:,k));
end
beta=betaE(:,k); K=k;