function [beta u]=NG(X,y,lam,beta0,r)
% nonnegative garrote
% r=1, there exists intercept  min|y-Xb|^2+2*lam*(u_2+...+u_p)
% other r: no intercept  min|y-Xb|^2+2*lam*(u_1+...+u_p)

[n p]=size(X); 
options=optimset('Display','none','LargeScale','off');
if r==1
    Beta=diag([1;beta0]); 
    H=Beta*X'*X*Beta; H=(H+H')/2; c=[0; lam*ones(p-1,1)]-Beta*X'*y;
    u=quadprog(H,c,[],[],[],[],[-inf; zeros(p-1,1)],inf*ones(p,1),[],options);
    beta=Beta*u;
else
    Beta=diag(beta0); 
    H=Beta*X'*X*Beta; H=(H+H')/2; c=lam*ones(p,1)-Beta*X'*y;
    u=quadprog(H,c,[],[],[],[],zeros(p,1),inf*ones(p,1),[],options);
    beta=Beta*u;
end