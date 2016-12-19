function k=Ka(n,p,h)
% correction factor for scale estimator based on LTS

a=(h-fix((p+1)/2))/n;
q=norminv(a/2+1/2);
M=5000; 
u=linspace(-q,q,M);
U=u.^2.*normpdf(u);
k=sqrt(h/n)/sqrt(mean(U));