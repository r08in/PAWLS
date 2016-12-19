function [beta sig]=LTS(X,y,h)

% fast-LTS (least trimmed squares) agorithm when n<=600
% X: the regression matrix (the first column is ones when intercept exists)
% y: response
% h: <=[n/2]+[(p+1)/2]


% Ref: Rousseeuw, P. J. and Van Driessen, K. (2006). computing LTS regression for large data sets

[n p]=size(X); R=500;
for r=1:R
    H0=randsample(1:n,p);
    X0=X(H0,:); eigen=eig(X0); reig=eigen(p)/eigen(1);
    while reig>1e14
        k0=randsample(setdiff(1:n,H0),1);
        H0(1)=k0; X0=X(H0,:); eigen=eig(X0); reig=eigen(p)/eigen(1);
    end
    b0=X0\y(H0);
    rsd=y-X*b0; [sar Ind]=sort(abs(rsd)); H1=Ind(1:h);
    X0=X(H1,:); y0=y(H1); b0=(X0'*X0)\(X0'*y0); 
    rsd=y-X*b0; [sar Ind]=sort(abs(rsd)); H2(r,:)=Ind(1:h);
    X0=X(H2(r,:),:); y0=y(H2(r,:)); b0=(X0'*X0)\(X0'*y0); 
    rsd=y-X*b0; ssr=sort(rsd.^2); Q0(r)=sum(ssr(1:h));
end
[sQ J]=sort(Q0); J=J(1:10); mQ=inf;
for r=1:10
    k=1; H(k,:)=H2(J(r),:); 
    X0=X(H(k,:),:); y0=y(H(k,:)); b0=(X0'*X0)\(X0'*y0); 
    rsd=y-X*b0; [ssr Ind]=sort(rsd.^2); H(k+1,:)=Ind(1:h); 
    while sort(H(k,:))~=sort(H(k+1,:))
        k=k+1;
        X0=X(H(k,:),:); y0=y(H(k,:)); b0=(X0'*X0)\(X0'*y0); 
        rsd=y-X*b0; [ssr Ind]=sort(rsd.^2); H(k+1,:)=Ind(1:h); 
    end
    Q=sum(ssr(1:h));
    if Q<mQ
        beta=b0;
    end
end

rsd=(y-X*beta).^2; srsd=sort(rsd);
sig=Ka(n,p,h)*sqrt(mean(srsd(1:h)));