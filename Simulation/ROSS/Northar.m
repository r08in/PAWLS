function [A c K]=Northar(X)  % new version (power method) for orthogonalizing a design matrix X by adding some rows. 
% A=Xa'*Xa, where Xa is the added matrix
% c is the largest eigenvalue of X'*X
[n p]=size(X); B=X'*X;
k=1; z0=rand(p,1);
z(:,1)=B*z0/sqrt(z0'*z0); sr=(z(:,1)-z0)'*(z(:,1)-z0);
while sr>=1e-10
    k=k+1;
    z(:,k)=B*z(:,k-1)/sqrt(z(:,k-1)'*z(:,k-1)); sr=(z(:,k)-z(:,k-1))'*(z(:,k)-z(:,k-1));
end
c=sqrt(z(:,k)'*z(:,k));
A=c*eye(p)-B;
K=k;