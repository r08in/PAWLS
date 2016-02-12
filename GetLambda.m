function [lamdas,deltas]=GetLambda(X,y)
% Get range of lambda for ridge penalty

[n p]=size(X); 
if nargin<3, numlam=min([20 n size(X,2)]); end
if nargin<4, cualcv=0; end
if nargin<5, showhist=0; end
if nargin<6, nkeep=5; end

%Normalize and center X and y
[Xnor ynor mux sigx muy]=prepara(X,y);
%Spherical Principal COmponents (no centering)
%privar, Beig= vector of robust "eigenvalues" and matrix of eigenvectors
%Xnor is now =PCA scores= "orthonormalized Xnor "
[privar Beig muspam Xnor]=SPC(Xnor,0); 
[n p]=size(Xnor);  %p is now the "actual" dimension
privar=privar*n; %Makes the robust eigenvalues of the same order as those of classical PCA used for LS
 nlam=min([p numlam]);  
pmax=min([p n/2]);   %edf<=n/2 to keep BDP >=0.25
pp=linspace(1,pmax,nlam);  %"candidate edf's"
lamdas=findlam(privar,pp); %find lambdas corresponding to the edf's
deltas=0.5*(1-pp/n);  %for the M-escale used with Peña-Yohai
