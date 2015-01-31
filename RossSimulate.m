function [betaROSS]=RossSimulate(X,y)
%intial
n=length(y);
d=size(X);
p=d(2);

betaLS=(X'*X)\(X'*y);
%mse_LS(r)=(betaLS-beta)'*(betaLS-beta);
h=fix(n/2)+fix((p+1)/2)
[betaLTS sigLTS]=LTS(X,y,h);
%mse_LTS(r)=(betaLTS-beta)'*(betaLTS-beta);
[betaROS w]=ROS_mc(X, y, betaLTS);
%mse_ROS(r)=(betaROS-beta)'*(betaROS-beta);    
betaNGROSS=NGROSS(X, y, w, betaROS, 0);
betaROSS=ROSS(X, y, w, betaROS,betaNGROSS,0);
