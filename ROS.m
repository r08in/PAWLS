function [beta w]=ROS(X,y,beta0,mu)
% ROS with tuning mu for model y=Z*theta+e
% initial estimator of delta is the residual of theta0

[n p]=size(X);
delta0=y-X*beta0; Del=diag(delta0);
w=mu./(delta0.^2+mu); Q=X'*diag(w)*X;
beta=Q\X'*diag(w)*y;