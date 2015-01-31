function z=MC(X, y, beta0, mu)

% MC value of a ROS with initial beta0 (LTS with robustest h)

[n p]=size(X); b=ROS(X,y,beta0,mu); cvv=(y-X*b).^2; scv=sort(cvv); 
h0=fix(n/2)+fix((p+1)/2);
for h=h0:n    
    oz(h-h0+1)=sum(scv(1:h))/nchoosek(h+8,9);
end
z=min(oz);