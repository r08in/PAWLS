function [X,y]=GenerateData()
beta=[3 2 1.5 0 0 0 0 0]'; tag=[1; 2; 3];
%beta=[3 2 1.5 1 1 1 1 1]';
n=50
p=length(beta); 
no=10
sig=2
rho=0.5
for i=1:p
    for j=1:p
        Sig(i,j)=rho^(abs(i-j));
    end
end
cS=chol(Sig);
    for i=1:n
        X(i,:)=cS'*normrnd(0,1,[p 1]);
    end
    y=X*beta+normrnd(0,sig,[n 1]);
    if no>0
        for i=1:no
            X(n-i+1,1)=X(n-i+1,1)+10;
            y(n-i+1)=y(n-i+1)+(-1)^randsample(2,1)*unifrnd(20,30);
            %X(n-i+1,:)=X(n-i+1,:)+3*rand(1,p);
        end
    end