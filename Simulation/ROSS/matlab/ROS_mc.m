function [betaROS w mu]=ROS_mc(X, y, beta0)  % ROS regression with tuning by minimizing MC

[n p]=size(X); 

Lower0=0.1; Upper0=n/10; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end

Lower0=n/10; Upper0=n/4; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

Lower0=n/4; Upper0=n/2; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

Lower0=n/2; Upper0=n*3/4; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

Lower0=n*3/4; Upper0=n; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

Lower0=n; Upper0=n^(3/2); rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

Lower0=n^(3/2); Upper0=n^2; rt=(sqrt(5)-1)/2; delta=1e-2;
Lower=Lower0; Upper=Upper0; phi=Lower+(1-rt)*(Upper-Lower); psi=Lower+rt*(Upper-Lower); k=1;
while k>=1
    if MC(X, y, beta0, phi)>MC(X, y, beta0, psi)
        if Upper-phi<=delta
            mu2=psi;
            break
        else
            Lower=phi; phi=psi; psi=Lower+rt*(Upper-Lower); k=k+1;
        end
    else
        if psi-Lower<=delta
            mu2=phi;
            break
        else
            Upper=psi; psi=phi; phi=Lower+(1-rt)*(Upper-Lower); k=k+1;
        end
    end
end
if MC(X, y, beta0, mu)>MC(X, y, beta0, mu2)
    mu=mu2;
end

[betaROS w]=ROS(X,y,beta0,mu); 