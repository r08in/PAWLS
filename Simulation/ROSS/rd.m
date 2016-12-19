function z=rd(a0,a1)   % relative difference between a0 and a1
p=length(a0); z=0;
for j=1:p
    if abs(a0(j))<=1e-10   % a0(j) appro 0
        if  abs(a1(j))>1e-10  
            z=inf;
            break
        end
    else
        z=max(z,abs(a1(j)-a0(j))/abs(a0(j)));
    end
end