function d=runpbindsim

k(1)=10^3;
k(2)=10^2;

km=10^2;

y0(1)=1;
y0(3)=0;

for i=1:9
    y0m=0.1*(2^(i-4));
    for j=1:9
        k(1)=km*(2^(j-4));
        data.k=k;

        tlim=100;
        opts=odeset('RelTol',1e-2,'AbsTol',1e-4,'MaxStep',1);
        
        for m=1:9;
            y0(2)=y0m*(2^(m-4));
            [T,Y] = ODE15s(@pbinds,0:1:tlim,y0,opts,data);
            dc(m)=Y(length(T),3);
        end
        %[p,r]=fit((2.^(1:9))',dc','a*x+b');
        d(i,j)=corr((2.^(1:9))',dc');
    end
end




end


function yp=pbinds(t,y,dat)
    k=dat.k;
    yp(1)=y(3)*k(1)-y(1)*y(2)*k(2);
    yp(2)=y(3)*k(1)-y(1)*y(2)*k(2);
    yp(3)=y(1)*y(2)*k(2)-y(3)*k(1);
    yp=yp';
end
