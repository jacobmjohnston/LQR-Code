%% Spring with Viscous Damper

x0 = [5; 0];
N = 100;
x = zeros(2,N+1);
u = zeros(1,N+1);
x(:,1) = x0;
for i = 2:(N+1)
    x(:,i) = fS(x(:,i-1),u(:,i-1));
end
xf = [0;0];

[x,u,runs,error] = ddpIter(@fS,@df,@dl,@dlf,N,x,u,xf,1e-4,100,.25);

ti = 1:N+1;
plot(ti,x(1,:))

%%
J = 0;
for i = 1:N
    J = J + x(1,i)^2 + u(i)^2;
end
J = J + 10*x(1,N+1)^2

%%
% Spring with Viscous Drag System Function and Derivatives
% 0 = x2' + bx2^2 + kx1 - u
% x1' = x1 + 1/t*(x2)
% x2' = x2 + 1/t*(F - bx2*|x2| - kx1 + u) + *w*/t

function[xn] = fS(x,u)
    k = .5;
    b = 0.5;
    dt = .1;
    
    xn(1) = x(1) + dt*x(2);
    xn(2) = x(2) + dt*(-b*x(2)*abs(x(2))-k*x(1) + u);
    xn = xn';
end

function[f] = df(x,u,partial)
    k = .5;
    b = 0.5;
    dt = .1;
    
    if partial == "x"
        f = [1 , dt; -k*dt , 1-dt*2*b*abs(x(2))];
    elseif partial == "u"
        f = [0; dt];
    elseif partial == "xx"
        f(:,:,1) = [0, 0; 0, 0];
        f(:,:,2) = [0, 0; 0, -dt*2*b*x(2)/abs(x(2))];
    elseif partial == "uu"
        f = [0;0];
    elseif partial == "ux"
        f = [0,0;0,0];  
    end
end

function[l] = dl(x,u,partial)
    if partial == "x"
        l = [2*x(1); 0];
    elseif partial == "u"
        l = 2*u;
    elseif partial == "xx"
        l = [2, 0; 0, 0];
    elseif partial == "uu"
        l = 2;
    elseif partial == "ux"
        l = [0,0];  
    end
end

function[lf] = dlf(x,partial)
    if partial == "x"
        lf = [20*x(1); 0];
    elseif partial == "xx"
        lf = [20, 0; 0, 0];
    end
end