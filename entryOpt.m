%inital point
x0 = [6478000 7450 -.5/180*pi 5*pi/180];

dt = 3.2;
N = 500;

%initial x and u trajectory
x = zeros(4,N+1);
u = zeros(1,N+1);
x(:,1) = x0;
for i = 2:(N+1)
    x(:,i) = x(:,i-1) + dt*entry(x(:,i-1),u(i-1));
end

%sets endpoint and cost matrices
xf = [6403000; 0.1; -10*pi/180; 0];
Qf = [.001 0 0 0;
      0 10 0 0;
      0 0 23 0;
      0 0 0 5];
Q = [.001 0 0 0;
     0 10 0 0;
     0 0 17 0;
     0 0 0 5];
R = 1000;

%Call to iLQR algorithm
[x,u,runs,error] = ilqrIter(@entry,@dA,@dB,Q,R,Qf,N,x,u,xf,2,1e-5,2,dt);


t = dt*(0:N);
%% Nonlinear function for trajectory propagation
function dy = entry(y,u)
    R0 = 6378000;
    g0 = 9.81;
    %w 7.2921159e-5;

    Aref = 391.22;
    density = 1.23*exp(-(y(1)-R0)/7990);
    m = 104305;
    
    if y(2) > 4570 
        a = 40;
    else
        a = 40 - 0.20705*(y(2)-4570)^2/340^2;
    end
    
    Cl = -0.041065 + 0.016292*a + 0.0002602*a^2;
    Cd = 0.080505 - 0.03026*Cl + 0.86495*Cl^2;
    
    L = 1/2*density*Aref*(y(2))^2*Cl/m;
    D = 1/2*density*Aref*(y(2))^2*Cd/m;
    
    dy(1) = y(2)*sin(y(3));
    dy(2) = -D - R0^2*g0*sin(y(3))/(y(1))^2;
    dy(3) = L*cos(y(4))/y(2) + (y(2)/y(1) - (R0^2*g0)/(y(2)*(y(1))^2))*cos(y(3));
    dy(4) = u;
    
    dy = dy';
end

% Linearized matrix for r,v,gamma,bank
function[A] = dA(x,u)
    R0 = 6378000;
    g0 = 9.81;
    Aref = 391.22;
    density = 1.23*exp(-(x(1)-R0)/7990);
    m = 104305;
    if x(2) > 4570 
        a = 40;
        Cl = -0.041065 + 0.016292*a + 0.0002602*a^2;
        Cd = 0.080505 - 0.03026*Cl + 0.86495*Cl^2;
        dcddv = 0;
        dcldv = 0;
    else
        a = 40 - 0.20705*(x(2)-4570)^2/340^2;
        Cl = -0.041065 + 0.016292*a + 0.0002602*a^2;
        Cd = 0.080505 - 0.03026*Cl + 0.86495*Cl^2;
        dadv = (0.4141)*(4570-x(2))/340^2;
        dcldv = (0.016292 + 0.0005204*a)*dadv;
        dcddv = (-0.03026 + 1.7299*Cl)*dcldv;
    end    
    L = 1/2*density*Aref*(x(2))^2*Cl/m;
    D = 1/2*density*Aref*(x(2))^2*Cd/m;
  
    A = [0, sin(x(3)), x(2)*cos(x(3)), 0; %r
        1/7990*D+2*R0^2*g0*sin(x(3))/(x(1)^3), -2*D/x(2)-D/Cd*dcddv, ... %v
        -R0^2*g0*cos(x(3))/(x(1)^2), 0;
        -1/7990*L*cos(x(4))/x(2) + (-x(2)/(x(1))^2 +2*R0^2*g0/x(2)/(x(1))^3)*cos(x(3)), ... %gamma
        -L*cos(x(4))/(x(2)^2) + (2*L/x(2) + dcldv*L/Cl)*cos(x(4))/x(2) + (1/x(1) + R0^2*g0/(x(1)*x(2))^2)*cos(x(3)), ...
        -(x(2)/x(1)-(R0^2*g0)/(x(2)*(x(1))^2))*sin(x(3)), -L*sin(x(4))/x(2);
        0, 0, 0, 0];       
end

%linearized input matrix
function[B] = dB(x,u)
    B = [0; 0; 0; 1];
end
