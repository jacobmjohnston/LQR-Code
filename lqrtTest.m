%% Linear Deterministic Example
N = 20;
A = zeros(2,2,N);
B = zeros(2,1,N);
for i = 1:N
    A(:,:,i) = [1 1; 0 1];
    B(:,:,i) = [0;1];
end
Q = [1 0; 0 1];
Qf = Q;
R = 10;
K = lqrSolve(A,B,Q,R,Qf,N);
[x,u] = lqrState([1;0],K,A,B,[0 0;0 0],20);
y = [1 0] * x;

t = 1:(N+1);

figure(1);
plot(t,y)
xlabel("Time");
ylabel("y");

figure(2);
plot(t,u)
xlabel("Time");
ylabel("u");
%% Linear Stochastic Example
Q = diag(ones(1,5));
Qf = 10*diag(ones(1,5));
R = diag(ones(1,2));
N = 30;
A = zeros(5,5,N);
B = zeros(5,2,N);
for i = 1:N
    A(:,:,i) = [1 1 0 0 0; 0 1 1 0 0; 0 0 1 0 1; 0 1 0 1 0; 0 0 0 0 1];
    B(:,:,i) = [0 0; 1 0; 0 1; 0 0; 1 1];
end

K = lqrSolve(A,B,Q,R,Qf,N);
W = diag(ones(1,5))/4;
[x,u] = lqrState([1 0 0 0 0], K, A, B, W, N);

figure(3)
t = 1:(N+1);
hold on;
plot(t,x(1,:));
plot(t,x(2,:));
plot(t,x(3,:));
plot(t,x(4,:));
plot(t,x(5,:));


%% LQR Example - Spring with Damper
Q = [1 0; 0 0];
Qf = 10*[1 0; 0 1];
R = 1;
N = 100;
A = zeros(2,2,N);
B = zeros(2,1,N);
t = 10;
k = .5;
b = 0.5;
for i = 1:N
    A(:,:,i) = [1 1/t; -k/t 1-b/t];
    B(:,:,i) = [0; 1/t];
end

K = lqrSolve(A,B,Q,R,Qf,N);
[x,u] = lqrState([5 0], K, A, B, [0 0; 0 1/t], N);

figure(4)
ti = 1:N+1;
hold on
plot(ti, x(1,:), "LineWidth", 1.3)
plot(ti, x(2,:), "LineWidth", 1.3)
plot(ti, u, "LineWidth", 1.3)


%% iLQR Example - Spring with viscous drag
x0 = [5; 0];
N = 100;
x = zeros(2,N+1);
u = zeros(1,N+1);
x(:,1) = x0;
for i = 2:(N+1)
    x(:,i) = fS(x(:,i-1),u(:,i-1));
end
xf = [0;0];
Q = [1 0; 0 1];
Qf = 10*Q;
R = 1;

dt = .1;
W = [0 0; 0 0]*dt;

[x,u,runs,error] = ilqrIter(@fS,@dA,@dB,Q,R,Qf,N,x,u,xf,W,2,1e-5,20,1);

figure(5)
ti = 1:N+1;
hold on
plot(ti, x(1,:), "LineWidth", 1.3)
plot(ti, x(2,:), "LineWidth", 1.3)
plot(ti, u, "LineWidth", 1.3)
legend("x1","x2","u");
xlabel("Time");
ylabel("Displacement");

%% iLQR example - Inverted Pendulum
x0 = [pi; 0];
N = 30;
x = zeros(2,N+1);
u = zeros(1,N+1);

x = zeros(2,N+1);
u = zeros(1,N+1);
x(:,1) = x0;
for i = 2:(N+1)
    x(:,i) = fP(x(:,i-1),u(:,i-1));
end
xf = [0;0];
Q = [1 0; 0 0];
Qf = 10*[1 0; 0 1];
R = 1;

dt = .1;
W = [0 0; 0 0]*dt;

figure(6)
[x,u,runs,error] = ilqrIter(@fP,@dAP,@dBP,Q,R,Qf,N,x,u,xf,W,1,1e-5,150,1);
ti = 1:N+1;
hold on;
plot(ti, x(1,:), "LineWidth", 1.3)
%plot(ti, u, "LineWidth", 1.3)


%% 
% Spring with Viscous Drag System Function and Derivatives
% 0 = x2' + bx2^2 + kx1 - F + u
% x1' = x1 + 1/t*(x2)
% x2' = x2 + 1/t*(F - bx2*|x2| - kx1 + u) + *w*/t
function[xn] = fS(x,u)
    k = .5;
    b = 0.5;
    dt = .1;
    xn(1) = (x(1) + x(2)*dt);
    xn(2) = (x(2) + (0 - b*x(2)*abs(x(2))-k*x(1) + u)*dt);
end
function[A] = dA(x,u)
    k = .5;
    b = 0.5;
    dt = .1;
    A = [1 1*dt; -k*dt 1-(b*2*abs(x(2)))*dt];
end
function[B] = dB(x,u)
    k = .5;
    b = 0.5;
    dt = .1;
    B = [0;1*dt];
end


% Inverted Pendulum
% 0 = v2' - g/l*sin(v1) + u   
% v1' = v1 + v2/t
% v2' = v2 + g/l*sin(v1)/t +u/t
function[x] = fP(x,u)
    g = 9.8;
    l = 1;
    dt = .1;
    x(1) = (x(1) + x(2)*dt);
    x(2) = (x(2) + g/l*sin(x(1))*dt + u*dt);
end
function[A] = dAP(x,u)
    g = 9.8;
    l = 1;
    dt = .1;
    A = [1 1*dt; g/l*cos(x(1)) 1];
end
function[B] = dBP(x,u)
    g = 9.8;
    l = 1;
    dt = .1;
    B = [0;1*dt];
end