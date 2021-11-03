%% Deterministic Example
A = [1 1; 0 1];
B = [0;1];
N = 20;
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

%% Stochastic Example
Q = diag(ones(1,5));
Qf = 10*diag(ones(1,5));
R = diag(ones(1,2));
N = 30;
A = [1 1 0 0 0; 0 1 1 0 0; 0 0 1 0 1; 0 1 0 1 0; 0 0 0 0 1];
B = [0 0; 1 0; 0 1; 0 0; 1 1];

K = lqrSolve(A,B,Q,R,Qf,N);
W = diag(ones(1,5))/2;
[x,u] = lqrState([1 0 0 0 0], K, A, B, W, N);

y = [1 0 0 0 0]*x;
t = 1:(N+1);
plot(t,y);
