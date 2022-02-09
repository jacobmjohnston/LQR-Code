function[K,Kv,Ku,v] = ilqrSolve(fA,fB,Q,R,Qf,N,x,u,xf)
    %Uses algorithm described in "Iterative Linear Quadratic Regulator Design 
    % for Nonlinear Biological Movement Systems" Li, Todorov

    [n,~] = size(Q);
    [m,~] = size(R);
    P = zeros(n,n,N);
    v = zeros(n,N);
    K = zeros(m,n,N);
    Ku = zeros(m,m,N);
    Kv = zeros(m,n,N);
    P(:,:,N+1) = Qf;
    v(:,N+1) = Qf*(x(:,N+1)-xf);
    for i = linspace(N,1,N)
        A = fA(x(:,i),u(:,i));
        B = fB(x(:,i),u(:,i));
        C = (B'*P(:,:,i+1)*B+R)^(-1);
        E = P(:,:,i+1)*A;
        Ku(:,:,i) = C*R;
        Kv(:,:,i) = C*B';
        K(:,:,i) = Kv(:,:,i)*E;
        D = (A-B*K(:,:,i));
        v(:,i) = D'*v(:,i+1)-K(:,:,i)'*R*u(:,i)+Q*x(:,i);
        P(:,:,i) = E'*D+Q;
    end
end