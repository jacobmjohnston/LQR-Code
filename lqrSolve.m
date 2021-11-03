function [K] = lqrSolve(A,B,Q,R,Qf,N)
    [n,m] = size(Q);
    P = zeros(n,m,N+1);
    [xl,~] = size(A);
    [~,ul] = size(B);
    K = zeros(ul,xl,N);
    P(:,:,N+1) = Qf;
    for i = linspace(N,1,N)
        Bi = B'*P(:,:,i+1)*B;
        ABi = A'*P(:,:,i+1)*B;
        K(:,:,i) = -(R+Bi)^(-1)*ABi';
        P(:,:,i) = Q + A'*P(:,:,i+1)*A + ABi*K(:,:,i);
    end
end