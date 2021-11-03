function [x,u] = lqrState(x0,K,A,B,W,N)
    x = zeros(length(x0),N+1);
    [~,ul] = size(B);
    u = zeros(ul,N+1);
    x(:,1) = x0;
    for i = 2:(N+1)
        u(:,i-1) = K(:,:,i-1)*x(:,i-1);
        x(:,i) = A*x(:,i-1) + B*u(:,i-1) + randomN(0,W);
    end
end

