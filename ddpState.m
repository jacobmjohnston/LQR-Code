function[xn,un] = ddpState(f,K,j,x,u,N,a)
    xn = zeros(length(x(:,1)),N+1);
    un = zeros(length(u(:,1)),N+1);
    xn(:,1) = x(:,1);

    for i = 2:(N+1)
        un(:,i-1) = u(:,i-1) + K(:,:,i-1)*(xn(:,i-1)-x(:,i-1)) + j(:,i-1);
        un(:,i-1) = a*un(:,i-1) + (1-a)*u(:,i-1);
        xn(:,i) = f(xn(:,i-1),un(:,i-1));
    end
end