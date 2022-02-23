function[K,j] = ddpSolve(df,dl,dlf,N,x,u,xf)
    n = length(xf);
    [m,~] = size(u);
    K = zeros(m,n,N);
    j = zeros(m,N);
    
    A = dlf(x(:,end),'xx');
    b = dlf(x(:,end),'x');
    
    for i = linspace(N,1,N)
        dfx = df(x(:,i),u(:,i),'x');
        dfu = df(x(:,i),u(:,i),'u');
        Qx = dl(x(:,i),u(:,i),'x') + dfx'*b;
        Qu = dl(x(:,i),u(:,i),'u') + dfu'*b;
        
        % calculates the product of a vector dot 3rd order tensor 
        dfxx = df(x(:,i),u(:,i),'xx');
        [~,~,o] = size(dfxx);
        c = zeros(n,n);
        for p = 1:o
            c = c + b'*dfxx(:,:,p);
        end
        
        Qxx = dl(x(:,i),u(:,i),'xx') + dfx'*A*dfx + c;
        Quu = dl(x(:,i),u(:,i),'uu') + dfu'*A*dfu + b'*(df(x(:,i),u(:,i),'uu'));
        Qux = dl(x(:,i),u(:,i),'ux') + dfu'*A*dfx + b'*(df(x(:,i),u(:,i),'ux'));
        
        K(:,:,i) = -(Quu^-1)*Qux;
        j(:,i) = -(Quu^-1)*Qu;
        
        A = Qxx + K(:,:,i)'*Quu*K(:,:,i) + Qux'*K(:,:,i) + K(:,:,i)'*Qux;
        b = Qx + K(:,:,i)'*Quu*j(:,i) + Qux'*j(:,i) + K(:,:,i)'*Qu;
    end
end