function [xn,un] = ilqrState(K,Ku,Kv,v,f,fA,fB,x,u,N,method,a)
    dx = zeros(length(x(:,1)),N+1);
    du = zeros(length(u(:,1)),N+1);
    
    if(method == 2)
        xn = zeros(length(x(:,1)),N+1);
        un = zeros(length(u(:,1)),N+1);
        xn(:,1) = x(:,1);
        un(:,1) = u(:,1);
    end
    
    for i = 2:(N+1)
        if(method == 2)
            A = fA(xn(:,i-1),un(:,i-1));
            B = fB(xn(:,i-1),un(:,i-1));
        end
        if (method == 1)
            A = fA(x(:,i-1),u(:,i-1));
            B = fB(x(:,i-1),u(:,i-1));
        end
        du(:,i-1) = -K(:,:,i-1)*dx(:,i-1) -Kv(:,:,i-1)*v(:,i) - Ku(:,:,i-1)*u(:,i-1);
        dx(:,i) = A*dx(:,i-1) + B*du(:,i-1);
        
        if(method == 2)
            un(:,i-1) = u(:,i-1) + a*du(:,i-1);
            xn(:,i) = f(xn(:,i-1),un(:,i-1));
        end
    end
    if(method == 1)
        xn = x + a*dx;
        un = u + a*du;
    end
end