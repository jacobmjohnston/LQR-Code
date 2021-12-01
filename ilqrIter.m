function[x,u,i,err] = ilqrIter(f,fA,fB,Q,R,Qf,N,x,u,xf,W,method,errmax,imax,a)
    err = 1;
    i = 1;
    while(err > errmax && i < imax)
        [K,Kv,Ku,v] = ilqrSolve(fA,fB,Q,R,Qf,N,x,u,xf);
        [xn,un] = ilqrState(K,Ku,Kv,v,f,fA,fB,x,u,N,method,a);
        xe = abs(x - xn);
        ue = abs(u(:,1:N)-un(:,1:N));
        err = max(max(xe(:)),max(ue(:)));
        x = xn;
        u = un;
        i = i+1;
    end
    if(max(W(:) > 0))
        for j = 1:N
            w = randomN(0,W);
            y = f(x(:,j),u(:,j));
            x(:,j+1) = y(:) + w;
        end
    end
end
