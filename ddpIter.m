function[x,u,i,err] = ddpIter(f,df,dl,dlf,N,x,u,xf,errmax,imax,a)
    err = 1;
    i = 1;
    while(err > errmax && i < imax)
        %Gets contorl feedback matrices 
        [K,j] = ddpSolve(df,dl,dlf,N,x,u,xf);
        
        [xn,un] = ddpState(f,K,j,x,u,N,a);
        
        xe = abs(x - xn);
        ue = abs(u(:,1:N)-un(:,1:N));
        err = max(max(xe(:)),max(ue(:)));
        
        x = xn;
        u = un;
        i = i+1;
    end
end