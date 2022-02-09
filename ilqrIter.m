function[x,u,i,err] = ilqrIter(f,fA,fB,Q,R,Qf,N,x,u,xf,method,errmax,imax,dt)
    err = 1;
    i = 1;
    while(err > errmax && i < imax)
        %Gets contorl feedback matrices 
        [K,Kv,Ku,v] = ilqrSolve(fA,fB,Q,R,Qf,N,x,u,xf);
        
        %Progataes new trajectory with new control input
        [xn,un] = ilqrState(K,Ku,Kv,v,f,fA,fB,x,u,N,method,dt);
        
        %error terms
        xe = abs(x - xn);
        ue = abs(u(:,1:N)-un(:,1:N));
        err = max(max(xe(:)),max(ue(:)));
        
        
        x = xn;
        u = un;
        i = i+1;
    end
end
