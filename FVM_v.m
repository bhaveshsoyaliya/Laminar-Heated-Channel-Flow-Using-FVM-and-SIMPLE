function FVM_v(Nx,Ny,dx,dy,dz,rho,Dx,Dy,Iv,jF,alphaU,u,v,pStar)

global Fw Fe Fn Fs DF aE aW aN aS aP bP dV

% Momentum equation for y-velocity component

% Setup coefficients
    Fw(Iv,jF) = rho*dy*dz*0.5*(u(Iv-1,jF) + u(Iv-1,jF+1));
    Fe(Iv,jF) = rho*dy*dz*0.5*(u(Iv  ,jF) + u(Iv,  jF+1));
    Fs(Iv,jF) = rho*dx*dz*0.5*(v(Iv  ,jF) + v(Iv,  jF-1));
    Fn(Iv,jF) = rho*dx*dz*0.5*(v(Iv  ,jF) + v(Iv,  jF+1));
    
    % Convective flux using hybrid upwinding
    for i = Iv
        for j = jF
            aW(i,j) = max([ Fw(i,j),(Dx + 0.5*Fw(i,j)),0]);
            aE(i,j) = max([-Fe(i,j),(Dx - 0.5*Fe(i,j)),0]);
            aS(i,j) = max([ Fs(i,j),(Dy + 0.5*Fs(i,j)),0]);
            aN(i,j) = max([-Fn(i,j),(Dy - 0.5*Fn(i,j)),0]);
        end
    end

    aE(Nx+1,:) = 0;     % Right boundary (outlet, dv/dx = 0)

    DF(Iv,jF) = Fe(Iv,jF) - Fw(Iv,jF) + Fn(Iv,jF) - Fs(Iv,jF);    
    
    aP(Iv,jF) = aE(Iv,jF) + aW(Iv,jF) + aN(Iv,jF) + aS(Iv,jF) + DF(Iv,jF);
    bP(Iv,jF) = dx*dz*(pStar(Iv,jF) - pStar(Iv,jF+1));
    dV(Iv,jF) = (dx*dz)./(aP(Iv,jF)/alphaU);

end