function FVM_u(Nx,Ny,dx,dy,dz,rho,Dx,Dy,iF,Ju,alphaU,u,v,pStar,BC_S_MO)

global Fw Fe Fn Fs DF aE aW aN aS aP bP dU

    % Momentum equation for x-velocity component

    % Setup coefficients
    Fw(iF,Ju) = rho*dy*dz*0.5*(u(iF,Ju)   + u(iF-1,Ju  ));
    Fe(iF,Ju) = rho*dy*dz*0.5*(u(iF,Ju)   + u(iF+1,Ju  ));
    Fs(iF,Ju) = rho*dx*dz*0.5*(v(iF,Ju-1) + v(iF+1,Ju-1));
    Fn(iF,Ju) = rho*dx*dz*0.5*(v(iF,Ju)   + v(iF+1,Ju  ));
     
    % Convective flux using hybrid upwinding
    for i = iF
        for j = Ju
            aW(i,j) = max([ Fw(i,j),(Dx + 0.5*Fw(i,j)),0]);
            aE(i,j) = max([-Fe(i,j),(Dx - 0.5*Fe(i,j)),0]);
            aS(i,j) = max([ Fs(i,j),(Dy + 0.5*Fs(i,j)),0]);
            aN(i,j) = max([-Fn(i,j),(Dy - 0.5*Fn(i,j)),0]);
        end
    end
    
    if BC_S_MO == 0     % Wall or symmetry at bottom
        aS(:,2) = 2*Dy;
    else
        aS(:,2)    = 0;     
    end
    aN(:,Ny+1) = 2*Dy;  % Wall at top
    aE(Nx+1,:) = 0;     % Right boundary (outlet, du/dx = 0)
    
    DF(iF,Ju) = Fe(iF,Ju) - Fw(iF,Ju) + Fn(iF,Ju) - Fs(iF,Ju);
    
    aP(iF,Ju) = aE(iF,Ju) + aW(iF,Ju) + aN(iF,Ju) + aS(iF,Ju) + DF(iF,Ju);
    bP(iF,Ju) = dy*dz*(pStar(iF,Ju) - pStar(iF+1,Ju));
    dU(iF,Ju) = (dy*dz)./(aP(iF,Ju)/alphaU);

end