function FVM_pcorr(Nx,Ny,dx,dy,dz,rho,Ip,Jp,uStar,vStar)

global aE aW aN aS aP bP dU dV

    % Setup coefficients
    aE(Ip,Jp) = rho*dU(Ip,  Jp  )*dy*dz;
    aW(Ip,Jp) = rho*dU(Ip-1,Jp  )*dy*dz;
    aN(Ip,Jp) = rho*dV(Ip,  Jp  )*dx*dz;
    aS(Ip,Jp) = rho*dV(Ip,  Jp-1)*dx*dz;

    aP(Ip,Jp) = aE(Ip,Jp) + aW(Ip,Jp)+ aN(Ip,Jp) + aS(Ip,Jp);
    
    bP(Ip,Jp) = rho*(uStar(Ip-1,Jp) - uStar(Ip,Jp))*dy*dz...
              + rho*(vStar(Ip,Jp-1) - vStar(Ip,Jp))*dx*dz;
    
    % Fix pressure at outlet to zero
    i = Nx+1;
    aE(i,Jp) = 0.0; aW(i,Jp) = 0.0;
    aN(i,Jp) = 0.0; aS(i,Jp) = 0.0;
    aP(i,Jp) = 1.0; bP(i,Jp) = 0.0;  

end