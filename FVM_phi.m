function FVM_phi(Nx,Ny,dx,dy,dz,rho,gamma,dphi,Ip,Jp,u,v,BC_S,BC_N)

global Fw Fe Fn Fs DF aE aW aN aS aP bP

% Coefficients for advection-diffusion equation for phi

% Setup diffusion coefficients
Dx = (gamma/dx)*(dy*dz);
Dy = (gamma/dy)*(dx*dz);

% Setup flow (or advection) coefficients
Fe(Ip,Jp) = rho*dy*dz*u(Ip  ,Jp  );
Fw(Ip,Jp) = rho*dy*dz*u(Ip-1,Jp  );
Fn(Ip,Jp) = rho*dx*dz*v(Ip  ,Jp  );
Fs(Ip,Jp) = rho*dx*dz*v(Ip  ,Jp-1);

% Convective flux using hybrid upwinding
for i = Ip
    for j = Jp
        aE(i,j) = max([-Fe(i,j),(Dx-0.5*Fe(i,j)),0]);
        aW(i,j) = max([ Fw(i,j),(Dx+0.5*Fw(i,j)),0]);
        aN(i,j) = max([-Fn(i,j),(Dy-0.5*Fn(i,j)),0]);
        aS(i,j) = max([ Fs(i,j),(Dy+0.5*Fs(i,j)),0]);
    end
end

bP(Ip,Jp) = 0;

if BC_N == 0                % Top: wall (no-slip) at Tw or qw
    aN(:,Ny+1) = 2*Dy;  
else
    aN(:,Ny+1) = 0;
    bP(:,Ny+1) = dphi*dx*dz;
end
if BC_S == 0                % Bottom: wall (no-slip) at Tw or symmetry
    aS(:,2)    = 2*Dy;      
else
    aS(:,2)    = 0;
end

aE(Nx+1,:) = 0;             % Right boundary (outlet, dT/dx = 0)

DF(Ip,Jp) = Fe(Ip,Jp) - Fw(Ip,Jp) + Fn(Ip,Jp) - Fs(Ip,Jp);
aP(Ip,Jp) = aE(Ip,Jp) + aW(Ip,Jp) + aN(Ip,Jp) + aS(Ip,Jp) + DF(Ip,Jp);

end