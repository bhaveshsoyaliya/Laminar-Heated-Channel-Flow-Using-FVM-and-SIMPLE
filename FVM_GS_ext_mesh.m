function [phi,res_sum] = FVM_GS_ext_mesh(Nx,Ny,alpha,Nit_max,tol,phi)

global aW aE aS aN aP bP

% Function that uses Gauss Seidel iteration with relaxation to solve:
%   (a*phi)_P = sum(a*phi)_nb + b

res = zeros(Nx-1,Ny-1);

for Nit = 1:Nit_max
    
    for J = 2:Ny
        for I = 2:Nx 
            
            phi(I,J) = (aW(I,J)*phi(I-1,J) + aE(I,J)*phi(I+1,J) + ...
                        aS(I,J)*phi(I,J-1) + aN(I,J)*phi(I,J+1) + ...
                        bP(I,J))/(aP(I,J)/alpha) + (1 - alpha)*phi(I,J);
        end
    end
    
    for J = 2:Ny
        for I = 2:Nx
          res(I-1,J-1) = (aW(I,J)*phi(I-1,J) + aE(I,J)*phi(I+1,J) + ...
                          aS(I,J)*phi(I,J-1) + aN(I,J)*phi(I,J+1) + bP(I,J))...
                        - aP(I,J)*phi(I,J);
        end
    end
    
    res_sum = sum(sum(abs(res)))/...
              sum(sum(abs(aP(2:Nx,2:Ny).*phi(2:Nx,2:Ny))));
    if res_sum < tol
        break
    end
    
end

end
