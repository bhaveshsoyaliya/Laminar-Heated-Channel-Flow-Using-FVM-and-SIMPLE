function FVM_Vplot(Nx,Ny,x,xu,y,Ly,u,v,p,um)

% Print out results

fprintf('\nx-velocity matrix: \n')
for J = Ny:-1:1
    fprintf('%8.4f', u(:,J))
    fprintf('\n')
end
fprintf('\n\n')

fprintf('\ny-velocity matrix: \n')
for J = Ny+1:-1:1
    fprintf('%8.4f', v(:,J))
    fprintf('\n')
end
fprintf('\n\n')

fprintf('\npressure matrix: \n')
for J = Ny:-1:1
    fprintf('%8.1e', p(:,J))
    fprintf('\n')
end
fprintf('\n\n')

%% Plot out results

% Contour plot of velocities

y_sym  = [-fliplr(y)  y];
u_sym  = [flipud(u') ; u'];
p_sym  = [flipud(p') ; p'];

figure('Name','Velocity and Pressure Contour Plots',...
       'Position',[200 200 500 500])

subplot(3,2,[1,2]);
contourf(xu,y_sym,u_sym); colormap 'jet'; colorbar
title('Contour Lines of x-Velocity (m/s)','FontSize',12)
xlabel('\itx \rm(m)','FontSize',12)
ylabel('\ity \rm(m)','FontSize',12)
axis equal

subplot(3,2,[3,4]);

contourf(x,y_sym,p_sym); colormap 'jet'; colorbar
title('Contour Lines of Pressure (Pa)','FontSize',12)
xlabel('\itx \rm(m)','FontSize',12)
ylabel('\ity \rm(m)','FontSize',12)
axis equal

% Dimensionless x-velocity versus y/Ly at several x/Dh locations

ystar   = 0:0.05:1;
ustarFD = (1 - (ystar).^2); % FD velocity profile

Ni = 6;
Nc = floor(linspace(1,Nx,Ni));
ustar = zeros(Ni,Ny);
for i = 1:Ni
    I = Nc(i);
    ustar(i,:) = u(I,:)/(3/2*um);
end

subplot(3,2,[5,6]);
style = ["-bo","-bx","-b*","-b+","-bx","-bs"];
hold on
for i = 1:Ni
    plot(ustar(i,:),y/Ly,style(i),...
         'DisplayName',num2str(x(Nc(i))/(2*Ly),'%3.1f'))
end
plot(ustarFD,ystar,'--k','DisplayName','FD','LineWidth',2)
title('Dimensionless Velocity Profiles')
xlabel('(\itu \rm) / (\itu_m \rm)','FontSize',12)
ylabel('\ity \rm/ \itL_y \rm','FontSize',12)
lgd = legend('Location','southwest');
title(lgd,'\itx \rm/ \itD_h')
hold off

end