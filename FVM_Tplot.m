function FVM_Tplot(Nx,Ny,x,y,Lx,Ly,rho,cp,u,T,um,TW,TN,qN,BC)

% Calculate mean, surface, dimensionless, and analytical FD temperatures

uCD = (u(1:Nx,:) + u(2:Nx+1,:))/2;
Tm = sum((uCD.*T),2)/(Ny*um);

if BC == 0
    Ts   = TN*ones(Nx,1); 
else
    Tm_e = (TW + qN*x/(rho*um*Ly*cp));
    Ts   = (1.5*T(:,Ny) - 0.5*T(:,Ny-1))';
end

Ni = 6;
Nc = floor(linspace(1,Nx,Ni));
Tstar = zeros(Ni,Ny);
for i = 1:length(Nc)
    I = Nc(i);
    Tstar(i,:) = (Ts(I) - T(I,:))./(Ts(I) - Tm(I));
end

ystar   = [0:0.05:1];
TstarFD = (35/136)*(5-6*(ystar).^2+(ystar).^4);

%% Print out results

fprintf('Temperature matrix: \n\n')
for J = Ny:-1:1
    fprintf('%8.2f', T(:,J))
    fprintf('\n')
end
fprintf('\n\n')

%% Plot out results

% Contour plot of temperatures

if BC == 0
    figure('Name','Constant Surface Temperature',...
           'Position',[300 300 500 500])
else
    figure('Name','Constant Surface Heat Flux',...
           'Position',[300 300 500 500])
end

subplot(2,2,[1,2]);
y_sym = [-fliplr(y) y];
T_sym = [flipud(T') ; T'];

contourf(x,y_sym,T_sym)
colormap 'jet'
colorbar
title('Contour Lines of Temperature (\circC)','FontSize',12)
xlabel('\itx \rm(m)','FontSize',12)
ylabel('\ity \rm(m)','FontSize',12)
axis equal

% Mean and surface temperatures versus x

subplot(2,2,3);
hold on
plot([0,x],[TW,Tm'],'--ro','DisplayName','mean')
if BC == 0
    plot([0,Lx],[TN,TN], '--k', 'DisplayName','surface',...
                                    'LineWidth',2)    
else
    plot([0,x],[TW,Tm_e],'--k', 'DisplayName','mean (exact)',...
                                    'LineWidth',2)
    plot([0,x],[TW,Ts  ],'--bo','DisplayName','surface')
end
title('Axial Temperature Variation')
xlabel('\itx \rm(m)','FontSize',12)
ylabel('Temperature (\circC)','FontSize',12)
legend('Location','northwest')
hold off

% Dimensionless temperature versus y/Ly at several x/Dh locations

subplot(2,2,4);
style = ["-bo","-bx","-b*","-b+","-bx","-bs"];
hold on
for i = 1:Ni
    plot(Tstar(i,:),y/Ly,style(i),...
         'DisplayName',num2str(x(Nc(i))/(2*Ly),'%3.1f'))
end
plot(TstarFD,ystar,'--k','DisplayName','FD','LineWidth',2)
title('Dimensionless Temperature Profiles')
xlabel('(\itT - T_s \rm) / (\itT_m - T_s \rm)','FontSize',12)
ylabel('\ity \rm/ \itL_y \rm','FontSize',12)
lgd = legend('Location','southwest');
title(lgd,'\itx \rm/ \itD_h')
hold off

end