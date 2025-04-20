% Use the Grid Class to solve the Poisson equation using the
% discrete sine transformation
function PotentialSin(Nx,Ny,Nit)
    % Set the coordinate
    Xmin = 0.0;
    Xmax = 1.5;
    Ymin = 0.0;
    Ymax = 1.0;
    Fgrid = Grid(Xmin,Xmax,Ymin,Ymax,Nx,Ny);
    %
    % Solver Iteration
    for l = 1:Nit
        Fgrid.SinSolver();
        res = Fgrid.residual();
        resnorm = norm(res)/sqrt(Fgrid.dx*Fgrid.dy);
        %
%         fprintf('|res| = %g\n',resnorm);
%         figure(1)
%         cla
%         axis([Xmin Xmax Ymin Ymax -2 2.4]);  
%         Fgrid.plotv();
%         view(-35,45)
%         figure(2)
%         cla
%         surface(Fgrid.xg,Fgrid.yg,res,'EdgeAlpha',0.15);
%         view(-35,45)
        %
    end
    fprintf('|res| = %g\n',resnorm);
    figure(1)
    cla
    axis([Xmin Xmax Ymin Ymax -2 2.4]);  
    Fgrid.plotv();
    view(-35,45)
    figure(2)
    cla
    surface(Fgrid.xg,Fgrid.yg,Fgrid.residual(),'EdgeAlpha',0.15);
    view(-35,45)
    figure(3)
    cla
    surface(Fgrid.xg,Fgrid.yg,Fgrid.rho,'EdgeAlpha',0.15);
    view(-35,45)
    figure(4)
    cla
    contour(Fgrid.xg,Fgrid.yg,Fgrid.v)
    [Ey,Ex] = gradient(Fgrid.v,Fgrid.dy,Fgrid.dx);
    Ex = - Ex;Ey = -Ey;
    figure(5)
    cla
    quiver(Fgrid.xg,Fgrid.yg,Ex,Ey)
end