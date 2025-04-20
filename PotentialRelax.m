% Use the Grid Class to solve the Poisson equation using various
% relaxation methods.
function PotentialRelax(Nx,Ny,eps)
    % Set the coordinate
    Xmin = 0.0;
    Xmax = 1.5;
    Ymin = 0.0;
    Ymax = 1.0;
    Fgrid = Grid(Xmin,Xmax,Ymin,Ymax,Nx,Ny);
    %
    % SOR Loop
    omega = 1.9; % 1.9 is quite good
    for l = 1:100000
%         Fgrid.Jacobi(omega);
%         Fgrid.SOR(omega);
        Fgrid.RedBlackSOR(omega);
        res = Fgrid.residual();
        resnorm = norm(res)/sqrt(Fgrid.dx*Fgrid.dy);
        %
        fprintf('|res| = %g\n',resnorm);
        figure(1)
        cla
        axis([Xmin Xmax Ymin Ymax -2 2.4]);  
        Fgrid.plotv();
        view(-35,45)
        figure(2)
        cla
        surface(Fgrid.xg,Fgrid.yg,res,'EdgeAlpha',0.15);
        view(-35,45)
        %
        if resnorm < eps
            fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm)
            break
        end
    end
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