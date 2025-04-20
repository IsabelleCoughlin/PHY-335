% Use the Gauss-Seidel method with SOR to solve a simple 2D 
% potential problem with Dirichlet boundary conditions.
% This simple example defines the computational grid inside of the
% main program 'PotentialSOR_simple'
%
% Subroutines are used to fill the sourse function 'rho',
% to set the boundary conditions on the computational grid,
% to compute the residual, and to compute the Jacobi update
%
function PotentialSOR_simple_1(Nx,Ny,eps)
    % Set the coordinate
    Xmin = 0.0;
    Xmax = 1.5;
    Ymin = 0.0;
    Ymax = 1.0;
    dx = (Xmax-Xmin)/(Nx-1);
    dy = (Ymax-Ymin)/(Ny-1);
    x = Xmin:dx:Xmax;
    y = Ymin:dy:Ymax;
    [yg,xg]=meshgrid(y,x);
    % Initialize solution and source grids
    v = zeros(Nx,Ny);
    rho = source(x,y);
    figure(3)
    surface(xg,yg,rho,'EdgeAlpha',0.25);
    view(-35,45)
    % set the Dirichlet BCs
    v = DirichletBC(v,x,y);
    figure(1)
    cla
    surface(xg,yg,v,'EdgeAlpha',0.25);
    view(-35,45)
    %
    % SOR Loop
    res = residual(v,rho,dx,dy);
    figure(2)
    cla
    surface(xg,yg,res,'EdgeAlpha',0.25);
    view(-35,45)
    omega = 1.0; % Simple Gauss-Seidel relaxation
%     omega = 1.9; % 1.9 is quite good for this problem
    for l = 1:10000
        v = SOR(v,rho,dx,dy,omega);
        res = residual(v,rho,dx,dy);
        resnorm = norm(res)/sqrt(dx*dy);
        fprintf('|res| = %g\n',resnorm);
        figure(1)
        cla
        surface(xg,yg,v,'EdgeAlpha',0.25);
        view(-35,45)
        figure(2)
        cla
        surface(xg,yg,res,'EdgeAlpha',0.25);
        view(-35,45)
        pause(0.01)
        if resnorm < eps
            fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm)
            break
        end
    end
    figure(4)
    contour(xg,yg,v)
    [Ey,Ex] = gradient(v,dy,dx);
    Ex = - Ex;Ey = -Ey;
    figure(5)
    quiver(xg,yg,Ex,Ey)
end

% This subroutine defines the source function as a function of the 
% grid coordinates.  'x' is a 1-D array with the x coordinate values
% on the Cartesial grid.  'y' is a 1-D array with the y coordinate
% values on the grid, and the subroutine returns the source function
% on the computational grid.
function rho = source(x,y)
    nx = length(x);
    ny = length(y);
    rho = zeros(nx,ny);
    for i=1:nx
        for j=1:ny
            rho(i,j) = -300.0*exp(-20.0*((x(i)-0.75)^2 + (y(j)-0.5)^2))*sin(4.0*pi*x(i)/3);
        end
    end
end

% This subroutine sets the Dirichlet boundary conditions on the
% boundary of the computational grid.  'x' is a 1-D array with 
% the x coordinate values on the Cartesial grid.  'y' is a 1-D 
% array with the y coordinate values on the grid, and v is the 
% solution grid.
function v = DirichletBC(v,x,y)
    v(:,1) = 4.0*x/3.0;
    v(:,length(y)) = 2.0*(x - 1.0);
    v(1,:) = -2.0*y;
    v(length(x),:) = 2.0 - y;
end

% This subroutine computes the residual associated with the current
% approximate solution 'v'.  'rho' is the source function, and 'dx'
% and 'dy' are the x and y grid spacings.
function res = residual(v,rho,dx,dy)
    vol = dx*dy;
    rx = dy/dx;
    ry = dx/dy;
    res = - 2.0*(rx+ry)*v - vol*rho;
    [nx,ny] = size(v);
    res(:,1) = 0.0; res(:,ny)=0.0;res(1,:)=0.0;res(nx,:)=0.0;
    for i=2:nx-1
        for j=2:ny-1
            res(i,j) = res(i,j) + rx*(v(i+1,j)+v(i-1,j)) + ry*(v(i,j+1)+v(i,j-1));
        end
    end
end

% This subroutine computes the Gauss-Seidel correction to the 
% approximate solution v.  'omega' is the over-relaxation parameter.
% 'omega'=1 yields a standard Gauss-Seidel correction to the grid.
% 1<'omega'<2 yields an SOR correction to the grid.
function v = SOR(v,rho,dx,dy,omega)
    vol = dx*dy;
    rx = dy/dx;
    ry = dx/dy;
    [nx,ny] = size(v);
    for i=2:nx-1
        for j=2:ny-1
            delta = rx*(v(i+1,j)-2.0*v(i,j) + v(i-1,j)) + ry*(v(i,j+1)-2.0*v(i,j) +v(i,j-1)) - vol*rho(i,j);
            v(i,j) = v(i,j) - omega*delta/(-2.0*(rx+ry));
        end
    end
end