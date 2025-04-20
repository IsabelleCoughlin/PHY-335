%Nx = 20; % Number of points on x
%Ny = 20; % Number of points on y
clear
T = 1; % LENGTH of C

Nc = 4;
dx = T/(Nc-1);

W = 2*T;
Nb = Nc*2;
dy = W/(Nb-1);

Na = 4*Nc;
La = dx*(Na - 1);
Ne = 4*Nc;
Le = dx*(Ne-1);
Lg = Le + T + La;
Nx = Ne + Nc + Na-2;
Ng = Nx;


Ny = 5*Nc;
Lf = dy*(Ny - 1);


% Nx = Ne + Nc + Na
% Ny = Nf


eps = 1.e-3;

% Set the coordinate
% Nc = 1;
% Nb = 3;
% Na = 6;
% Ne = 5;
% Nf = 9;
% dx = (Ne + Nc + Na)/(Nx-1);
% dy = (Nf)/(Ny-1);

% Initialize solution and source grids
psi = zeros(Ny,Nx);
xi = zeros(Ny, Nx);
res_x = zeros(Ny,Nx);
res_p = zeros(Ny, Nx);
v0 = 1;
nu = 0.5;
for j = 1:Ny
    for i = 1:Nx
        %psi(j, :) = v0*y(j); 
        psi(i, j) = v0*dy*(j-1)/Ny;
    end
end
omega = 1.0; % Simple Gauss-Seidel relaxation
test = zeros(Ny,Nx);

test = test_psi_boundary_conditions(test, Nx, Ny, Na, Nb, Nc, Ne, Ny, v0, dx, dy);

[psi, xi] = boundary_conditions(psi, xi, Nx, Ny, Na, Nb, Nc, Ne, Ny, v0, dx, dy);
figure(1)
cla
contour(psi);


for l = 1:10000

    [psi, res_p] = SOR_Psi(psi, xi, dx, dy, omega, Na, Nb, Nc, Ne, Ny, Nx, res_p);
    [xi, res_x] = SOR_Xi(xi, psi, dx, dy, omega, Na, Nb, Nc, Ne, Ny, Nx, nu, res_x);
    
    res_p = abs(res_p);
    res_x = abs(res_x);
    max_p = max(res_p);
    
    max_x = max(res_x);
    resnorm_p = norm(max_p)/sqrt(dx*dy);
    resnorm_x = norm(max_x)/sqrt(dx*dy);

    figure(1)
    cla
    contour(psi);

    figure(2)
    cla
    contour(xi);
    
    % figure(3)
    % cla
    % contour(res_p);
    % 
    % figure(4)
    % cla
    % contour(res_x);

    fprintf('PSI: at l=%d, |res| = %g\n',l,resnorm_p)
    fprintf('XI: at l=%d, |res| = %g\n',l,resnorm_x)

    

    if (resnorm_x < eps) && (resnorm_p < eps)
        fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm_p)
        fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm_x)
        break
    end
   
end

function [psi, xi] = boundary_conditions(psi, xi, Nx, Ny, Na, Nb, Nc, Ne, Nf, v0, dx, dy)

    % All things not touching the boundary



    % Left side
    xi(:, 1) = 0;
    psi(:, 1) = psi(:, 2);

    % Left Test


    % Top
    xi(Nf, :) = 0;
    psi(Nf, :) = v0*dy + psi(Nf-1, :);
    %psi(Nf, :) = v0;

    % Right
    xi(:, end) = xi(:, end-1);
    psi(:, end) = psi(:, end-1);

    % For things touching the boundary

    % E
    for i = 1:Ne-1
        xi(1,i) = 0;
        psi(1, i) = 0;
    end

    % D
    for j = 1:Nb
        psi(j, Ne) = 0;
        psi(j, Ne) = psi(j, Ne-1);
        xi(j, Ne) = (-2*(psi(j, Ne-1)))/(dx^2);
    end

    % C
    for i = Ne+1:Ne+Nc-1
        psi(Nb, i) = xi(Nb, i-1);
        psi(Nb, i) = 0;
        xi(Nb,i) = (-2*(psi(Nb+1, i)))/(dy^2);
    end

    % B
    for j = 1:Nb
        psi(j, Ne+Nc) = psi(j, Ne + Nc + 1);
        psi(j, Ne + Nc) = 0;
        xi(j, Ne + Nc) = (-2*(psi(j, Ne+Nc+1)))/(dx^2);
    end

    % A
    for i = (Ne + Nc):(Nx)
        xi(1, i) = 0;
        psi(1,i) = 0;
    end
end

function [f,res_psi]  = SOR_Psi(f, source, dx, dy, omega, Na, Nb, Nc, Ne, Nf, Nx, res_psi)
    vol = dx*dy; % good
    rx = dy/dx; % good
    ry = dx/dy; % good
    
    % f is psi, source is xi (or i name them that?
    
    % Psi first
    for i=2:Ne-1 
        for j=2:(Nf-1)
            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            delta = omega * res / (-2*(rx + ry));
            %res_psi = max(res_psi, abs(res));
            res_psi(j,i) = res;
            f(j,i) = f(j,i) - delta;
        end
    end
    
    for i = Ne:(Ne+Nc)
        for j = (Nb + 1): (Nf-1)
            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            %res_psi = max(res_psi, abs(res));
            delta = omega * res / (-2*(rx + ry));
            res_psi(j,i) = res;
            f(j,i) = f(j,i) - delta;
        end
    end
    for i = (Ne + Nc):(Nx)
        for j = 2:(Nf-1)
            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            %res_psi = max(res_psi, abs(res));
            delta = omega * res / (-2*(rx + ry));
            res_psi(j,i) = res;
            f(j,i) = f(j,i) - delta;
        end
    end
end


function [xi,res_xi]  = SOR_Xi(xi, psi, dx, dy, omega, Na, Nb, Nc, Ne, Nf,Nx, nu, res_xi)
    vol = dx*dy; % good
    rx = dy/dx; % good
    ry = dx/dy; % good
    
    
    
    for i=2:Ne-1 
        for j=2:(Nf-1)
            
            dpsi_dx = (psi(j,i+1) - psi(j,i-1)) / (2*dx);
            dpsi_dy = (psi(j+1,i) - psi(j-1,i)) / (2*dy);
            dxi_dx  = (xi(j,i+1) - xi(j,i-1)) / (2*dx);
            dxi_dy  = (xi(j+1,i) - xi(j-1,i)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(j,i+1) - 2*xi(j,i) + xi(j,i-1)) + ry*(xi(j+1,i) - 2*xi(j,i) + xi(j-1,i)) - vol * source;
            %res_xi = max(res_xi, abs(res));
            delta = omega * res / (-2*(rx + ry));
            res_xi(j,i) = res;
            xi(j,i) = xi(j,i) - delta;
        end
    end
    
    for i = Ne:(Ne+Nc)
        for j = (Nb + 1): (Nf-1)
            dpsi_dx = (psi(j,i+1) - psi(j,i-1)) / (2*dx);
            dpsi_dy = (psi(j+1,i) - psi(j-1,i)) / (2*dy);
            dxi_dx  = (xi(j,i+1) - xi(j,i-1)) / (2*dx);
            dxi_dy  = (xi(j+1,i) - xi(j-1,i)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(j,i+1) - 2*xi(j,i) + xi(j,i-1)) + ry*(xi(j+1,i) - 2*xi(j,i) + xi(j-1,i)) - vol * source;
            %res_xi = max(res_xi, abs(res));
            delta = omega * res / (-2*(rx + ry));
            res_xi(j,i) = res;
            xi(j,i) = xi(j,i) - delta;
        end
    end
    for i = (Ne + Nc):(Nx)
        for j = 2:(Nf-1)
            dpsi_dx = (psi(j,i+1) - psi(j,i-1)) / (2*dx);
            dpsi_dy = (psi(j+1,i) - psi(j-1,i)) / (2*dy);
            dxi_dx  = (xi(j,i+1) - xi(j,i-1)) / (2*dx);
            dxi_dy  = (xi(j+1,i) - xi(j-1,i)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(j,i+1) - 2*xi(j,i) + xi(j,i-1)) + ry*(xi(j+1,i) - 2*xi(j,i) + xi(j-1,i)) - vol * source;
            %res_xi = max(res_xi, abs(res));
            delta = omega * res / (-2*(rx + ry));
            res_xi(j,i) = res;
            xi(j,i) = xi(j,i) - delta;
        end
    end
end


% function [psi] = test_xi_boundary_conditions(test, Nx, Ny, Na, Nb, Nc, Ne, Nf, v0, dx, dy)
% 
%     % All things not touching the boundary
% 
% 
% 
%     % Left side
%     xi(:, 1) = 0;
% 
% 
%     % Left Test
% 
% 
%     % Top
%     xi(Nf, :) = 0;
%     %psi(Nf, :) = v0*dy + psi(Nf-1, :);
% 
%     % Right
%     xi(:, end) = xi(:, end-1);
% 
%     % For things touching the boundary
% 
%     % E
%     for i = 1:Ne-1
%         xi(1,i) = 0;
%         psi(1, i) = 0;
%     end
% 
%     % D
%     for j = 1:Nb
%         psi(j, Ne) = 0;
%         psi(j, Ne) = psi(j, Ne-1);
%         xi(j, Ne) = (-2*(psi(j, Ne-1)))/(dx^2);
%     end
% 
%     % C
%     for i = Ne+1:Ne+Nc-1
%         psi(Nb, i) = xi(Nb, i-1);
%         psi(Nb, i) = 0;
%         xi(Nb,i) = (-2*(psi(Nb+1, i)))/(dy^2);
%     end
% 
%     % B
%     for j = 1:Nb
%         psi(j, Ne+Nc) = psi(j, Ne + Nc + 1);
%         psi(j, Ne + Nc) = 0;
%         xi(j, Ne + Nc) = (-2*(psi(j, Ne+Nc+1)))/(dx^2);
%     end
% 
%     % A
%     for i = (Ne + Nc + 1):(Ne + Nc + Na)
%         xi(1, i) = 0;
%         psi(1,i) = 0;
%     end
% end

function test = test_psi_boundary_conditions(test, Nx, Ny, Na, Nb, Nc, Ne, Nf, v0, dx, dy)

    % All things not touching the boundary

    
    
    % Left side
    
    %psi(:, 1) = psi(:, 2);
    test(:, 1) = 1;

    % Left Test

    
    % Top
    test(Nf, :) = 2;
    
    % Right
    
    test(:, end) = 3;
    
    % For things touching the boundary
    
    % E
    for i = 1:Ne-1
        
        test(1, i) = 4;
    end
    
    % D
    for j = 1:Nb
        test(j, Ne) = 5;
        
    end
    
    % C
    for i = Ne+1:Ne+Nc-1
        test(Nb, i) = 6;
    end
    
    % B
    for j = 1:Nb
        test(j, Ne+Nc) = 7;
        
    end
    
    % A
    for i = (Ne + Nc):(Nx)
        test(1,i) = 8;
    end
end
