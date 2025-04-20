function HW10(Nx,Ny,eps)
    % Set the coordinate
    Nc = 1;
    Nb = 3;
    Na = 5;
    Ne = 4;
    Nf = 10;
    dx = (Ne + Nc + Na)/(Nx-1);
    dy = (Nf)/(Ny-1);

    % Initialize solution and source grids
    psi = zeros(Ny,Nx);
    xi = zeros(Ny, Nx);
    v0 = 1;
    nu = 0.5;
    for j = 1:Ny
        for i = 1:Nx
            %psi(j, :) = v0*y(j); 
            psi(i, j) = v0*dy*(j-1)/Ny;
        end
    end
    omega = 1.0; % Simple Gauss-Seidel relaxation
  

    [psi, xi] = boundary_conditions(psi, xi, Nx, Ny, Na, Nb, Nc, Ne, Nf, v0, dx, dy);
    
    figure(1)
    cla
    contour(psi);

    % figure(2)
    % contour(xi);
    
   

    for l = 1:10000
   
        [psi, res_p] = SOR_Psi(psi, xi, dx, dy, omega, Na, Nb, Nc, Ne, Nf);
        [xi, res_x] = SOR_Xi(xi, psi, dx, dy, omega, Na, Nb, Nc, Ne, Nf, nu);
        
        resnorm_p = norm(res_p)/sqrt(dx*dy);
        resnorm_x = norm(res_x)/sqrt(dx*dy);

        figure(2)
        cla
        contour(psi);

        figure(3)
        cla
        contour(xi);
    
        

        fprintf('|res| = %g\n',resnorm_p);
        % figure(3)
        % cla
        % surface(xg,yg,v,'EdgeAlpha',0.25);
        % view(-35,45)
        % figure(2)
        % cla
        % surface(xg,yg,res,'EdgeAlpha',0.25);
        % view(-35,45)
        % pause(0.01)

        if (resnorm_p < eps)
            fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm)
            break
        end

        % if (resnorm_x < eps) && (resnorm_p < eps)
        %     fprintf('Stopped at l=%d, |res| = %g\n',l,resnorm)
        %     break
        % end
    end
    % figure(4)
    % contour(xg,yg,v)
    % [Ey,Ex] = gradient(v,dy,dx);
    % Ex = - Ex;Ey = -Ey;
    % figure(5)
    % quiver(xg,yg,Ex,Ey)
end

function [psi, xi] = boundary_conditions(psi, xi, Nx, Ny, Na, Nb, Nc, Ne, Nf, v0, dx, dy)
    
    % All things not touching the boundary
    
    % Left side
    xi(:, 1) = 0;
    psi(:, 1) = psi(:, 2);
    
    % Top
    xi(Nf, :) = 0;
    psi(Nf, :) = v0*dy + psi(Nf-1, :);

    
    % Right
    xi(:, end) = xi(:, end-1);
    psi(:, end) = psi(:, end-1);

    % For things touching the boundary

    % E
    for i = 1:1:Ne-1
        xi(1,i) = 0;
        psi(1, i) = 0;
    end

    % D
    for j = 1:1:Nb
        xi(j, Ne) = 0;
        psi(j, Ne) = psi(j, Ne-1);
    end

    % C
    for i = Ne+1:1:Ne+Nc-1
        xi(Nb, i) = xi(Nb, i-1);
        psi(Nb, i) = 0;
    end

    % B
    for j = 1:1:Nb
        xi(j, Ne+Nc) = xi(j, Ne + Nc + 1);
        psi(j, Ne + Nc) = 0;
    end

    % A
    for i = (Ne + Nc + 1):1:(Ne + Nc + Na)
        xi(1, i) = 0;
        psi(1,i) = 0;
    end
end

function [f,res_psi]  = SOR_Psi(f, source, dx, dy, omega, Na, Nb, Nc, Ne, Nf)
    vol = dx*dy; % good
    rx = dy/dx; % good
    ry = dx/dy; % good
    
    % f is psi, source is xi (or i name them that?
    res_psi = 0;

    % Psi first
    for i=2:Ne-1 
        for j=2:(Nf-1)
            % res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            % res_psi = max(res_psi, abs(res));
            % delta = omega * res / (2*(rx + ry));
            % f(i,j) = f(i,j) - delta;
            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(j,i) = f(j,i) - delta;
        end
    end

    for i = Ne:(Ne+Nc)
        for j = (Nb + 1): (Nf-1)
            % res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            % res_psi = max(res_psi, abs(res));
            % delta = omega * res / (2*(rx + ry));
            % f(i,j) = f(i,j) - delta;
            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(j,i) = f(j,i) - delta;
        end
    end
    for i = (Ne + Nc +1):(Ne + Nc + Na-1)
        for j = 2:(Nf-1)
            % res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            % res_psi = max(res_psi, abs(res));
            % delta = omega * res / (2*(rx + ry));
            % f(i,j) = f(i,j) - delta;

            res = rx*(f(j,i+1) - 2*f(j, i) + f(j, i-1)) + ry*(f(j+1,i) - 2*f(j,i) + f(j-1,i)) + vol * source(j,i);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(j,i) = f(j,i) - delta;
        end
    end
end


function [xi,res_xi]  = SOR_Xi(xi, psi, dx, dy, omega, Na, Nb, Nc, Ne, Nf, nu)
    vol = dx*dy; % good
    rx = dy/dx; % good
    ry = dx/dy; % good

    % % Define derivatives
    % dpsi_dx = (psi(i+1,j) - psi(i-1,j)) / (2*dx);
    % dpsi_dy = (psi(i,j+1) - psi(i,j-1)) / (2*dy);
    % dxi_dx  = (xi(i+1,j) - xi(i-1,j)) / (2*dx);
    % dxi_dy  = (xi(i,j+1) - xi(i,j-1)) / (2*dy);
    % 
    % % Define source function
    % source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
    % res = rx*(xi(i+1,j) - 2*xi(i,j) + xi(i-1,j)) + ry*(xi(i,j+1) - 2*xi(i,j) + xi(i,j-1)) - vol * source;
    % delta = omega * res / (2*(rx + ry));
    % xi(i,j) = xi(i,j) + delta;

    res_xi = 0;
   
    for i=2:Ne-1 
        for j=2:(Nf-1)
            dpsi_dx = (psi(i+1,j) - psi(i-1,j)) / (2*dx);
            dpsi_dy = (psi(i,j+1) - psi(i,j-1)) / (2*dy);
            dxi_dx  = (xi(i+1,j) - xi(i-1,j)) / (2*dx);
            dxi_dy  = (xi(i,j+1) - xi(i,j-1)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(i+1,j) - 2*xi(i,j) + xi(i-1,j)) + ry*(xi(i,j+1) - 2*xi(i,j) + xi(i,j-1)) - vol * source;
            res_xi = max(res_xi, abs(res));
            delta = omega * res / (2*(rx + ry));
            xi(i,j) = xi(i,j) - delta;
        end
    end

    for i = Ne:(Ne+Nc)
        for j = (Nb + 1): (Nf-1)
            % res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ...
            %       ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            % delta = omega * res / (2*(rx + ry));
            % f(i,j) = f(i,j) - delta;
            dpsi_dx = (psi(i+1,j) - psi(i-1,j)) / (2*dx);
            dpsi_dy = (psi(i,j+1) - psi(i,j-1)) / (2*dy);
            dxi_dx  = (xi(i+1,j) - xi(i-1,j)) / (2*dx);
            dxi_dy  = (xi(i,j+1) - xi(i,j-1)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(i+1,j) - 2*xi(i,j) + xi(i-1,j)) + ry*(xi(i,j+1) - 2*xi(i,j) + xi(i,j-1)) - vol * source;
            res_xi = max(res_xi, abs(res));
            delta = omega * res / (2*(rx + ry));
            xi(i,j) = xi(i,j) - delta;
        end
    end
    for i = (Ne + Nc +1):(Ne + Nc + Na-1)
        for j = 2:(Nf-1)
            % res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ...
            %       ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            % delta = omega * res / (2*(rx + ry));
            % f(i,j) = f(i,j) - delta;
            dpsi_dx = (psi(i+1,j) - psi(i-1,j)) / (2*dx);
            dpsi_dy = (psi(i,j+1) - psi(i,j-1)) / (2*dy);
            dxi_dx  = (xi(i+1,j) - xi(i-1,j)) / (2*dx);
            dxi_dy  = (xi(i,j+1) - xi(i,j-1)) / (2*dy);
        
            % Define source function
            source = (1/nu) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
            res = rx*(xi(i+1,j) - 2*xi(i,j) + xi(i-1,j)) + ry*(xi(i,j+1) - 2*xi(i,j) + xi(i,j-1)) - vol * source;
            res_xi = max(res_xi, abs(res));
            delta = omega * res / (2*(rx + ry));
            xi(i,j) = xi(i,j) - delta;
        end
    end
end
