Nx = 10;
Ny = 10;
Nc = 1; Nb = 3; Na = 5; Ne = 4; Nf = 10;
dx = (Ne + Nc + Na)/(Nx-1);
dy = (Nf)/(Ny-1);
psi = zeros(Ny, Nx);
xi = zeros(Ny, Nx);
omega = 1.0;
[psi, res_p] = SOR_Psi(psi, xi, dx, dy, omega, Na, Nb, Nc, Ne, Nf);


function [f,res_psi]  = SOR_Psi(f, source, dx, dy, omega, Na, Nb, Nc, Ne, Nf)
    vol = dx*dy; % good
    rx = dy/dx; % good
    ry = dx/dy; % good
    res_psi = 0;
    
    % f is psi, source is xi (or i name them that?

    % Psi first
    for i=2:Ne-1 
        for j=2:(Nf-1)
            res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(i,j) = f(i,j) - delta;
        end
    end

    for i = Ne:(Ne+Nc)
        for j = (Nb + 1): (Nf-1)
            res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(i,j) = f(i,j) - delta;
        end
    end
    for i = (Ne + Nc +1):(Ne + Nc + Na-1)
        for j = 2:(Nf-1)
            res = rx*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) + ry*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) + vol * source(i,j);
            res_psi = max(res_psi, abs(res));
            delta = omega * res / (2*(rx + ry));
            f(i,j) = f(i,j) - delta;
        end
    end
end
