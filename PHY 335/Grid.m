% The class Grid defines a computational grid along with a set of
% methods used to solve some problem on that grid.
%
% The constructor takes either 1 or 6 arguements.
% 6 arguements:  grid = Grid(Xmin,Xmax,Ymin,Ymax,Nx,Ny)
%                These define the Cartesian computational domain
%                and the discretization
% 1 arguement: grid2 = Grid(grid1)
%              Given an instance of Grid called grid1, create and
%              new grid with the discretization length in each
%              direction cut in half.  (A 'linked-list' of grids
%              is created.)
classdef Grid < handle
    properties
        level % Level of the current grid in the linked list
        cGrid % The next coarser grid in the linked list
        fGrid % The next finer grid in the linked list
        v % The current estimated solution
        vguess % An initial guess solution saved for future use
        rho % The source function
        nx % The number of discrete points in the x direction
        ny % The number of discrete points in the y direction
        dx % The discretization length in the x direction
        dy % The discretization length in the y direction
        x % 1D vector of x grid locations
        y % 1D vector of y grid locations
        xg % 2D array of x grid locations (meshgrid)
        yg % 2D array of y grid locations (meshgrid)
    end
    methods
        function obj = Grid(varargin)
            if nargin == 6
                % varargin{1}=Xmin, varargin{2}=Xmax
                % varargin{3}=Ymin, varargin{4}=Ymax
                % varargin{5} = Nxc, varargin{6}=Nyc
                obj.level = 1;
                obj.nx = varargin{5};
                obj.ny = varargin{6};
                obj.dx = (varargin{2}-varargin{1})/(varargin{5}-1);
                obj.dy = (varargin{4}-varargin{3})/(varargin{6}-1);
                obj.x = varargin{1}:obj.dx:varargin{2};
                obj.y = varargin{3}:obj.dy:varargin{4};
                obj.cGrid = 0;
                obj.fGrid = 0;
                obj.vguess = zeros(varargin{5},varargin{6});
            else
                % varargin{1} = Parent
                obj.level = varargin{1}.level+1;
                obj.nx = 2*varargin{1}.nx-1;
                obj.ny = 2*varargin{1}.ny-1;
                obj.dx = varargin{1}.dx/2;
                obj.dy = varargin{1}.dy/2;
                obj.x = varargin{1}.x(1):obj.dx:varargin{1}.x(end);
                obj.y = varargin{1}.y(1):obj.dy:varargin{1}.y(end);
                obj.cGrid = varargin{1};
                obj.fGrid = 0;
                varargin{1}.fGrid = obj;
                obj.vguess = obj.cGrid.prolong2(obj.cGrid.v);
            end
            [obj.yg,obj.xg] = meshgrid(obj.y,obj.x);
            obj.DirichletBC();
            obj.v = obj.vguess;
            obj.rho = obj.source();
        end
        % Compute the result of the discrete elliptic operator
        % (representing the 2D Laplacian) acting on the discrete
        % solution 'u'
        function Lu = Elliptic(obj,u)
            rx = obj.dy/obj.dx;
            ry = obj.dx/obj.dy;
            Lu = - 2.0*(rx+ry)*u;
            Lu(:,1) = 0.0; Lu(:,end)=0.0;Lu(1,:)=0.0;Lu(end,:)=0.0;
            for j=2:obj.ny-1
                for i=2:obj.nx-1
                    Lu(i,j) = Lu(i,j) + rx*(u(i+1,j)+u(i-1,j)) + ry*(u(i,j+1)+u(i,j-1));
                end
            end
        end
        % Compute the residual based on the current stored solution
        % 'obj.v'. 
        function res = residual(obj)
            vol = obj.dx*obj.dy;
            res = obj.Elliptic(obj.v);
            res(2:end-1,2:end-1)=res(2:end-1,2:end-1)...
                - vol*obj.rho(2:end-1,2:end-1);
        end
        % Compute the source function for the elliptic problem 
        % being solved.
        function rho = source(obj)
            rho = zeros(obj.nx,obj.ny);
            for j=1:obj.ny
                for i=1:obj.nx
                    rho(i,j) = -300.0*exp(-20.0*((obj.x(i)-0.75)^2 + (obj.y(j)-0.5)^2))*sin(4.0*pi*obj.x(i)/3);
                end
            end
        end
        % Set the Dirichlet boundary conditions on the currently
        % stored guess solution 'obj.vguess'
        function DirichletBC(obj)
            obj.vguess(:,1) = 4.0*obj.x/3.0;
            obj.vguess(:,end) = 2.0*(obj.x - 1.0);
            obj.vguess(1,:) = -2.0*obj.y;
            obj.vguess(end,:) = 2.0 - obj.y;
        end
        % Compute a Jacobi update to the stored solution 'obj.v'
        % using the relaxation parameter 'omega'.  
        % 'omega'=1 corresponds to standar Jacobi relaxation
        function  Jacobi(obj,omega)
            res = obj.residual();
            rx = obj.dy/obj.dx;
            ry = obj.dx/obj.dy;
            dv = res/(-2.0*(rx+ry));
            obj.v = obj.v - omega*dv;
        end
        % Perform Gauss-Seidel relaxation on the currently stored
        % solution 'obj.v' using the relaxation parameter 'omega'.
        % 1<'omega'<2 gives successive-over-relaxation
        function SOR(obj,omega)
            vol = obj.dx*obj.dy;
            rx = obj.dy/obj.dx;
            ry = obj.dx/obj.dy;
            lnx = obj.nx;
            lny = obj.ny;
            corfac = omega/(-2*(rx+ry));
            for i=2:lnx-1
               for j=2:lny-1
                   delta = rx*(obj.v(i+1,j)-2.0*obj.v(i,j) + obj.v(i-1,j))...
                       + ry*(obj.v(i,j+1)-2.0*obj.v(i,j) + obj.v(i,j-1))...
                       - vol*obj.rho(i,j);
                   obj.v(i,j) = obj.v(i,j) - corfac*delta;
               end
            end
        end
        % Perform Red-Black SOR relaxation on the currently stored
        % solution 'obj.v' using the relaxation parameter 'omega'.
        % 1<'omega'<2 gives successive-over-relaxation
        function RedBlackSOR(obj,omega)
            vol = obj.dx*obj.dy;
            rx = obj.dy/obj.dx;
            ry = obj.dx/obj.dy;
            corfac = omega/(-2.0*(rx+ry));
            % Red squares
            vc = obj.v;
            delta = -2.0*(rx+ry)*vc(2:2:end-1,2:2:end-1)...
                + rx*(vc(3:2:end,2:2:end-1)...
                      + vc(1:2:end-2,2:2:end-1))...
                + ry*(vc(2:2:end-1,3:2:end)...
                      + vc(2:2:end-1,1:2:end-2))...
                - vol*obj.rho(2:2:end-1,2:2:end-1);
            obj.v(2:2:end-1,2:2:end-1) = ...
                vc(2:2:end-1,2:2:end-1) - corfac*delta;
            delta = -2.0*(rx+ry)*vc(3:2:end-1,3:2:end-1)...
                + rx*(vc(4:2:end,3:2:end-1)...
                      + vc(2:2:end-2,3:2:end-1))...
                + ry*(vc(3:2:end-1,4:2:end)...
                      + vc(3:2:end-1,2:2:end-2))...
                - vol*obj.rho(3:2:end-1,3:2:end-1);
            obj.v(3:2:end-1,3:2:end-1) = ...
                vc(3:2:end-1,3:2:end-1) - corfac*delta;
            % Black squares
            vc = obj.v;
            delta = -2.0*(rx+ry)*vc(3:2:end-1,2:2:end-1)...
                + rx*(vc(4:2:end,2:2:end-1)...
                      + vc(2:2:end-2,2:2:end-1))...
                + ry*(vc(3:2:end-1,3:2:end)...
                      + vc(3:2:end-1,1:2:end-2))...
                - vol*obj.rho(3:2:end-1,2:2:end-1);
            obj.v(3:2:end-1,2:2:end-1) = ...
                vc(3:2:end-1,2:2:end-1) - corfac*delta;
            delta = -2.0*(rx+ry)*vc(2:2:end-1,3:2:end-1)...
                + rx*(vc(3:2:end,3:2:end-1)...
                      + vc(1:2:end-2,3:2:end-1))...
                + ry*(vc(2:2:end-1,4:2:end)...
                      + vc(2:2:end-1,2:2:end-2))...
                - vol*obj.rho(2:2:end-1,3:2:end-1);
            obj.v(2:2:end-1,3:2:end-1) = ...
                vc(2:2:end-1,3:2:end-1) - corfac*delta;
        end
        % Use the sine transformation to update the current guess
        % solution stored in 'obj.v'
        function SinSolver(obj)
            lnx = obj.nx;
            lny = obj.ny;
            rx = obj.dy/obj.dx;
            ry = obj.dx/obj.dy;
            % compute the residual
            res = obj.residual();
            % Compute the 2D sine transform of res
            reshat = dst(res(2:lnx-1,2:lny-1)');
            reshat = dst(reshat');
            % Solve the PDE
            vhat = reshat;
            for i=1:lnx-2
                tmp = rx*(cos(pi*i/(lnx-1)) - 1.0);
                for j=1:lny-2
                    vhat(i,j) = vhat(i,j)/(2.0*(tmp + ry*(cos(pi*j/(lny-1))-1.0)));
                end
            end
            vhat = idst(vhat');
            obj.v(2:lnx-1,2:lny-1) = obj.v(2:lnx-1,2:lny-1) - idst(vhat');
        end
        % Plot the currently stored solution 'obj.v'
        function plotv(obj)
            surface(obj.xg,obj.yg,obj.v,'EdgeAlpha',0.15);
        end
    end
end