classdef lattice < handle
    % Class to store a 2-D lattice of "spins" and methods to perform
    % lattice update "sweeps", plot the lattice, and accumulate statistics.
    % 
    properties
        T    % Lattice temperature
        data % The 2-D lattice data (stores +1 or -1 in each location)
        ran  % The random number generator to use
        m    % The size of the lattice (m rows and m columns)
    end
    methods
        function obj = lattice(T,m,seed)
            % Constructor
            obj.T = T;
            obj.m = m;
            obj.data = zeros(m,m);
            obj.ran = NumericalRecipes.Ran(seed);
        end
        %
        function M = Magnetization(obj)
            % Compute and return the magnetization/spin of the current
            % lattice 
         end
        %
        function E = Energy(obj)
            % Compute and return the energy/spin of the current
            % lattice 
        end
        %
        function e = energy(obj)
            % Compute and return the total energy of the current lattice
        end
        %
        function sweep(obj)
            % Apply the Metropolis algorithm to m*m randomly chosen spins.
        end
        %
        function de = DeltaEnergy(obj,i,j)
            % Compute the change in the energy if the spin at lattice site
            % (i,j) flips
        end
        %
        function stat(obj)
            % Collect statistics for the current state of the lattice
        end
        %
        function resetstats(obj)
            % Reset the accumulated statistics to zero.
        end         
        %
        function [M,MM,E,EE] = CollectData(obj)
            % return the expectation values <M>, <M^2>, <E>, and <E^2>
        end
        %
        function image(obj)
            % Produce a visualization of the lattice
            img = uint8(floor(0.5*(obj.data+1.0)));
            colormap([0 0 1;1 0 0]);
            image(img);
            pause(0.05);
        end
    end
end