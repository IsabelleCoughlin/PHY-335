classdef lattice < handle
    % Class to store a 2-D lattice of "spins" and methods to perform
    % lattice update "sweeps", plot the lattice, and accumulate statistics.
    
    properties
        T    % Lattice temperature
        data % The 2-D lattice data (stores +1 or -1 in each location)
        ran  % The random number generator to use
        m    % The size of the lattice (m rows and m columns)
        Mag  % Sum of magnetization values
        M2   % Sum of magnetization squared
        E    % Sum of energy values
        E2   % Sum of energy squared
        samples % Number of samples collected
    end
    methods
        function obj = lattice(T,m,seed)
            % Constructor
            obj.T = T;
            obj.m = m;
            %obj.data = (randi([0, 1], [m,m])*2)-1;
            obj.data = ones([m,m]);
            obj.ran = NumericalRecipes.Ran(seed);
            obj.Mag = 0;      
            obj.M2 = 0;     
            obj.E = 0;      
            obj.E2 = 0;     
            obj.samples = 0;    
        end
        %
        function M = Magnetization(obj)
            % Compute and return the magnetization/spin of the current
            M = 0;
            for i = 1:obj.m
                for j = 1:obj.m
                    M = M  + obj.data(i,j);
                end
            end
            M = abs(M/obj.m^2); % Divided by lattice size

         end
        
        function E = Energy(obj)
            % Compute and return the energy/spin of the current lattice
            %M = Magnetization(obj);
            e = energy(obj);
            E = e/(obj.m^2);
        end
        %
        function e = energy(obj)
            % Compute and return the total energy of the current lattice
            % Current data object
            S2 = 0; % Down
            S3 = 0; % Right
            e = 0;

            % Very redundant bounds, but logically correct
            for i = 1:obj.m
                for j = 1:obj.m
                    S1 = obj.data(i,j);
                    if i == 1
                        if j == 1
                            
                            S2 = obj.data(i+1, j);
                            S3 = obj.data(i, j+1);
                            
                        elseif j == obj.m
                            S2 = obj.data(i+1, j);
                            S3 = obj.data(i, 1);
                        else
                            S2 = obj.data(i+1, j);
                            S3 = obj.data(i, j+1);
                            
                        end
                    elseif i == obj.m
                        if j == 1
                            S2 = obj.data(1, j);
                            S3 = obj.data(i, j+1);
                            
                        elseif j == obj.m
                           S2 = obj.data(1, j);
                           S3 = obj.data(i, 1);
                            
                        else
                            S2 = obj.data(1, j);
                            S3 = obj.data(i, j+1);
                            
                        end
                    elseif j == 1
                        
                        S2 = obj.data(i+1, j);
                        S3 = obj.data(i, j+1);
                        
                    elseif j == obj.m
                        S2 = obj.data(i+1, j);
                        S3 = obj.data(i, 1);
                    else
                        %Normal case
                        S2 = obj.data(i+1, j);
                        S3 = obj.data(i, j+1);
                    end
                        
                    e = e + S1*(S2 + S3);
                end
            end
            e = -1*e;
        end
        %
        function sweep(obj)
            
            % Apply the Metropolis algorithm to m*m randomly chosen spins.
            % Randomly picking a point in the lattice and picking it NxN
            % times. 

            %Pick a random value between 0 and m (ceiling)
            sweeps = 0;
            for f = 1:obj.m*obj.m
                % Call stat every ten sweeps
                if mod(sweeps, 10) == 0
                    stat(obj)
                end
                sweeps = sweeps + 1;
                xi = ceil(obj.ran.doub*obj.m);
                yi = ceil(obj.ran.doub*obj.m);
               
                % Calculate dE_s
                
                de = DeltaEnergy(obj, xi, yi);
                
                if exp(-de/obj.T) > obj.ran.doub
                    obj.data(xi,yi) = -1*obj.data(xi, yi);
                end
            end
        end
        
        function de = DeltaEnergy(obj,i,j)
            og_spin = obj.data(i,j);
            S1 = 0;
            S2 = 0;
            S3 = 0;
            S4 = 0;
            if i == 1
                if j == 1
                    S1 = obj.data(i, obj.m);
                    S2 = obj.data(i+1, j);
                    S3 = obj.data(i, j+1);
                    S4 = obj.data(obj.m, j);
                    
                elseif j == obj.m
                    S1 = obj.data(i, j-1);
                    S2 = obj.data(obj.m, j);
                    S3 = obj.data(i, 1);
                    S4 = obj.data(i+1, j);
                else
                    S1 = obj.data(i, j-1);
                    S2 = obj.data(obj.m, j);
                    S3 = obj.data(i, j+1);
                    S4 = obj.data(i+1, j);
                end
            elseif i == obj.m
                if j == 1
                    S1 = obj.data(i, obj.m);
                    S2 = obj.data(i-1, j);
                    S3 = obj.data(i, j+1);
                    S4 = obj.data(1, j);
                elseif j == obj.m
                    S1 = obj.data(i, j-1);
                    S2 = obj.data(i-1, j);
                    S3 = obj.data(i, 1);
                    S4 = obj.data(1, j);
                else
                    S1 = obj.data(i, j-1);
                    S2 = obj.data(i-1, j);
                    S3 = obj.data(i, j+1);
                    S4 = obj.data(1, j);
                end
            elseif j == 1
                S1 = obj.data(i, obj.m);
                S2 = obj.data(i-1, j);
                S3 = obj.data(i, j+1);
                S4 = obj.data(i+1, j);
            elseif j == obj.m
                S1 = obj.data(i, j-1);
                S2 = obj.data(i-1, j);
                S3 = obj.data(i, 1);
                S4 = obj.data(i+1, j);
            else
                %Normal case
                S1 = obj.data(i, j-1);
                S2 = obj.data(i-1, j);
                S3 = obj.data(i, j+1);
                S4 = obj.data(i+1, j);
            end
            % Compute the change in the energy if the spin at lattice site
            % (i,j) flips
            de = 2*og_spin*(S1+S2+S3+S4);
        end

        function stat(obj)
            %Sampling the system every N_s (10) sweeps
            
            % Collect statistics for the current state of the lattice
            M = Magnetization(obj);
            E = Energy(obj);

            obj.Mag = obj.Mag + M;
            obj.M2 = obj.M2 + M*M;
            obj.E = obj.E + E;
            obj.E2 = obj.E2 + E*E;
            obj.samples = obj.samples + 1;
        end
        %
        function resetstats(obj)
            % Reset the accumulated statistics to zero.
            obj.Mag = 0;
            obj.M2 = 0;
            obj.E = 0;
            obj.E2 = 0;
            obj.samples = 0;
        end         
        %
        function [M,MM,E,EE] = CollectData(obj)
            % return the expectation values <M>, <M^2>, <E>, and <E^2>
           
            M = obj.Mag/obj.samples;
            MM = obj.M2/obj.samples;
            E = obj.E/obj.samples;
            EE = obj.E2/obj.samples;
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
