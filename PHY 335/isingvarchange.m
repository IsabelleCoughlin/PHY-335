% Varying the lattice size
function isingvarchange(m,seed)
digit=0;
    % Set the Temperature values to compute
while m<=16
    digit=digit+1;
    Tstart = 1;
    Tend = 5;
    Tstep = 0.2;
    Temps = Tstart:Tstep:Tend;
    % Storage for various expectation values at each temperature
    % M = zeros(length(Temps),1);
    % MM = zeros(length(Temps),1);
    % E = zeros(length(Temps),1);
    % EE = zeros(length(Temps),1);
    % Cv = zeros(length(Temps),1);
    % Ms = zeros(length(Temps),1);
    % Loop over all temperatures
    for k=1:length(Temps)
        % Create the Class for the Lattice
        if k==1 
            lat = lattice(Temps(k),m,seed);
        end
        % Call the routine to carry out the relaxation and 
        % measurement sweeps
        Tising(lat,Temps(k),m);
        % Gather the data and save into the arrays.
        [M(k,digit),MM(k,digit),E(k,digit),EE(k,digit)] = lat.CollectData();
      	Ms(k,digit) = m^2*(MM(k,digit) -(M(k,digit)*M(k,digit)))/Temps(k);
        Cv(k,digit) = m^2*(EE(k,digit)-(E(k,digit)*E(k,digit)))/(Temps(k)^2);
        % Display the information
        fprintf(' T= %f: |M|= %f, E = %f: Ms = %f, C_v = %f \n',Temps(k),abs(M(k)),E(k),Ms(k),Cv(k) );
        lat.image();
        title(sprintf('T = %f\n',Temps(k)));
        lat.resetstats(Tstep); % using this to keep the same lattice but to reset
        %the energy and magnatic parameters see more in lattice->reset
        %stats
        %pause(0.05);
        
    end
    m=m*2;
end
    % Plot the various physical quantities vs temperature
    % Energy
    figure(2);
    plot(Temps,E);
    xlabel('T');
    ylabel('\langle{E}\rangle');
    title('Energy per spin');
    legend("4","8","16")
    % Magnetization
    figure(3);
    plot(Temps,abs(M));
    xlabel('T');
    ylabel('|\langle{M}\rangle{}|');
    title('Magnetization per spin');
   legend("4","8","16")
    % Specific heat
    figure(4);
    plot(Temps,Cv);
    xlabel('T');
    ylabel('C_v');
    title('Specific Heat');
    legend("4","8","16")
    % Magnetic Susceptibility
    figure(5);
    plot(Temps,Ms);
    xlabel('T');
    ylabel('\chi');
    title('Magnetic Susceptibility');
    legend("4","8","16");
end


function Tising(lat,T,m)
    figure(1);
    hold on
    lat.image;
    title(sprintf('T = %f\n',T));
    % Sweep to allow the system to come to equilibrium
    fprintf('Thermalize----\n');
    for k=1:1000
        lat.sweep;
        if mod(k,10)==0
            lat.image;
            fprintf('%5i: |M|= %f, E = %f\n',k,abs(lat.Magnetization),lat.Energy);
        end
    end
    figure(1);
    lat.image();
    % Sweep to collect the data
    fprintf('Collect Data\n');
    for k=1:1000
        lat.sweep;
        if mod(k,10)==0
            lat.stat();
        end
        if mod(k,50)==0
            figure(1);
            lat.image;
            fprintf('%5i: |M|= %f, E = %f\n',k,abs(lat.Magnetization),lat.Energy);
        end
    end
    hold off
end
