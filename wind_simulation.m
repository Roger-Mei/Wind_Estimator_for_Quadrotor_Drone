function v_w = wind_simulation(tspan, dt, n, L_v, sigma_v)
    %% Input:
    % Wind velocity v_w = [North, East, Down]
    % Time span corresponds to tspan
    % n is the nunmber of sinusoids
    % L_v is vertical gust length scales
    % sigma_v is vertical turbulence intensities
    % wind_dir is the static dominant wind direction
    %% Output:
    % v_w = [North(x), East(y), Down(z)]
    % size(v_w) = 3xlength(tspan)
    %%  parameter initialization
    % Output
    v_w = zeros(3, length(tspan));
    % Dominant static wind velocity (Tunable)
    v_w_0 = [0.5,0.5,0.1]; 
    % L_h is the horizontal gust length scales
    L_h = 1/power((0.177 + 0.000823 * L_v), 1.2) * L_v;
    % sigma_h is the horizontal turbulence intensities
    sigma_h = 1/power((0.177 + 0.00823 * L_v), 0.4) * sigma_v;
    % index length
    l = (tspan(2)-tspan(1))/dt + 1;
    %% Build the Dryden wind gust model:
    % Set the north, east, down wind velocity
    for k = 1:3
        % Parameter Initialization:
        omega_i = [];
        psi_i = [];
        delta_omega_i = [];
        Phi_h = [];
        Phi_v = [];
        Phi = [];
        a_i = [];
        for j = 1:n
            % Generate frequency omega_i ranging from 0.1 to 1.5 rad/s
            omega_i(j) = 0.1 + rand * (1.5 - 0.1);
            % Generate phase psi_i ranging from -pi to pi
            psi_i(j) = -pi + rand * 2 * pi;
            if j == 1
                delta_omega_i(j) = omega_i(j);
            else
                delta_omega_i(j) = omega_i(j) - omega_i(j-1);
            end
            % Horizontal Power Spectral Density
            Phi_h(j) = sigma_h^2*(2*L_h/pi)*(1/(1+(L_h*omega_i(j))^2));
            % Vertical Power Spectral Density
            Phi_v(j) = sigma_v^2*(2*L_v/pi)*(1/(1+(L_v*omega_i(j))^2));
            % Combine Phi_h and Phi_v into Phi
            Phi(j) = sqrt(Phi_h(j)^2 + Phi_v(j)^2);
            % Magnitude of sinudoide
            a_i(j) = sqrt(delta_omega_i(j)*Phi(j));
        end
    
        for i = 1:l
            % time stamp
            t = dt * (i-1); 
            % Initialize summation
            sum = 0;
            % Add sinusoid with certain magnitude into the summation
            for j = 1:n
                sum  = sum + a_i(j) * sin(omega_i(j)*t + psi_i(j));
            end
            v_w(k,i) = v_w_0(k) + sum;
            if k == 3
                v_w(k,i) = 0.15 * v_w(k,i);
%                   v_w(k,i) = 0 * v_w(k,i);
            end
        end
    end
end

