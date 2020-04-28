% Kayla Ormiston
% 4/27/2020
% Solar Sail Controls - CS 595a Final Project
% This file runs the orbital and attitude dynamics for a solar sail in
% Low-Earth Orbit (LEO), using a PD conroller to control the angle (alpha)
% between the sunline and sail normal, thus controlling the amount of
% thrust provided by Solar Radiation Pressure (SRP).


% This function runs the equations to calculate the satellite dynamics
function [heuristic, Rx, Ry, Rz, Rx2, Ry2, Rz2, re] = simulation(Kp,Kd)

% Inputs:
    % Kp - controller gain driving the proportial term of the PD controller
    % Kd - controller gain driving the derivative term of the PD controller
    
% Outputs:
    % heurisitc - value based on how much the orbit has raised from the
    % previous rotation around the Earth
    % Rx - Array of the X-coordinate position of the satellite using the predicted
    % orbit without SRP effects
    % Ry - Array of the Y-coordinate position of the satellite using the predicted
    % orbit without SRP effects
    % Rz - Array of the Z-coordinate position of the satellite using the predicted
    % orbit without SRP effects
    % R2x - Array of the X-coordinate position of the satellite using
    % calculated orbit with SRP effects
    % R2y - Array of the Y-coordinate position of the satellite using
    % calculated orbit with SRP effects
    % R2z - Array of the Z-coordinate position of the satellite using
    % calculated orbit with SRP effects
    % re - Radius of the Earth
    
    
    % Define constants
    altitude = 500000;      % meters
    re= 6371000;            % meters
    G = 6.67430*10^(-11);   % Nm2/kg2
    M_earth = 5.972*10^24;  % kg
    mu = G*M_earth;         % Earth's gravitational constant        

    % Initialization of attitude angle alpha
    alpha = 0;              % rad - angle between sun and normal to sail
    alpha_dot = 0.1;        % rad/s - initial rate of change of alpha

    % Parameters needed to calculate angle of the sun's rays
    n = juliandate(datetime)-2451545;               % number of days since J2000
    L = (280.459+0.98564736*n)*pi/180;              % rad - mean longitude of the sun
    L = mod(L, 2*pi);
    M = (357.529+0.98560023*n)*pi/180;              % rad - mean anomaly of the sun
    M = mod(M, 2*pi);
    gamma = (L+1.915*sin(M)+0.02*sin(2*M))*pi/180;  % solar ecliptic longitude
    gamma = mod(gamma, 2*pi);
    epsilon = (23.439-3.56*10^7*n)*pi/180;          % obliquity of the sun

    r = re + altitude;          % meters - distance from center of earth to satellite
    a = r;                      % meters - semimajor axis of the orbit
    v = sqrt(mu*(2/r-1/a));     % m/s - magnitude of the velocity of the satellite
    R = r * [1,0,0];            % meters - initial position vector of the satellite
    V = v * [0,0,1];            % m/s - initial velocity vector of the satellite
    H = cross(R, V);            % initial angular momentum vector of the satellite
    h = sqrt(dot(H, H));        % initial magnitude of the angular momentum
    e = sqrt(-(h^2/(a*mu))+1);  % initial eccentricity of the orbit
    theta = 0;                  % initial true anomaly of the orbit
    Me = theta;                 % initial mean anomaly of the orbit
    E = Me + e*sin(Me) + ...    % initial eccentric anomaly
        e^2/2*sin(2*Me) + ...
        e^3/8*(3*sin(3*Me) - sin(Me));
    T_initial = 2*pi*r^(3/2)/sqrt(mu); % initial period of the orbit
    T = T_initial;
    
    % Initialize calculated orbit with SRP
    r2(1) = a*(1-e*cos(E));     % magnitude of the radius
    Rx2(1) = r2(1)*cos(theta);  % X-component of position
    Ry2(1) = 0;                 % Y-component of position
    Rz2(1) = r2(1)*sin(theta);  % Z-component of position
    
    dt = 0.1;                   % define the time step of the simulation

    % initialize counter and flag variables
    i = 2;
    valid_run = true;
    % Loop for a full revolution of an orbit
    while(abs(theta) <= 2*pi)
        % Get the values from the next time step
        [alpha, R, V, Me, e, T, theta, r, alpha_dot] = next_step([Rx2(i-1), Ry2(i-1), Rz2(i-1)], ...
                                                        V, Me, e, T, dt, mu, gamma, epsilon, ...
                                                        alpha, alpha_dot, Kp, Kd);

        % Store the new position of the satellite
        Rx2(i) = R(1);
        Ry2(i) = R(2);
        Rz2(i) = R(3);
        r2(i) = r;

        % Store the new attitude of the satellite
        alpha2(i) = alpha;
        i = i+1;

        % If any of these variables gain imaginary components, stop the
        % simulation
        if (isreal(theta) == false || isreal(alpha) == false || ...
            isreal(R(1)) == false || isreal(R(3)) == false || ...
            isreal(V(1)) == false || isreal(V(3)) == false || isreal(Me) == false || ...
            isreal(e) == false || isreal(T) == false || isreal(r) == false || i > 100000)
            
            valid_run = false;
            break
        end
    end

    % Evaluate whether the simulation was run successfully and return the
    % heuristic
    if (valid_run)
        heuristic = r2(end)-r2(1);
    else
        heuristic = 0;
    end

    % Calculate projected orbit without SRP
    t = 0:10000;
    th = (2*pi)/T_initial*t;
    Rx = r * cos(th);
    Ry = zeros(1,10001);
    Rz = r * sin(th);    

    % This function calculates each iteration of the simulation
    function [alpha, R, V, e, Me, T, theta, r, alpha_dot] = next_step(R, V, e,...
                                                           Me, T, dt, mu, gamma,...
                                                           epsilon, alpha, alpha_dot, Kp, Kd)
    % Inputs:
        % R - Position vector of the satellite
        % V - velocity vector of the satellite
        % e - eccentricity of the orbit
        % Me - Mean anomaly of the orbit
        % T - Period of the Orbit
        % dt - time step of the simulation
        % mu - Earth's gravitational constant
        % gamma - solar ecliptic longitude
        % epsilon - obliquity of the sun
        % alpha - angle between sun and normal to sail
        % alpha_dot - angular velocity of alpha
        % Kp - controller gain driving the proportial term of the PD controller
        % Kd - controller gain driving the derivative term of the PD controller
        
    % Outputs:
        % alpha - angle between sun and normal to sail
        % R - Position vector of the satellite
        % V - velocity vector of the satellite
        % e - eccentricity of the orbit
        % Me - Mean anomaly of the orbit
        % T - Period of the Orbit
        % theta - true anomaly of the orbit
        % r - magnitude of the position vector
        % alpha_dot - angular velocity of alpha
        
        
        % Define constants needed for calculation of solar radiation
        % pressure
        ur = cos(gamma);                % transformation parameter from RSW frame to body frame
        us = cos(epsilon)*sin(gamma);   % transformation parameter from RSW frame to body frame
        c = 2.998*10^8;                 % m/s - speed of light
        m = 10;                         % kg - mass of satelite
        Cr = 2;                         % radiation pressure coeff (2=reflect 1=absorb)
        Rsun = 149.6*10^9;              % m - dist from center of the sun
        area = 100;                     % m^2 - Area of the sail
        So = 63.15*10^6;                % W/m^2 - radiated power intensity
        Ro = 696000000;                 % m - radius of the photosphere
        S = So*(Ro/(Rsun))^2;           % electromagnetic radiation
        N = 15000;                      % SRP factor
        Psr = (S*area*cos(alpha)*Cr)/(m*c)*N;  % Solar Radiation Pressure
        
        % Calculate the next step in the mean anomaly
        Me = mod(Me + 2*pi/T * dt, 2*pi);
        
        % Calculate the next step in the eccentric anomaly
        % This is an iterative process; the first guess for E is
        % initialized here
        if (Me < pi)
            E = Me + e/2;
        else
            E = Me - e/2;
        end
        precision = 1;
        while (precision > 0.0001)
            E_old = E;
            E = E_old - (E_old - e*sin(E_old)-Me)/(1-e*cos(E_old));
            precision = abs(E-E_old);
        end
        
        % Calculate the next step in the true anomaly
        if (Me < pi)
            theta = acos((e-cos(E))/(e*cos(E)-1));
        else
            theta = -acos((e-cos(E))/(e*cos(E)-1));
        end
        
        % Calculate the new angular momentum
        H = cross(R, V);
        h = sqrt(dot(H, H));
        
        % Calculate the new magnitude of the position
        r = sqrt(dot(R,R));
        
        % Calculate the change in theta due to SRP
        theta_srp = h/r^2-(Psr/(e*h))*(h^2/mu*cos(theta)*ur-(r+(h^2/mu)*sin(theta)*us));
        if (theta_srp == inf)
            theta_srp=0;
        end
        theta = theta + theta_srp*dt;
        
        % Calculate the change in angular momentum due to SRP
        h_srp = -Psr*r*us;
        h = h + h_srp*dt;
        
        % Calculate the change in eccentricity due to SRP
        e_srp = -Psr*(h/mu*sin(theta)*ur+1/(mu*h)*((h^2+mu*r)*cos(theta)+mu*e*r)*us);
        
        % Define the principle moments of inertia for the satellite
        J1 = 6000;              %kg*m^2
        J2 = 3000;              %kg*m^2
        J3 = 3000;              %kg*m^2
       
        % Define the ideal (commanded) alpha for the controller
        alphaC = 0;
        
        % Define the control torque
        u2 = -Kp*(alpha-alphaC)-Kd*alpha_dot;
        
        % Calculate the acceleration of alpha
        alpha_ddot = (u2+3*mu/r^3*(J1-J3)*sin(alpha-theta)*cos(alpha-theta))/J2;
        
        % Obtain new alpha
        alpha = alpha + alpha_dot*dt + 1/2*alpha_ddot*dt^2;
        alpha_dot = alpha_dot + alpha_ddot*dt;
        
        % Get new R using new theta
        r = h^2/(mu*(1+e*cos(theta)));
        R = [r*cos(theta);
            0;
            r*sin(theta)];
        
        % Get new V
        p = r*(1+e*cos(theta));
        V = [-sqrt(mu/p)*sin(theta);
            0;
            sqrt(mu/p)*(e+cos(theta))];
        v = sqrt(dot(V,V));
        
        % Get new H
        H = cross(R,V);
        h = sqrt(dot(H,H));
        
        % Get new e
        e = (E-Me)/sin(E)+e_srp*dt;
        
        % Get new T
        T = 2*pi/mu^2*(h/sqrt(1-e^2))^3;
    end
end
