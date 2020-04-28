% Kayla Ormiston
% 4/27/2020
% Solar Sail Controls - CS 595a Final Project
% This file runs a simulated annealing algorithm to determine the optimal
% PD controller gains for on a simulation of a solar sail spacecraft  

clear
clc

inital_probability = 0.6;
final_probability = 0.1;

interations = 200;

T0 = -1/log(inital_probability);
TN = -1/log(final_probability);
Tstep = (TN/T0)^(1/interations);

i = 2;
T = T0;
while T>TN
    current = 1+ (10000-1).*rand(1,2);
    Kp = current(1);
    Kd = current(2);
    
    % Get the heuristic for a new simulation
    heuristic(i) = simulation(Kp, Kd);
    
    % Set goal for heuristic success
    if heuristic(i) >= 700000
        break
    end

    %Select random successor
    successors(1) = 1+ (10000-1).*rand(1,1);
    successors(2) = 1+ (10000-1).*rand(1,1);
    next(1) = successors(1);
    next(2) = successors(2);    
 
    %Calculate deltaE - difference in estimated cost of new state vs. current state
    deltaE = heuristic(i) - heuristic(i-1);

    %If deltaE is zero, i.e. improved move, then set current state to next.
    if deltaE < 0
        current = next;
            
    %Otherwise, if a random number is less than the probability value for moving to a 
    %worse state, then move to worse state
    elseif randn < exp(-(deltaE/T))
        current = next;
    end

    % Update temperature
    T = T * Tstep;
    i = i+1;
end

% Print the results
if heuristic(end) < 700000
    fprintf('Solution not found, please run again\n')
else
    fprintf('Heuristic: %d, Kp: %d, Kd: %d, i: %f\n', heuristic(end), current(1), current(2), i)

    % Get the values of the best heuristic result
    [~, Rx, Ry, Rz, Rx2, Ry2, Rz2, re] = simulation(current(1), current(2)); 

    % Plot the heuristic evolution
    f1 = figure();
    non_zero_heuristic = heuristic(heuristic>0);
    plot(1:length(non_zero_heuristic),non_zero_heuristic)


    % Plot the best resulting orbit
    f2 = figure();
    for i = 2:200:length(Rx2)
        % Plot the base orbit without any disturbances
        plot3(Rx, Ry, Rz);
        hold on

        % Plot the new orbit
        scatter3(Rx2(i), Ry2(i), Rz2(i), 100)

        % Resize the axis
        axis([-re re -re re -re re ]*1.5)
        pause(0.1)
        hold off
    end
end

