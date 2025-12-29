clear; close all; clc;

%% Excavation Dewatering Calculator
% This MATLAB script calculates groundwater discharge (dewatering/pumping quantity)
% required during excavation construction using multiple methods:
% - Thiem's formula (confined and unconfined aquifers)
% - Dupuit's assumptions for unconfined flow
% - Analytical solution for circular excavation
% - Wells method for rectangular excavation

%% INPUT PARAMETERS
% Aquifer properties
aquifer.hydraulic_conductivity = 1e-4; % m/s (typical for sand)
aquifer.thickness = 10; % m (aquifer thickness)
aquifer.porosity = 0.3; % dimensionless
aquifer.specific_yield = 0.2; % dimensionless (for unconfined)

% Excavation geometry
exca.shape = 'rectangular'; % 'circular' or 'rectangular'
exca.width = 20; % m (for rectangular)
exca.length = 30; % m (for rectangular)
exca.radius = 10; % m (for circular)
exca.depth = 5; % m (excavation depth below water table)

% Groundwater conditions
gw.initial_water_table = 2; % m (depth from ground surface)
gw.drawdown_required = 3; % m (required drawdown at excavation)
gw.well_radius = 0.15; % m (typical well radius)

% Pumping system
pump.well_diameter = 0.3; % m
pump.well_efficiency = 0.75; % dimensionless
pump.max_capacity = 100; % m³/h (per well)

% Simulation parameters
sim.safety_factor = 1.2; % factor of safety
sim.num_wells = 4; % initial estimate for number of wells
sim.time_phases = [1, 3, 7]; % days (construction phases)

%% UNIT CONVERSIONS
mm = 1e-3; % mm to m
m3_per_s_to_m3_per_h = 3600; % conversion factor
m3_per_s_to_m3_per_day = 86400; % conversion factor

%% CALCULATIONS

% 1. Thiem's Formula (Confined Aquifer)
% Q = (2 * π * K * b * (h2 - h1)) / ln(r2/r1)
% For dewatering, we calculate required pumping rate to achieve drawdown

% Equivalent radius for rectangular excavation
if strcmp(exca.shape, 'rectangular')
    % Using the method of equivalent radius
    exca.equivalent_radius = sqrt((exca.width * exca.length) / pi);
else
    exca.equivalent_radius = exca.radius;
end

% Radius of influence (approximate)
R = 300 * sqrt(aquifer.hydraulic_conductivity * aquifer.thickness * sim.safety_factor);

% Drawdown at well (s_w) and at excavation perimeter (s_e)
s_w = gw.drawdown_required * sim.safety_factor;
s_e = gw.drawdown_required;

% Thiem's formula for confined aquifer
Q_thiem_confined = (2 * pi * aquifer.hydraulic_conductivity * aquifer.thickness * (s_w - s_e)) / ...
                   log(R / pump.well_radius);

% 2. Thiem's Formula (Unconfined Aquifer)
% Q = (π * K * (h2² - h1²)) / ln(r2/r1)
% For unconfined, we use the Dupuit approximation

% Initial and final water table elevations
H_initial = aquifer.thickness + gw.initial_water_table;
H_final = H_initial - s_w;

Q_thiem_unconfined = (pi * aquifer.hydraulic_conductivity * (H_initial^2 - H_final^2)) / ...
                     log(R / pump.well_radius);

% 3. Dupuit's Assumptions for Unconfined Flow
% Q = (π * K * (H² - h²)) / ln(R/r_w)
% This is essentially the same as Thiem's unconfined formula
Q_dupuit = Q_thiem_unconfined;

% 4. Analytical Solution for Circular Excavation
% Using the formula for steady-state flow to a well in a confined aquifer
% Q = (2 * π * K * b * s) / ln(R/r_w)

Q_circular = (2 * pi * aquifer.hydraulic_conductivity * aquifer.thickness * s_w) / ...
             log(R / pump.well_radius);

% 5. Wells Method for Rectangular Excavation
% Using the equivalent radius approach
% Q = (2 * π * K * b * s) / ln(R/r_w)

Q_rectangular = (2 * pi * aquifer.hydraulic_conductivity * aquifer.thickness * s_w) / ...
                log(R / pump.well_radius);

% Adjust for multiple wells
Q_total = Q_rectangular * sim.num_wells;

% Convert to practical units
Q_thiem_confined_m3h = Q_thiem_confined * m3_per_s_to_m3_per_h;
Q_thiem_confined_m3day = Q_thiem_confined * m3_per_s_to_m3_per_day;

Q_thiem_unconfined_m3h = Q_thiem_unconfined * m3_per_s_to_m3_per_h;
Q_thiem_unconfined_m3day = Q_thiem_unconfined * m3_per_s_to_m3_per_day;

Q_dupuit_m3h = Q_dupuit * m3_per_s_to_m3_per_h;
Q_dupuit_m3day = Q_dupuit * m3_per_s_to_m3_per_day;

Q_circular_m3h = Q_circular * m3_per_s_to_m3_per_h;
Q_circular_m3day = Q_circular * m3_per_s_to_m3_per_day;

Q_rectangular_m3h = Q_rectangular * m3_per_s_to_m3_per_h;
Q_rectangular_m3day = Q_rectangular * m3_per_s_to_m3_per_day;

Q_total_m3h = Q_total * m3_per_s_to_m3_per_h;
Q_total_m3day = Q_total * m3_per_s_to_m3_per_day;

% Calculate number of wellpoints needed
well_spacing = 2 * R / sim.num_wells;
actual_wells_needed = ceil(Q_total_m3h / (pump.max_capacity * pump.well_efficiency));

% Factor of safety calculations
FOS_pumping = pump.max_capacity * actual_wells_needed * pump.well_efficiency / Q_total_m3h;

% Sensitivity analysis (vary hydraulic conductivity)
K_range = linspace(0.5e-4, 2e-4, 5);
Q_sensitivity = zeros(size(K_range));

for i = 1:length(K_range)
    aquifer_temp = aquifer;
    aquifer_temp.hydraulic_conductivity = K_range(i);
    
    Q_temp = (2 * pi * aquifer_temp.hydraulic_conductivity * aquifer_temp.thickness * s_w) / ...
             log(R / pump.well_radius);
    Q_sensitivity(i) = Q_temp * m3_per_s_to_m3_per_h;
end

% Time-dependent pumping (simplified)
% Assume pumping rate decreases as construction progresses
Q_time_dependent = Q_total_m3h * exp(-0.1 * sim.time_phases);

%% OUTPUT AND REPORTING

fprintf('\n=== EXCAVATION DEWATERING CALCULATIONS ===\n');
fprintf('Excavation Shape: %s\n', exca.shape);
fprintf('Excavation Dimensions: %.1f m x %.1f m\n', exca.width, exca.length);
fprintf('Excavation Depth: %.1f m\n', exca.depth);
fprintf('Initial Water Table: %.1f m\n', gw.initial_water_table);
fprintf('Required Drawdown: %.1f m\n', gw.drawdown_required);
fprintf('\n');

fprintf('=== AQUIFER PROPERTIES ===\n');
fprintf('Hydraulic Conductivity: %.2e m/s\n', aquifer.hydraulic_conductivity);
fprintf('Aquifer Thickness: %.1f m\n', aquifer.thickness);
fprintf('Porosity: %.2f\n', aquifer.porosity);
fprintf('Specific Yield: %.2f\n', aquifer.specific_yield);
fprintf('\n');

fprintf('=== CALCULATION RESULTS ===\n');
fprintf('Radius of Influence: %.1f m\n', R);
fprintf('Drawdown at Excavation: %.2f m\n', s_e);
fprintf('\n');

fprintf('Thiem Confined: %.2f m³/h (%.2f m³/day)\n', Q_thiem_confined_m3h, Q_thiem_confined_m3day);
fprintf('Thiem Unconfined: %.2f m³/h (%.2f m³/day)\n', Q_thiem_unconfined_m3h, Q_thiem_unconfined_m3day);
fprintf('Dupuit Method: %.2f m³/h (%.2f m³/day)\n', Q_dupuit_m3h, Q_dupuit_m3day);
fprintf('Circular Excavation: %.2f m³/h (%.2f m³/day)\n', Q_circular_m3h, Q_circular_m3day);
fprintf('Rectangular Excavation: %.2f m³/h (%.2f m³/day)\n', Q_rectangular_m3h, Q_rectangular_m3day);
fprintf('Total Pumping Rate: %.2f m³/h (%.2f m³/day)\n', Q_total_m3h, Q_total_m3day);
fprintf('\n');

fprintf('=== SYSTEM DESIGN ===\n');
fprintf('Number of Wells Required: %d\n', actual_wells_needed);
fprintf('Well Spacing: %.1f m\n', well_spacing);
fprintf('Factor of Safety: %.2f\n', FOS_pumping);
fprintf('\n');

% Create summary table
results_table = table(
    {'Thiem Confined'; 'Thiem Unconfined'; 'Dupuit'; 'Circular'; 'Rectangular'; 'Total'},
    [Q_thiem_confined_m3h; Q_thiem_unconfined_m3h; Q_dupuit_m3h; Q_circular_m3h; Q_rectangular_m3h; Q_total_m3h],
    [Q_thiem_confined_m3day; Q_thiem_unconfined_m3day; Q_dupuit_m3day; Q_circular_m3day; Q_rectangular_m3day; Q_total_m3day],
    'VariableNames', {'Method', 'm3_per_hour', 'm3_per_day'});

disp(results_table);

%% VISUALIZATION

figure('Position', [100, 100, 1200, 800]);

% 1. Drawdown Profile
subplot(2, 2, 1);
r = linspace(pump.well_radius, R, 100);

% For confined aquifer
s_confined = s_w * log(r / pump.well_radius) / log(R / pump.well_radius);

% For unconfined aquifer (Dupuit)
s_unconfined = H_initial - sqrt(H_initial^2 - (Q_thiem_unconfined / (pi * aquifer.hydraulic_conductivity)) * log(r / R));

plot(r, s_confined, 'b-', 'LineWidth', 2);
hold on;
plot(r, s_unconfined, 'r--', 'LineWidth', 2);
xlabel('Distance from Well (m)');
ylabel('Drawdown (m)');
title('Drawdown Profile');
legend('Confined Aquifer', 'Unconfined Aquifer');
grid on;

% Mark excavation perimeter
if strcmp(exca.shape, 'circular')
    xline(exca.radius, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Excavation');
else
    xline(exca.equivalent_radius, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Excavation');
end

% 2. Water Table Before/After
subplot(2, 2, 2);
x = linspace(-R, R, 100);

% Before dewatering
water_table_before = gw.initial_water_table * ones(size(x));

% After dewatering (simplified)
water_table_after = zeros(size(x));
for i = 1:length(x)
    if abs(x(i)) < exca.equivalent_radius
        water_table_after(i) = gw.initial_water_table - gw.drawdown_required;
    else
        % Linear interpolation between excavation edge and radius of influence
        if abs(x(i)) <= R
            drawdown = gw.drawdown_required * (R - abs(x(i))) / (R - exca.equivalent_radius);
            water_table_after(i) = gw.initial_water_table - drawdown;
        else
            water_table_after(i) = gw.initial_water_table;
        end
    end
end

plot(x, water_table_before, 'b-', 'LineWidth', 2);
hold on;
plot(x, water_table_after, 'r-', 'LineWidth', 2);
xlabel('Distance from Center (m)');
ylabel('Water Table Elevation (m)');
title('Water Table Before/After Dewatering');
legend('Before Dewatering', 'After Dewatering');
grid on;

% Mark excavation area
if strcmp(exca.shape, 'rectangular')
    rectangle('Position', [-exca.width/2, gw.initial_water_table - exca.depth, exca.width, exca.depth], ...
              'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1);
else
    rectangle('Position', [-exca.radius, gw.initial_water_table - exca.depth, 2*exca.radius, exca.depth], ...
              'Curvature', [1,1], 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1);
end

% 3. Excavation Geometry
subplot(2, 2, 3);
if strcmp(exca.shape, 'rectangular')
    % Draw rectangular excavation
    rectangle('Position', [-exca.width/2, -exca.depth, exca.width, exca.depth], ...
              'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;
    
    % Draw water table
    plot([-R R], [0 0], 'b-', 'LineWidth', 2);
    plot([-R R], [-gw.drawdown_required -gw.drawdown_required], 'r-', 'LineWidth', 2);
    
    % Draw wells
    for i = 1:sim.num_wells
        theta = 2 * pi * i / sim.num_wells;
        well_x = exca.equivalent_radius * cos(theta);
        well_y = -gw.drawdown_required;
        rectangle('Position', [well_x - pump.well_radius/2, well_y - pump.well_radius/2, ...
                               pump.well_radius, pump.well_radius], ...
                  'Curvature', [1,1], 'FaceColor', 'r', 'EdgeColor', 'k');
    end
    
    xlim([-R R]);
    ylim([-exca.depth - 2, 2]);
    xlabel('Distance (m)');
    ylabel('Elevation (m)');
    title('Excavation Geometry with Wells');
    grid on;
    axis equal;
else
    % Draw circular excavation
    rectangle('Position', [-exca.radius, -exca.depth, 2*exca.radius, exca.depth], ...
              'Curvature', [1,1], 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;
    
    % Draw water table
    plot([-R R], [0 0], 'b-', 'LineWidth', 2);
    plot([-R R], [-gw.drawdown_required -gw.drawdown_required], 'r-', 'LineWidth', 2);
    
    % Draw central well
    rectangle('Position', [-pump.well_radius/2, -gw.drawdown_required - pump.well_radius/2, ...
                           pump.well_radius, pump.well_radius], ...
              'Curvature', [1,1], 'FaceColor', 'r', 'EdgeColor', 'k');
    
    xlim([-R R]);
    ylim([-exca.depth - 2, 2]);
    xlabel('Distance (m)');
    ylabel('Elevation (m)');
    title('Excavation Geometry with Wells');
    grid on;
    axis equal;
end

% 4. Sensitivity Analysis
subplot(2, 2, 4);
plot(K_range * 1e4, Q_sensitivity, 'bo-', 'LineWidth', 2);
xlabel('Hydraulic Conductivity (×10^{-4} m/s)');
ylabel('Pumping Rate (m³/h)');
title('Sensitivity Analysis: Hydraulic Conductivity vs Pumping Rate');
grid on;

%% EXAMPLE CASES

% Example 1: Small residential excavation
% Inputs: 
% aquifer.hydraulic_conductivity = 5e-5; % m/s (silty sand)
% aquifer.thickness = 8; % m
% exca.shape = 'rectangular';
% exca.width = 10; % m
% exca.length = 15; % m
% exca.depth = 3; % m
% gw.drawdown_required = 2; % m
% Expected output: ~50-100 m³/h

% Example 2: Large commercial basement
% Inputs:
% aquifer.hydraulic_conductivity = 2e-4; % m/s (clean sand)
% aquifer.thickness = 12; % m
% exca.shape = 'rectangular';
% exca.width = 40; % m
% exca.length = 60; % m
% exca.depth = 6; % m
% gw.drawdown_required = 4; % m
% Expected output: ~200-400 m³/h

% Example 3: Deep circular shaft
% Inputs:
% aquifer.hydraulic_conductivity = 1.5e-4; % m/s
% aquifer.thickness = 15; % m
% exca.shape = 'circular';
% exca.radius = 8; % m
% exca.depth = 8; % m
% gw.drawdown_required = 5; % m
% Expected output: ~150-300 m³/h

%% NESTED HELPER FUNCTIONS

function Q = calculate_thiem_confined(K, b, s_w, s_e, R, r_w)
    % Thiem's formula for confined aquifer
    Q = (2 * pi * K * b * (s_w - s_e)) / log(R / r_w);
end

function Q = calculate_thiem_unconfined(K, H1, H2, R, r_w)
    % Thiem's formula for unconfined aquifer
    Q = (pi * K * (H1^2 - H2^2)) / log(R / r_w);
end

function Q = calculate_dupuit(K, H, h, R, r_w)
    % Dupuit's formula for unconfined flow
    Q = (pi * K * (H^2 - h^2)) / log(R / r_w);
end

function Q = calculate_circular_well(K, b, s, R, r_w)
    % Analytical solution for circular excavation
    Q = (2 * pi * K * b * s) / log(R / r_w);
end

function Q = calculate_rectangular_well(K, b, s, R, r_w, num_wells)
    % Wells method for rectangular excavation
    Q_single = (2 * pi * K * b * s) / log(R / r_w);
    Q = Q_single * num_wells;
end