%% RETAINING WALL GLOBAL STABILITY ANALYSIS
% This script analyzes the global stability of a cantilever retaining wall
% Checks: Sliding, Overturning, Bearing Capacity, and Global Stability

clear; clc; close all;

%% INPUT PARAMETERS

% Wall Geometry (all dimensions in meters)
H = 6.0;              % Total wall height (m)
B = 4.0;              % Base width (m)
t_stem = 0.4;         % Stem thickness at top (m)
t_base = 0.6;         % Base thickness (m)
t_toe = 1.2;          % Toe length (m)
t_heel = B - t_toe;   % Heel length (m)

% Soil Properties
gamma_soil = 18.0;    % Unit weight of backfill soil (kN/m?)
phi = 30;             % Internal friction angle of soil (degrees)
c = 0;                % Cohesion of soil (kPa)
delta = 2/3 * phi;    % Wall friction angle (degrees)

% Foundation Soil Properties
gamma_foundation = 19.0;  % Unit weight of foundation soil (kN/m?)
phi_f = 28;               % Friction angle of foundation soil (degrees)
c_f = 10;                 % Cohesion of foundation soil (kPa)
q_ult_foundation = 300;   % Ultimate bearing capacity (kPa) - if known

% Material Properties
gamma_concrete = 24.0;    % Unit weight of concrete (kN/m?)

% Surcharge Load
q = 10.0;             % Surcharge load on backfill (kPa)

% Safety Factors (Minimum Required)
FS_sliding_min = 1.5;
FS_overturning_min = 2.0;
FS_bearing_min = 3.0;

%% CALCULATIONS

fprintf('='*60); fprintf('\n');
fprintf('RETAINING WALL GLOBAL STABILITY ANALYSIS\n');
fprintf('='*60); fprintf('\n\n');

% Convert angles to radians
phi_rad = deg2rad(phi);
delta_rad = deg2rad(delta);
phi_f_rad = deg2rad(phi_f);

%% 1. ACTIVE EARTH PRESSURE COEFFICIENT (Coulomb's Theory)

% Active earth pressure coefficient
Ka = cos(phi_rad)^2 / (cos(delta_rad) * (1 + sqrt(sin(phi_rad + delta_rad) * ...
     sin(phi_rad) / cos(delta_rad)))^2);

fprintf('Active Earth Pressure Coefficient (Ka) = %.4f\n\n', Ka);

%% 2. LATERAL EARTH PRESSURE FORCES

% Active pressure from soil
Pa_soil = 0.5 * Ka * gamma_soil * H^2;  % kN/m

% Active pressure from surcharge
Pa_surcharge = Ka * q * H;  % kN/m

% Total horizontal force
Pa_total = Pa_soil + Pa_surcharge;  % kN/m

% Vertical component (due to wall friction)
Pv = Pa_total * tan(delta_rad);  % kN/m
Ph = Pa_total;  % Horizontal component

fprintf('Lateral Earth Pressure Forces:\n');
fprintf('  Pa (soil) = %.2f kN/m\n', Pa_soil);
fprintf('  Pa (surcharge) = %.2f kN/m\n', Pa_surcharge);
fprintf('  Total Horizontal Force (Ph) = %.2f kN/m\n', Ph);
fprintf('  Vertical Force (Pv) = %.2f kN/m\n\n', Pv);

%% 3. WEIGHT OF WALL COMPONENTS

% Stem weight (trapezoidal section)
W_stem = gamma_concrete * t_stem * (H - t_base);  % kN/m

% Base slab weight
W_base = gamma_concrete * B * t_base;  % kN/m

% Soil on heel weight
W_soil_heel = gamma_soil * t_heel * (H - t_base);  % kN/m

% Total vertical load
W_total = W_stem + W_base + W_soil_heel + Pv;  % kN/m

fprintf('Vertical Forces (Weights):\n');
fprintf('  Stem weight = %.2f kN/m\n', W_stem);
fprintf('  Base slab weight = %.2f kN/m\n', W_base);
fprintf('  Soil on heel = %.2f kN/m\n', W_soil_heel);
fprintf('  Vertical component of Pa = %.2f kN/m\n', Pv);
fprintf('  Total Vertical Force = %.2f kN/m\n\n', W_total);

%% 4. MOMENT ARMS (about toe)

% Distance from toe to center of each force
x_stem = t_toe + t_stem/2;
x_base = B/2;
x_soil_heel = t_toe + t_stem + t_heel/2;
x_Pv = B;  % Acts at heel
y_Ph = H/3;  % Active pressure acts at H/3 from base

fprintf('Moment Arms (from toe):\n');
fprintf('  Stem: %.2f m\n', x_stem);
fprintf('  Base: %.2f m\n', x_base);
fprintf('  Soil on heel: %.2f m\n', x_soil_heel);
fprintf('  Horizontal force height: %.2f m\n\n', y_Ph);

%% 5. STABILITY AGAINST OVERTURNING

% Resisting moments
MR_stem = W_stem * x_stem;
MR_base = W_base * x_base;
MR_soil = W_soil_heel * x_soil_heel;
MR_Pv = Pv * x_Pv;
MR_total = MR_stem + MR_base + MR_soil + MR_Pv;

% Overturning moment
MO = Ph * y_Ph;

% Factor of safety against overturning
FS_overturning = MR_total / MO;

fprintf('OVERTURNING STABILITY:\n');
fprintf('  Resisting Moment = %.2f kN?m/m\n', MR_total);
fprintf('  Overturning Moment = %.2f kN?m/m\n', MO);
fprintf('  FS (Overturning) = %.2f', FS_overturning);
if FS_overturning >= FS_overturning_min
    fprintf(' [OK] ?\n');
else
    fprintf(' [UNSAFE] ? (Required: %.2f)\n', FS_overturning_min);
end
fprintf('\n');

%% 6. STABILITY AGAINST SLIDING

% Passive resistance (usually neglected for conservative design)
Kp = tan(deg2rad(45 + phi_f/2))^2;
Pp = 0;  % Conservatively neglected

% Resisting force (friction + cohesion)
Fr = W_total * tan(phi_f_rad) + c_f * B;

% Driving force
Fd = Ph;

% Factor of safety against sliding
FS_sliding = Fr / Fd;

fprintf('SLIDING STABILITY:\n');
fprintf('  Resisting Force = %.2f kN/m\n', Fr);
fprintf('  Driving Force = %.2f kN/m\n', Fd);
fprintf('  FS (Sliding) = %.2f', FS_sliding);
if FS_sliding >= FS_sliding_min
    fprintf(' [OK] ?\n');
else
    fprintf(' [UNSAFE] ? (Required: %.2f)\n', FS_sliding_min);
end
fprintf('\n');

%% 7. BEARING CAPACITY CHECK

% Resultant force location
x_resultant = (MR_total - MO) / W_total;

% Eccentricity
e = B/2 - x_resultant;

fprintf('BEARING CAPACITY:\n');
fprintf('  Resultant location from toe = %.2f m\n', x_resultant);
fprintf('  Eccentricity (e) = %.2f m\n', e);
fprintf('  B/6 = %.2f m\n', B/6);

% Check if resultant is within middle third
if abs(e) <= B/6
    fprintf('  Resultant within middle third [OK] ?\n');
    middle_third_ok = true;
else
    fprintf('  Resultant outside middle third [WARNING] ?\n');
    middle_third_ok = false;
end

% Pressure distribution
if abs(e) <= B/6
    % No tension - trapezoidal distribution
    q_max = (W_total / B) * (1 + 6*e/B);
    q_min = (W_total / B) * (1 - 6*e/B);
else
    % Tension exists - triangular distribution
    L_eff = 3 * (B/2 - e);  % Effective length
    q_max = 2 * W_total / L_eff;
    q_min = 0;
end

fprintf('  Maximum bearing pressure = %.2f kPa\n', q_max);
fprintf('  Minimum bearing pressure = %.2f kPa\n', q_min);

% Calculate allowable bearing capacity using Terzaghi's equation
% Bearing capacity factors
Nq = exp(pi * tan(phi_f_rad)) * tan(deg2rad(45 + phi_f/2))^2;
Nc = (Nq - 1) / tan(phi_f_rad);
Ngamma = 2 * (Nq - 1) * tan(phi_f_rad);

% Ultimate bearing capacity (simplified - for strip footing)
q_ult_calc = c_f * Nc + gamma_foundation * (B - 2*abs(e)) * 0.5 * Ngamma + ...
             gamma_foundation * t_base * Nq;

fprintf('  Calculated ultimate bearing capacity = %.2f kPa\n', q_ult_calc);

% Factor of safety for bearing capacity
FS_bearing = q_ult_calc / q_max;

fprintf('  FS (Bearing Capacity) = %.2f', FS_bearing);
if FS_bearing >= FS_bearing_min
    fprintf(' [OK] ?\n');
else
    fprintf(' [UNSAFE] ? (Required: %.2f)\n', FS_bearing_min);
end
fprintf('\n');

%% 8. GLOBAL STABILITY SUMMARY

fprintf('='*60); fprintf('\n');
fprintf('STABILITY SUMMARY\n');
fprintf('='*60); fprintf('\n');

fprintf('Check                    FS Required    FS Actual    Status\n');
fprintf('-'*60); fprintf('\n');

fprintf('Overturning             %.2f          %.2f        ', ...
        FS_overturning_min, FS_overturning);
if FS_overturning >= FS_overturning_min
    fprintf('[PASS]\n');
else
    fprintf('[FAIL]\n');
end

fprintf('Sliding                 %.2f          %.2f        ', ...
        FS_sliding_min, FS_sliding);
if FS_sliding >= FS_sliding_min
    fprintf('[PASS]\n');
else
    fprintf('[FAIL]\n');
end

fprintf('Bearing Capacity        %.2f          %.2f        ', ...
        FS_bearing_min, FS_bearing);
if FS_bearing >= FS_bearing_min
    fprintf('[PASS]\n');
else
    fprintf('[FAIL]\n');
end

fprintf('Middle Third Rule       N/A           N/A         ');
if middle_third_ok
    fprintf('[PASS]\n');
else
    fprintf('[FAIL]\n');
end

fprintf('='*60); fprintf('\n\n');

% Overall stability verdict
all_checks_pass = (FS_overturning >= FS_overturning_min) && ...
                  (FS_sliding >= FS_sliding_min) && ...
                  (FS_bearing >= FS_bearing_min) && ...
                  middle_third_ok;

if all_checks_pass
    fprintf('OVERALL VERDICT: RETAINING WALL IS STABLE ?\n');
else
    fprintf('OVERALL VERDICT: RETAINING WALL REQUIRES REDESIGN ?\n');
end

%% 9. VISUALIZATION

figure('Position', [100, 100, 1200, 600]);

% Subplot 1: Wall Geometry
subplot(1,2,1);
hold on; grid on; axis equal;

% Draw wall
wall_x = [t_toe, t_toe, t_toe+t_stem, t_toe+t_stem, B, B, 0, 0, t_toe];
wall_y = [0, t_base, t_base, H, H, 0, 0, t_base, 0];
fill(wall_x, wall_y, [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 2);

% Draw backfill soil
soil_x = [t_toe+t_stem, t_toe+t_stem, B+2, B+2, B, B, t_toe+t_stem];
soil_y = [t_base, H, H, 0, 0, t_base, t_base];
fill(soil_x, soil_y, [0.8 0.6 0.4], 'EdgeColor', 'k', 'LineWidth', 1);

% Draw forces
scale = 0.05;  % Scale for force arrows

% Horizontal force
quiver(t_toe+t_stem, y_Ph, Ph*scale, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(t_toe+t_stem+Ph*scale+0.3, y_Ph, sprintf('Pa=%.1f kN/m', Ph), 'Color', 'r');

% Weights (simplified - show resultant)
quiver(x_resultant, H+1, 0, -0.8, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(x_resultant, H+1.3, sprintf('W=%.1f kN/m', W_total), 'Color', 'b', 'HorizontalAlignment', 'center');

% Mark toe (origin)
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(0, -0.5, 'Toe (Origin)', 'HorizontalAlignment', 'center');

% Dimensions
plot([0, B], [-1, -1], 'k-', 'LineWidth', 1);
text(B/2, -1.5, sprintf('B = %.1f m', B), 'HorizontalAlignment', 'center');

xlabel('Distance (m)'); ylabel('Height (m)');
title('Retaining Wall Geometry and Forces');
xlim([-1, B+3]); ylim([-2, H+2]);

% Subplot 2: Bearing Pressure Distribution
subplot(1,2,2);
hold on; grid on;

% Base of wall
plot([0, B], [0, 0], 'k-', 'LineWidth', 3);

% Pressure distribution
if middle_third_ok
    pressure_x = [0, B, B, 0, 0];
    pressure_y = [0, 0, -q_max*0.02, -q_min*0.02, 0];
    fill(pressure_x, pressure_y, [0.3 0.6 0.9], 'EdgeColor', 'b', 'LineWidth', 2);
else
    L_eff = 3 * (B/2 - e);
    x_start = B/2 - L_eff/2;
    pressure_x = [x_start, x_start+L_eff, x_start+L_eff, x_start, x_start];
    pressure_y = [0, 0, -q_max*0.02, 0, 0];
    fill(pressure_x, pressure_y, [0.9 0.6 0.3], 'EdgeColor', 'r', 'LineWidth', 2);
end

% Resultant location
plot(x_resultant, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(x_resultant, 0.5, 'Resultant', 'HorizontalAlignment', 'center', 'Color', 'r');

% Middle third boundaries
plot([B/6, B/6], [-5, 1], 'g--', 'LineWidth', 1.5);
plot([5*B/6, 5*B/6], [-5, 1], 'g--', 'LineWidth', 1.5);
text(B/6, 1.5, 'B/6', 'Color', 'g', 'HorizontalAlignment', 'center');
text(5*B/6, 1.5, '5B/6', 'Color', 'g', 'HorizontalAlignment', 'center');

% Annotations
text(B/2, -6, sprintf('q_{max} = %.2f kPa', q_max), 'HorizontalAlignment', 'center');
text(B/2, -7, sprintf('q_{min} = %.2f kPa', q_min), 'HorizontalAlignment', 'center');
text(B/2, -8, sprintf('e = %.2f m', e), 'HorizontalAlignment', 'center');

xlabel('Base Width (m)'); ylabel('Pressure (kPa)');
title('Bearing Pressure Distribution');
xlim([-0.5, B+0.5]); ylim([-10, 2]);

sgtitle('Retaining Wall Global Stability Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Analysis complete. Check the figure for visualization.\n');
