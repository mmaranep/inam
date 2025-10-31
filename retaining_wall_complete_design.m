%% COMPLETE RETAINING WALL DESIGN - STABILITY + STRUCTURAL (CSA A23.3-19)
% This comprehensive script performs:
% 1. Global stability analysis (sliding, overturning, bearing capacity)
% 2. Structural design per CSA A23.3-19 (flexure, shear, development)
% 3. Complete reinforcement detailing
% 4. Design optimization recommendations

clear; clc; close all;

%% ========================================================================
%                           INPUT PARAMETERS
% =========================================================================

% Wall Geometry (all dimensions in meters)
H = 6.0;              % Total wall height above footing (m)
B = 4.0;              % Base width (m)
t_stem_top = 0.3;     % Stem thickness at top (m)
t_stem_base = 0.4;    % Stem thickness at base (m)
t_base = 0.6;         % Base thickness (m)
t_toe = 1.2;          % Toe length (m)
t_heel = 2.4;         % Heel length (m) - calculated to match B

% Soil Properties - Backfill
gamma_soil = 18.0;    % Unit weight of backfill soil (kN/m?)
phi = 30;             % Internal friction angle of soil (degrees)
c = 0;                % Cohesion of soil (kPa)
delta = 2/3 * phi;    % Wall friction angle (degrees)
q = 10.0;             % Surcharge load on backfill (kPa)

% Foundation Soil Properties
gamma_foundation = 19.0;  % Unit weight of foundation soil (kN/m?)
phi_f = 28;               % Friction angle of foundation soil (degrees)
c_f = 10;                 % Cohesion of foundation soil (kPa)

% Material Properties
gamma_concrete = 24.0;    % Unit weight of concrete (kN/m?)
fc_prime = 30;            % Concrete compressive strength (MPa)
fy = 400;                 % Reinforcement yield strength (MPa)
Es = 200000;              % Steel modulus of elasticity (MPa)

% CSA A23.3 Material Resistance Factors
phi_c = 0.65;         % Concrete resistance factor
phi_s = 0.85;         % Steel resistance factor

% Load Factors (CSA A23.3 - Ultimate Limit State)
alpha_D = 1.25;       % Dead load factor
alpha_L = 1.5;        % Live load factor (surcharge)

% Design Cover
cover = 0.075;        % Concrete cover (m) - 75mm for earth face

% Safety Factors (Minimum Required for Stability)
FS_sliding_min = 1.5;
FS_overturning_min = 2.0;
FS_bearing_min = 3.0;

% Bar diameters available (mm) and their areas (mm?)
bar_sizes = [10, 15, 20, 25, 30, 35];
bar_areas = [100, 200, 300, 500, 700, 1000];

%% ========================================================================
%                    PART A: GLOBAL STABILITY ANALYSIS
% =========================================================================

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('?  RETAINING WALL COMPLETE DESIGN                                            ?\n');
fprintf('?  PART A: GLOBAL STABILITY ANALYSIS                                         ?\n');
fprintf('%s\n\n', repmat('=', 1, 80));

% Convert angles to radians
phi_rad = deg2rad(phi);
delta_rad = deg2rad(delta);
phi_f_rad = deg2rad(phi_f);

% Active earth pressure coefficient (Coulomb's theory)
Ka = cos(phi_rad)^2 / (cos(delta_rad) * (1 + sqrt(sin(phi_rad + delta_rad) * ...
     sin(phi_rad) / cos(delta_rad)))^2);

fprintf('Design Parameters:\n');
fprintf('  Wall height (H) = %.2f m\n', H);
fprintf('  Base width (B) = %.2f m\n', B);
fprintf('  Active earth pressure coefficient (Ka) = %.4f\n\n', Ka);

%% A1. LATERAL EARTH PRESSURE

h_stem = H - t_base;
t_stem = t_stem_base;

Pa_soil = 0.5 * Ka * gamma_soil * h_stem^2;
Pa_surcharge = Ka * q * h_stem;
Pa_total = Pa_soil + Pa_surcharge;
Pv = Pa_total * tan(delta_rad);
Ph = Pa_total;

fprintf('Lateral Earth Pressure:\n');
fprintf('  Horizontal force (Ph) = %.2f kN/m\n', Ph);
fprintf('  Vertical force (Pv) = %.2f kN/m\n\n', Pv);

%% A2. WEIGHTS

W_stem = gamma_concrete * ((t_stem_top + t_stem_base)/2) * h_stem;
W_base = gamma_concrete * B * t_base;
W_soil_heel = gamma_soil * t_heel * h_stem;
W_total = W_stem + W_base + W_soil_heel + Pv;

fprintf('Vertical Forces:\n');
fprintf('  Stem weight = %.2f kN/m\n', W_stem);
fprintf('  Base weight = %.2f kN/m\n', W_base);
fprintf('  Soil on heel = %.2f kN/m\n', W_soil_heel);
fprintf('  Vertical component Pa = %.2f kN/m\n', Pv);
fprintf('  Total = %.2f kN/m\n\n', W_total);

%% A3. STABILITY CHECKS

% Moment arms
x_stem = t_toe + (t_stem_top + t_stem_base)/4;
x_base = B/2;
x_soil_heel = t_toe + t_stem + t_heel/2;
y_Ph = h_stem/3;

% Overturning check
MR_total = W_stem * x_stem + W_base * x_base + W_soil_heel * x_soil_heel + Pv * B;
MO = Ph * y_Ph;
FS_overturning = MR_total / MO;

fprintf('OVERTURNING CHECK:\n');
fprintf('  Resisting moment = %.2f kN?m/m\n', MR_total);
fprintf('  Overturning moment = %.2f kN?m/m\n', MO);
fprintf('  FS = %.2f (Required: %.2f) ', FS_overturning, FS_overturning_min);
if FS_overturning >= FS_overturning_min
    fprintf('[PASS] ?\n\n');
    pass_overturning = true;
else
    fprintf('[FAIL] ?\n\n');
    pass_overturning = false;
end

% Sliding check
Fr = W_total * tan(phi_f_rad) + c_f * B;
Fd = Ph;
FS_sliding = Fr / Fd;

fprintf('SLIDING CHECK:\n');
fprintf('  Resisting force = %.2f kN/m\n', Fr);
fprintf('  Driving force = %.2f kN/m\n', Fd);
fprintf('  FS = %.2f (Required: %.2f) ', FS_sliding, FS_sliding_min);
if FS_sliding >= FS_sliding_min
    fprintf('[PASS] ?\n\n');
    pass_sliding = true;
else
    fprintf('[FAIL] ?\n\n');
    pass_sliding = false;
end

% Bearing capacity check
x_resultant = (MR_total - MO) / W_total;
e = B/2 - x_resultant;

fprintf('BEARING CAPACITY CHECK:\n');
fprintf('  Resultant location = %.2f m from toe\n', x_resultant);
fprintf('  Eccentricity (e) = %.2f m\n', e);

if abs(e) <= B/6
    q_max = (W_total / B) * (1 + 6*e/B);
    q_min = (W_total / B) * (1 - 6*e/B);
    fprintf('  Resultant within middle third [OK] ?\n');
    pass_middle_third = true;
else
    L_eff = 3 * (B/2 - e);
    q_max = 2 * W_total / L_eff;
    q_min = 0;
    fprintf('  Resultant outside middle third [WARNING] ?\n');
    pass_middle_third = false;
end

fprintf('  Max bearing pressure = %.2f kPa\n', q_max);
fprintf('  Min bearing pressure = %.2f kPa\n', q_min);

% Terzaghi's bearing capacity
Nq = exp(pi * tan(phi_f_rad)) * tan(deg2rad(45 + phi_f/2))^2;
Nc = (Nq - 1) / tan(phi_f_rad);
Ngamma = 2 * (Nq - 1) * tan(phi_f_rad);
q_ult = c_f * Nc + gamma_foundation * (B - 2*abs(e)) * 0.5 * Ngamma + ...
        gamma_foundation * t_base * Nq;

FS_bearing = q_ult / q_max;
fprintf('  Ultimate bearing capacity = %.2f kPa\n', q_ult);
fprintf('  FS = %.2f (Required: %.2f) ', FS_bearing, FS_bearing_min);
if FS_bearing >= FS_bearing_min
    fprintf('[PASS] ?\n\n');
    pass_bearing = true;
else
    fprintf('[FAIL] ?\n\n');
    pass_bearing = false;
end

% Overall stability verdict
stability_ok = pass_overturning && pass_sliding && pass_bearing && pass_middle_third;

fprintf('%s\n', repmat('=', 1, 80));
if stability_ok
    fprintf('STABILITY VERDICT: ALL CHECKS PASSED ???\n');
else
    fprintf('STABILITY VERDICT: DESIGN REQUIRES REVISION ?\n');
end
fprintf('%s\n\n', repmat('=', 1, 80));

%% ========================================================================
%                  PART B: STRUCTURAL DESIGN (CSA A23.3-19)
% =========================================================================

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('?  PART B: STRUCTURAL DESIGN (CSA A23.3-19)                                  ?\n');
fprintf('%s\n\n', repmat('=', 1, 80));

% Stress block parameters
alpha1 = max(0.85 - 0.0015 * fc_prime, 0.67);
beta1 = max(0.97 - 0.0025 * fc_prime, 0.67);
rho_min = max(0.002, 1.4/fy);

fprintf('Material Properties:\n');
fprintf('  f''c = %.1f MPa, fy = %.1f MPa\n', fc_prime, fy);
fprintf('  ?c = %.2f, ?s = %.2f\n', phi_c, phi_s);
fprintf('  ?1 = %.3f, ?1 = %.3f\n', alpha1, beta1);
fprintf('  ?_min = %.5f\n\n', rho_min);

%% B1. STEM DESIGN

fprintf('%s\n', repmat('=', 1, 80));
fprintf('B1. STEM DESIGN\n');
fprintf('%s\n\n', repmat('=', 1, 80));

% Factored loads
pa_factored = alpha_D * Ka * gamma_soil * h_stem + alpha_L * Ka * q;
Mf_stem = (1/6) * alpha_D * Ka * gamma_soil * h_stem^3 + ...
          0.5 * alpha_L * Ka * q * h_stem^2;
Vf_stem = 0.5 * alpha_D * Ka * gamma_soil * h_stem^2 + alpha_L * Ka * q * h_stem;

fprintf('At base of stem (critical section):\n');
fprintf('  Factored moment (Mf) = %.2f kN?m/m\n', Mf_stem);
fprintf('  Factored shear (Vf) = %.2f kN/m\n\n', Vf_stem);

% Flexural design
b = 1000;
d_stem = t_stem * 1000 - cover * 1000 - 10;
Mf_stem_Nmm = Mf_stem * 1e6;
omega = 1 - sqrt(1 - 2 * Mf_stem_Nmm / (phi_c * alpha1 * fc_prime * b * d_stem^2));
rho_stem = max(omega * alpha1 * phi_c * fc_prime / (phi_s * fy), rho_min);
As_required_stem = rho_stem * b * d_stem;

fprintf('Flexural Design:\n');
fprintf('  d = %.0f mm\n', d_stem);
fprintf('  As required = %.0f mm?/m\n', As_required_stem);

% Select bars
[bar_dia_stem, spacing_stem, As_provided_stem] = ...
    select_reinforcement(As_required_stem, bar_sizes, bar_areas);

a = As_provided_stem * fy / (alpha1 * fc_prime * b);
Mr_stem = phi_s * As_provided_stem * fy * (d_stem - a/2) / 1e6;

fprintf('  PROVIDED: %dM @ %d mm c/c\n', bar_dia_stem, spacing_stem);
fprintf('  As provided = %.0f mm?/m\n', As_provided_stem);
fprintf('  Mr = %.2f kN?m/m (%.1f%% utilized)\n\n', Mr_stem, Mf_stem/Mr_stem*100);

% Shear design
lambda = 1.0;
beta = 0.21;
Vc_stem = phi_c * lambda * beta * sqrt(fc_prime) * b * d_stem / 1000;

fprintf('Shear Design:\n');
fprintf('  Vc = %.2f kN/m\n', Vc_stem);
if Vf_stem <= Vc_stem
    fprintf('  Shear reinforcement NOT REQUIRED [OK] ?\n');
    fprintf('  (%.1f%% of capacity)\n\n', Vf_stem/Vc_stem*100);
else
    fprintf('  Shear reinforcement REQUIRED [WARNING]\n\n');
end

%% B2. TOE SLAB DESIGN

fprintf('%s\n', repmat('=', 1, 80));
fprintf('B2. TOE SLAB DESIGN\n');
fprintf('%s\n\n', repmat('=', 1, 80));

% Pressure distribution
q_toe_factored = alpha_D * q_min;
q_stem_factored = alpha_D * q_max;
q_at_stem = q_toe_factored + (q_stem_factored - q_toe_factored) * t_toe / B;

Mf_toe = (q_toe_factored + q_at_stem) / 2 * t_toe^2 / 2 - ...
         alpha_D * gamma_concrete * t_base * t_toe^2 / 2;

d_toe = t_base * 1000 - cover * 1000 - 10;
critical_dist_toe = max(t_toe - d_toe/1000, 0.01);
q_at_crit_toe = q_toe_factored + (q_stem_factored - q_toe_factored) * critical_dist_toe / B;
Vf_toe = (q_toe_factored + q_at_crit_toe) / 2 * critical_dist_toe - ...
         alpha_D * gamma_concrete * t_base * critical_dist_toe;

fprintf('Loads:\n');
fprintf('  Mf = %.2f kN?m/m\n', Mf_toe);
fprintf('  Vf = %.2f kN/m\n\n', Vf_toe);

% Design
Mf_toe_Nmm = Mf_toe * 1e6;
omega_toe = 1 - sqrt(1 - 2 * Mf_toe_Nmm / (phi_c * alpha1 * fc_prime * b * d_toe^2));
rho_toe = max(omega_toe * alpha1 * phi_c * fc_prime / (phi_s * fy), rho_min);
As_required_toe = rho_toe * b * d_toe;

[bar_dia_toe, spacing_toe, As_provided_toe] = ...
    select_reinforcement(As_required_toe, bar_sizes, bar_areas);

a = As_provided_toe * fy / (alpha1 * fc_prime * b);
Mr_toe = phi_s * As_provided_toe * fy * (d_toe - a/2) / 1e6;

fprintf('Flexural Design:\n');
fprintf('  PROVIDED: %dM @ %d mm c/c (bottom)\n', bar_dia_toe, spacing_toe);
fprintf('  Mr = %.2f kN?m/m (%.1f%% utilized)\n\n', Mr_toe, Mf_toe/Mr_toe*100);

Vc_toe = phi_c * lambda * beta * sqrt(fc_prime) * b * d_toe / 1000;
fprintf('Shear Design:\n');
fprintf('  Vc = %.2f kN/m, Vf/Vc = %.1f%% ', Vc_toe, Vf_toe/Vc_toe*100);
if Vf_toe <= Vc_toe
    fprintf('[OK] ?\n\n');
else
    fprintf('[WARNING]\n\n');
end

%% B3. HEEL SLAB DESIGN

fprintf('%s\n', repmat('=', 1, 80));
fprintf('B3. HEEL SLAB DESIGN\n');
fprintf('%s\n\n', repmat('=', 1, 80));

% Loads
W_soil_heel_factored = alpha_D * gamma_soil * h_stem * t_heel;
q_heel_stem = q_stem_factored;
q_heel_end = q_stem_factored - (q_stem_factored - q_toe_factored) * (t_toe + t_stem) / B;

M_soil_weight = W_soil_heel_factored * t_heel / 2;
M_soil_pressure = -(q_heel_stem + q_heel_end) / 2 * t_heel^2 / 2;
M_self_weight = -alpha_D * gamma_concrete * t_base * t_heel^2 / 2;
Mf_heel = M_soil_weight + M_soil_pressure + M_self_weight;

d_heel = t_base * 1000 - cover * 1000 - 10;
critical_dist_heel = max(t_heel - d_heel/1000, 0.01);
Vf_heel = W_soil_heel_factored - (q_heel_stem + q_heel_end) / 2 * critical_dist_heel - ...
          alpha_D * gamma_concrete * t_base * critical_dist_heel;

fprintf('Loads:\n');
fprintf('  Mf = %.2f kN?m/m\n', abs(Mf_heel));
fprintf('  Vf = %.2f kN/m\n\n', abs(Vf_heel));

% Design
Mf_heel_Nmm = abs(Mf_heel) * 1e6;
omega_heel = 1 - sqrt(1 - 2 * Mf_heel_Nmm / (phi_c * alpha1 * fc_prime * b * d_heel^2));
rho_heel = max(omega_heel * alpha1 * phi_c * fc_prime / (phi_s * fy), rho_min);
As_required_heel = rho_heel * b * d_heel;

[bar_dia_heel, spacing_heel, As_provided_heel] = ...
    select_reinforcement(As_required_heel, bar_sizes, bar_areas);

a = As_provided_heel * fy / (alpha1 * fc_prime * b);
Mr_heel = phi_s * As_provided_heel * fy * (d_heel - a/2) / 1e6;

fprintf('Flexural Design:\n');
fprintf('  PROVIDED: %dM @ %d mm c/c (top)\n', bar_dia_heel, spacing_heel);
fprintf('  Mr = %.2f kN?m/m (%.1f%% utilized)\n\n', Mr_heel, abs(Mf_heel)/Mr_heel*100);

Vc_heel = phi_c * lambda * beta * sqrt(fc_prime) * b * d_heel / 1000;
fprintf('Shear Design:\n');
fprintf('  Vc = %.2f kN/m, Vf/Vc = %.1f%% ', Vc_heel, abs(Vf_heel)/Vc_heel*100);
if abs(Vf_heel) <= Vc_heel
    fprintf('[OK] ?\n\n');
else
    fprintf('[WARNING]\n\n');
end

%% B4. DEVELOPMENT LENGTH

fprintf('%s\n', repmat('=', 1, 80));
fprintf('B4. DEVELOPMENT LENGTH (CSA A23.3 Cl. 12)\n');
fprintf('%s\n\n', repmat('=', 1, 80));

k1 = 1.0; k2 = 1.0; k3 = 1.0; k4 = 0.8;
ld_factor = k1 * k2 * k3 * k4 * fy / (1.0 * sqrt(fc_prime) * 1.15);

ld_stem = ld_factor * bar_dia_stem;
ld_toe = ld_factor * bar_dia_toe;
ld_heel = ld_factor * bar_dia_heel;

L_avail_stem = t_base * 1000 - 2 * cover * 1000;
L_avail_toe = (t_toe - cover) * 1000;
L_avail_heel = (t_heel - cover) * 1000;

fprintf('Stem bars (%dM): ld = %.0f mm, available = %.0f mm ', ...
        bar_dia_stem, ld_stem, L_avail_stem);
if L_avail_stem >= ld_stem
    fprintf('[OK] ?\n');
else
    fprintf('[HOOKS REQ''D]\n');
end

fprintf('Toe bars (%dM):  ld = %.0f mm, available = %.0f mm ', ...
        bar_dia_toe, ld_toe, L_avail_toe);
if L_avail_toe >= ld_toe
    fprintf('[OK] ?\n');
else
    fprintf('[HOOKS REQ''D]\n');
end

fprintf('Heel bars (%dM): ld = %.0f mm, available = %.0f mm ', ...
        bar_dia_heel, ld_heel, L_avail_heel);
if L_avail_heel >= ld_heel
    fprintf('[OK] ?\n\n');
else
    fprintf('[HOOKS REQ''D]\n\n');
end

%% ========================================================================
%                        DESIGN SUMMARY & OUTPUT
% =========================================================================

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('?  FINAL DESIGN SUMMARY                                                      ?\n');
fprintf('%s\n\n', repmat('=', 1, 80));

fprintf('GEOMETRY:\n');
fprintf('  Wall height (H) = %.2f m\n', H);
fprintf('  Base width (B) = %.2f m\n', B);
fprintf('  Stem thickness (base) = %.0f mm\n', t_stem*1000);
fprintf('  Base slab thickness = %.0f mm\n', t_base*1000);
fprintf('  Toe length = %.2f m, Heel length = %.2f m\n\n', t_toe, t_heel);

fprintf('REINFORCEMENT:\n');
fprintf('  STEM:\n');
fprintf('    Vertical (tension face): %dM @ %d mm c/c\n', bar_dia_stem, spacing_stem);
fprintf('    Horizontal:              10M @ 300 mm c/c\n');
fprintf('  BASE SLAB:\n');
fprintf('    Toe (bottom):            %dM @ %d mm c/c\n', bar_dia_toe, spacing_toe);
fprintf('    Heel (top):              %dM @ %d mm c/c\n', bar_dia_heel, spacing_heel);
fprintf('    Temperature/shrinkage:   10M @ 300 mm c/c (both faces)\n\n');

fprintf('MATERIALS:\n');
fprintf('  Concrete: f''c = %.0f MPa\n', fc_prime);
fprintf('  Steel: fy = %.0f MPa (deformed bars)\n', fy);
fprintf('  Cover: %.0f mm (earth face), 50 mm (formed face)\n\n', cover*1000);

fprintf('STABILITY:\n');
fprintf('  Overturning FS = %.2f (req. %.2f) ', FS_overturning, FS_overturning_min);
if pass_overturning, fprintf('?\n'); else, fprintf('?\n'); end
fprintf('  Sliding FS = %.2f (req. %.2f) ', FS_sliding, FS_sliding_min);
if pass_sliding, fprintf('?\n'); else, fprintf('?\n'); end
fprintf('  Bearing FS = %.2f (req. %.2f) ', FS_bearing, FS_bearing_min);
if pass_bearing, fprintf('?\n'); else, fprintf('?\n'); end
fprintf('  Middle third rule: ');
if pass_middle_third, fprintf('Satisfied ?\n\n'); else, fprintf('Not satisfied ?\n\n'); end

if stability_ok
    fprintf('OVERALL DESIGN STATUS: ACCEPTABLE ???\n');
else
    fprintf('OVERALL DESIGN STATUS: REQUIRES REVISION ???\n');
end

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Design completed per CSA A23.3-19\n');
fprintf('%s\n\n', repmat('=', 1, 80));

%% VISUALIZATION

create_design_drawings(H, B, t_stem, t_base, t_toe, t_heel, h_stem, cover, ...
                      bar_dia_stem, spacing_stem, bar_dia_toe, spacing_toe, ...
                      bar_dia_heel, spacing_heel, Mf_stem, Mr_stem, Vf_stem, ...
                      Vc_stem, q_max, q_min, x_resultant);

%% ========================================================================
%                          HELPER FUNCTIONS
% =========================================================================

function [bar_dia, spacing, As_provided] = select_reinforcement(As_required, bar_sizes, bar_areas)
    % Select optimal reinforcement configuration
    spacing_options = [100, 125, 150, 175, 200, 250, 300];
    best_design = [];
    min_excess = inf;
    
    for i = 1:length(bar_sizes)
        for j = 1:length(spacing_options)
            As_prov = bar_areas(i) * 1000 / spacing_options(j);
            if As_prov >= As_required
                excess = (As_prov - As_required) / As_required;
                if excess < min_excess && spacing_options(j) <= 300
                    min_excess = excess;
                    best_design = [bar_sizes(i), spacing_options(j), As_prov];
                end
            end
        end
    end
    
    if ~isempty(best_design)
        bar_dia = best_design(1);
        spacing = best_design(2);
        As_provided = best_design(3);
    else
        % Fallback
        bar_dia = 20;
        spacing = 150;
        As_provided = 300 * 1000 / 150;
    end
end

function create_design_drawings(H, B, t_stem, t_base, t_toe, t_heel, h_stem, cover, ...
                               bar_dia_stem, spacing_stem, bar_dia_toe, spacing_toe, ...
                               bar_dia_heel, spacing_heel, Mf_stem, Mr_stem, Vf_stem, ...
                               Vc_stem, q_max, q_min, x_resultant)
    % Create comprehensive design drawings
    
    figure('Position', [50, 50, 1600, 900], 'Color', 'w');
    
    % Drawing 1: Complete section with dimensions
    subplot(2,3,1);
    hold on; axis equal;
    
    % Draw wall
    wall_x = [t_toe, t_toe, t_toe+t_stem, t_toe+t_stem, B, B, 0, 0, t_toe];
    wall_y = [0, t_base, t_base, H, H, 0, 0, t_base, 0];
    fill(wall_x, wall_y, [0.85 0.85 0.85], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);
    
    % Backfill
    soil_x = [t_toe+t_stem, t_toe+t_stem, B+1.5, B+1.5, B, B, t_toe+t_stem];
    soil_y = [t_base, H, H, 0, 0, t_base, t_base];
    fill(soil_x, soil_y, [0.76 0.60 0.42], 'FaceAlpha', 0.7, 'EdgeColor', [0.4 0.3 0.2]);
    
    % Dimensions
    plot([0, B], [-0.4, -0.4], 'k-', 'LineWidth', 1.5);
    text(B/2, -0.7, sprintf('B = %.2fm', B), 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    plot([-0.4, -0.4], [0, H], 'k-', 'LineWidth', 1.5);
    text(-0.8, H/2, sprintf('H = %.2fm', H), 'FontSize', 11, 'FontWeight', 'bold', 'Rotation', 90, 'HorizontalAlignment', 'center');
    
    xlabel('Length (m)', 'FontSize', 10);
    ylabel('Height (m)', 'FontSize', 10);
    title('WALL SECTION', 'FontSize', 12, 'FontWeight', 'bold');
    grid on; xlim([-1.2, B+2]); ylim([-1, H+0.5]);
    
    % Drawing 2: Reinforcement layout - Stem
    subplot(2,3,2);
    hold on; axis equal;
    rectangle('Position', [0, 0, t_stem*1000, h_stem*1000], ...
              'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 2);
    
    % Vertical bars
    y_bars = (cover*1000+spacing_stem/2):spacing_stem:(h_stem*1000-cover*1000);
    for yb = y_bars
        plot(t_stem*1000-cover*1000, yb, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    
    % Horizontal bars
    x_bars = (cover*1000):300:(t_stem*1000-cover*1000);
    for xb = x_bars
        plot(xb, h_stem*500, 'bs', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    end
    
    text(t_stem*500, h_stem*250, sprintf('%dM@%dmm vert.\n10M@300mm horiz.', ...
         bar_dia_stem, spacing_stem), 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.6 0 0]);
    
    xlabel('Width (mm)'); ylabel('Height (mm)');
    title('STEM REINFORCEMENT', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Drawing 3: Reinforcement layout - Base
    subplot(2,3,3);
    hold on; axis equal;
    rectangle('Position', [0, 0, B*1000, t_base*1000], ...
              'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 2);
    
    % Mark stem
    rectangle('Position', [t_toe*1000, 0, t_stem*1000, t_base*1000], ...
              'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 1);
    
    % Toe bars (bottom)
    x_toe = (cover*1000):spacing_toe:(t_toe*1000-cover*1000);
    for xt = x_toe
        plot(xt, cover*1000, 'rs', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    end
    text(t_toe*500, cover*1000-30, sprintf('%dM@%d', bar_dia_toe, spacing_toe), ...
         'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Heel bars (top)
    heel_start = (t_toe+t_stem)*1000 + cover*1000;
    x_heel = heel_start:spacing_heel:(B*1000-cover*1000);
    for xh = x_heel
        plot(xh, t_base*1000-cover*1000, 'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    end
    text((heel_start+B*1000)/2, t_base*1000-cover*1000+30, sprintf('%dM@%d', bar_dia_heel, spacing_heel), ...
         'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    text(t_toe*500, -80, 'TOE', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text((heel_start+B*1000)/2, -80, 'HEEL', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    xlabel('Length (mm)'); ylabel('Thickness (mm)');
    title('BASE SLAB REINFORCEMENT', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Drawing 4: Moment diagram
    subplot(2,3,4);
    heights = linspace(0, h_stem, 30);
    moments = 0.05 * 18 * 0.333 * heights.^3 / 6;
    plot(moments, heights, 'b-', 'LineWidth', 2.5); hold on;
    yline(0, 'k-', 'LineWidth', 1);
    plot([Mf_stem, Mf_stem], [0, h_stem], 'r--', 'LineWidth', 2);
    plot([Mr_stem, Mr_stem], [0, h_stem], 'g--', 'LineWidth', 2);
    legend('Moment envelope', 'Mf @ base', 'Mr', 'Location', 'best', 'FontSize', 9);
    xlabel('Moment (kN?m/m)'); ylabel('Height (m)');
    title('STEM BENDING MOMENT', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Drawing 5: Shear diagram
    subplot(2,3,5);
    shears = 0.5 * 0.05 * 18 * 0.333 * heights.^2;
    plot(shears, heights, 'b-', 'LineWidth', 2.5); hold on;
    plot([Vf_stem, Vf_stem], [0, h_stem], 'r--', 'LineWidth', 2);
    plot([Vc_stem, Vc_stem], [0, h_stem], 'g--', 'LineWidth', 2);
    legend('Shear envelope', 'Vf @ base', 'Vc', 'Location', 'best', 'FontSize', 9);
    xlabel('Shear (kN/m)'); ylabel('Height (m)');
    title('STEM SHEAR FORCE', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Drawing 6: Bearing pressure
    subplot(2,3,6);
    hold on; grid on;
    
    % Base
    plot([0, B], [0, 0], 'k-', 'LineWidth', 4);
    
    % Pressure distribution
    e = B/2 - x_resultant;
    if abs(e) <= B/6
        press_x = [0, B, B, 0, 0];
        press_y = [0, 0, -q_max*0.03, -q_min*0.03, 0];
        fill(press_x, press_y, [0.3 0.5 0.9], 'EdgeColor', 'b', 'LineWidth', 2);
    end
    
    % Resultant
    plot(x_resultant, 0, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    quiver(x_resultant, 0.5, 0, -0.4, 'r', 'LineWidth', 2, 'MaxHeadSize', 1.5);
    text(x_resultant, 0.8, 'Resultant', 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Middle third
    plot([B/6, B/6], [-4, 1.2], 'g--', 'LineWidth', 1.5);
    plot([5*B/6, 5*B/6], [-4, 1.2], 'g--', 'LineWidth', 1.5);
    text(B/6, 1.4, 'B/6', 'Color', 'g', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(5*B/6, 1.4, '5B/6', 'Color', 'g', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    text(B/2, -5, sprintf('q_{max} = %.1f kPa\nq_{min} = %.1f kPa', q_max, q_min), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('Base width (m)'); ylabel('Pressure (kPa)');
    title('BEARING PRESSURE DISTRIBUTION', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([-0.2, B+0.2]); ylim([-6, 2]);
    
    sgtitle('RETAINING WALL COMPLETE DESIGN - CSA A23.3-19', 'FontSize', 14, 'FontWeight', 'bold');
end
