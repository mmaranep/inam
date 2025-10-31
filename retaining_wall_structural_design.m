%% RETAINING WALL STRUCTURAL DESIGN - CSA A23.3-19 (CANADIAN CODE)
% This script performs structural design of a cantilever retaining wall
% following CSA A23.3-19 Design of Concrete Structures
% Includes: Flexural design, Shear design, Crack control, Development length

clear; clc; close all;

%% INPUT PARAMETERS

% Wall Geometry (all dimensions in meters)
H = 6.0;              % Total wall height (m)
B = 4.0;              % Base width (m)
t_stem_top = 0.3;     % Stem thickness at top (m)
t_stem_base = 0.4;    % Stem thickness at base (m)
t_base = 0.6;         % Base thickness (m)
t_toe = 1.2;          % Toe length (m)
t_heel = B - t_toe - t_stem_base;   % Heel length (m)
cover = 0.075;        % Concrete cover (m) - 75mm for earth face

% Soil Properties
gamma_soil = 18.0;    % Unit weight of backfill soil (kN/m?)
phi = 30;             % Internal friction angle of soil (degrees)
delta = 2/3 * phi;    % Wall friction angle (degrees)
q = 10.0;             % Surcharge load on backfill (kPa)

% Foundation Soil Properties
gamma_foundation = 19.0;  % Unit weight of foundation soil (kN/m?)
phi_f = 28;               % Friction angle of foundation soil (degrees)

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
alpha_E = 1.0;        % Earth pressure factor (used with alpha_D)

% Bar diameters available (mm)
bar_sizes = [10, 15, 20, 25, 30, 35];
bar_areas = [100, 200, 300, 500, 700, 1000]; % mm?

%% PRELIMINARY CALCULATIONS

fprintf('='*70); fprintf('\n');
fprintf('RETAINING WALL STRUCTURAL DESIGN - CSA A23.3-19\n');
fprintf('='*70); fprintf('\n\n');

% Convert angles to radians
phi_rad = deg2rad(phi);
delta_rad = deg2rad(delta);
phi_f_rad = deg2rad(phi_f);

% Active earth pressure coefficient (Coulomb's theory)
Ka = cos(phi_rad)^2 / (cos(delta_rad) * (1 + sqrt(sin(phi_rad + delta_rad) * ...
     sin(phi_rad) / cos(delta_rad)))^2);

fprintf('MATERIAL PROPERTIES:\n');
fprintf('  Concrete strength (f''c) = %.1f MPa\n', fc_prime);
fprintf('  Steel yield strength (fy) = %.1f MPa\n', fy);
fprintf('  Active earth pressure coefficient (Ka) = %.4f\n\n', Ka);

%% PART 1: STEM DESIGN

fprintf('='*70); fprintf('\n');
fprintf('STEM DESIGN\n');
fprintf('='*70); fprintf('\n\n');

% Design at base of stem (critical section)
h_stem = H - t_base;  % Height of stem
t_stem = t_stem_base; % Use base thickness for design

% Factored lateral earth pressure at base
pa_soil_base = Ka * gamma_soil * h_stem;  % kPa
pa_surcharge = Ka * q;  % kPa
pa_total_base = pa_soil_base + pa_surcharge;  % kPa

% Factored pressure (ULS)
pa_factored = alpha_D * pa_total_base;  % kPa

% Factored moment at base of stem (per meter width)
M_soil = (1/6) * alpha_D * Ka * gamma_soil * h_stem^3;  % kN?m/m
M_surcharge = 0.5 * alpha_L * Ka * q * h_stem^2;  % kN?m/m
Mf_stem = M_soil + M_surcharge;  % Total factored moment

% Factored shear at base of stem
V_soil = 0.5 * alpha_D * Ka * gamma_soil * h_stem^2;  % kN/m
V_surcharge = alpha_L * Ka * q * h_stem;  % kN/m
Vf_stem = V_soil + V_surcharge;  % Total factored shear

fprintf('STEM - At base (critical section):\n');
fprintf('  Height of stem = %.2f m\n', h_stem);
fprintf('  Stem thickness = %.3f m\n', t_stem);
fprintf('  Factored pressure at base = %.2f kPa\n', pa_factored);
fprintf('  Factored moment (Mf) = %.2f kN?m/m\n', Mf_stem);
fprintf('  Factored shear (Vf) = %.2f kN/m\n\n', Vf_stem);

%% STEM FLEXURAL DESIGN (CSA A23.3)

fprintf('STEM FLEXURAL DESIGN:\n');

% Effective depth
d_stem = t_stem * 1000 - cover * 1000 - 10;  % mm (assume 20mm bar)

% Calculate required steel ratio
% Using: Mf = phi_s * As * fy * (d - a/2)
% where a = As * fy / (alpha1 * phi_c * f'c * b)
% Iterative solution

b = 1000;  % Width = 1000 mm (per meter)
alpha1 = 0.85 - 0.0015 * fc_prime;  % Stress block parameter
if alpha1 < 0.67
    alpha1 = 0.67;
end
beta1 = 0.97 - 0.0025 * fc_prime;  % Stress block parameter
if beta1 < 0.67
    beta1 = 0.67;
end

% Calculate required reinforcement
Mf_stem_Nmm = Mf_stem * 1e6;  % Convert to N?mm/m
omega = 1 - sqrt(1 - 2 * Mf_stem_Nmm / (phi_c * alpha1 * fc_prime * b * d_stem^2));
rho_required = omega * alpha1 * phi_c * fc_prime / (phi_s * fy);

% Check minimum reinforcement (CSA A23.3-19 Cl. 10.5.1.2)
rho_min = max(0.002, 1.4/fy);  % Minimum reinforcement ratio

% Check maximum reinforcement (balanced condition)
epsilon_c = 0.0035;  % Concrete strain at ultimate
epsilon_y = fy / Es;  % Steel yield strain
c_b = epsilon_c / (epsilon_c + epsilon_y) * d_stem;  % Balanced neutral axis depth
a_b = beta1 * c_b;
rho_max = alpha1 * phi_c * fc_prime * a_b / (phi_s * fy * d_stem);

rho_stem = max(rho_required, rho_min);

if rho_stem > rho_max
    fprintf('  WARNING: Required reinforcement exceeds balanced condition!\n');
    fprintf('  Consider increasing stem thickness.\n');
    rho_stem = 0.7 * rho_max;  % Limit to 70% of balanced
end

As_required_stem = rho_stem * b * d_stem;  % mm?/m

fprintf('  Effective depth (d) = %.1f mm\n', d_stem);
fprintf('  Required steel ratio (?) = %.5f\n', rho_stem);
fprintf('  Minimum steel ratio (?_min) = %.5f\n', rho_min);
fprintf('  Required steel area = %.1f mm?/m\n', As_required_stem);

% Select reinforcement
spacing_options = [100, 125, 150, 175, 200, 250, 300];
best_design = [];
min_excess = inf;

for i = 1:length(bar_sizes)
    for j = 1:length(spacing_options)
        As_provided = bar_areas(i) * 1000 / spacing_options(j);
        if As_provided >= As_required_stem
            excess = (As_provided - As_required_stem) / As_required_stem;
            if excess < min_excess && spacing_options(j) <= 300
                min_excess = excess;
                best_design = [bar_sizes(i), spacing_options(j), As_provided];
            end
        end
    end
end

if ~isempty(best_design)
    bar_dia_stem = best_design(1);
    spacing_stem = best_design(2);
    As_provided_stem = best_design(3);
    fprintf('  PROVIDED: %dM @ %d mm c/c\n', bar_dia_stem, spacing_stem);
    fprintf('  Provided steel area = %.1f mm?/m\n', As_provided_stem);
else
    fprintf('  ERROR: Cannot find suitable reinforcement arrangement!\n');
    bar_dia_stem = 20;
    spacing_stem = 150;
    As_provided_stem = 300 * 1000 / 150;
end

% Design capacity check
a = As_provided_stem * fy / (alpha1 * fc_prime * b);
Mr_stem = phi_s * As_provided_stem * fy * (d_stem - a/2) / 1e6;  % kN?m/m

fprintf('  Moment resistance (Mr) = %.2f kN?m/m\n', Mr_stem);
fprintf('  Utilization ratio = %.2f%%\n\n', Mf_stem/Mr_stem * 100);

%% STEM SHEAR DESIGN (CSA A23.3 Cl. 11)

fprintf('STEM SHEAR DESIGN:\n');

% Concrete shear resistance (without shear reinforcement)
lambda = 1.0;  % Normal density concrete
phi_c_shear = phi_c;

% Calculate beta (for members without shear reinforcement)
% Simplified method (CSA A23.3 Cl. 11.3.4)
beta = 0.21;  % Conservative simplified method

Vc = phi_c_shear * lambda * beta * sqrt(fc_prime) * b * d_stem / 1000;  % kN

fprintf('  Factored shear (Vf) = %.2f kN/m\n', Vf_stem);
fprintf('  Concrete shear resistance (Vc) = %.2f kN/m\n', Vc);

if Vf_stem <= Vc
    fprintf('  Shear reinforcement: NOT REQUIRED [OK] ?\n');
    fprintf('  Utilization ratio = %.2f%%\n\n', Vf_stem/Vc * 100);
else
    fprintf('  Shear reinforcement: REQUIRED\n');
    % Calculate required shear reinforcement
    Vs_required = Vf_stem / phi_s - Vc;
    % Av/s = Vs / (phi_s * fy * dv)
    dv = max(0.9*d_stem, 0.72*h_stem*1000);
    Av_s_required = Vs_required * 1000 / (fy * dv);  % mm?/mm
    fprintf('  Required Av/s = %.3f mm?/mm\n', Av_s_required);
    fprintf('  Recommend 10M stirrups @ 200mm c/c\n\n');
end

%% PART 2: BASE SLAB DESIGN (TOE)

fprintf('='*70); fprintf('\n');
fprintf('TOE SLAB DESIGN\n');
fprintf('='*70); fprintf('\n\n');

% Calculate bearing pressure distribution
% First, find resultant and eccentricity from stability analysis

% Weights (unfactored)
W_stem = gamma_concrete * t_stem * h_stem;
W_base = gamma_concrete * B * t_base;
W_soil_heel = gamma_soil * t_heel * h_stem;
Pa_total = 0.5 * Ka * gamma_soil * h_stem^2 + Ka * q * h_stem;
Pv = Pa_total * tan(delta_rad);
W_total = W_stem + W_base + W_soil_heel + Pv;

% Moment arms
x_stem = t_toe + t_stem/2;
x_base = B/2;
x_soil_heel = t_toe + t_stem + t_heel/2;
y_Ph = h_stem/3;

% Resisting and overturning moments
MR_total = W_stem * x_stem + W_base * x_base + W_soil_heel * x_soil_heel + Pv * B;
MO = Pa_total * y_Ph;

% Resultant location
x_resultant = (MR_total - MO) / W_total;
e = B/2 - x_resultant;

% Bearing pressure (unfactored)
if abs(e) <= B/6
    q_max_unfactored = (W_total / B) * (1 + 6*e/B);
    q_toe_unfactored = (W_total / B) * (1 - 6*e/B);
else
    L_eff = 3 * (B/2 - e);
    q_max_unfactored = 2 * W_total / L_eff;
    q_toe_unfactored = 0;
end

% Factored bearing pressure
q_max_factored = alpha_D * q_max_unfactored;
q_toe_factored = alpha_D * q_toe_unfactored;

fprintf('TOE SLAB:\n');
fprintf('  Toe length = %.2f m\n', t_toe);
fprintf('  Base thickness = %.3f m\n', t_base);
fprintf('  Bearing pressure at toe (factored) = %.2f kPa\n', q_toe_factored);
fprintf('  Bearing pressure at stem (factored) = %.2f kPa\n', q_max_factored);

% Calculate moment and shear at critical section (face of stem)
% Pressure varies linearly across toe
q_at_stem_face = q_toe_factored + (q_max_factored - q_toe_factored) * t_toe / B;

% Moment at face of stem (critical section for toe)
Mf_toe = (q_toe_factored + q_at_stem_face) / 2 * t_toe^2 / 2 - ...
         alpha_D * gamma_concrete * t_base * t_toe^2 / 2;

% Shear at distance d from face of stem
d_toe = t_base * 1000 - cover * 1000 - 10;  % mm
critical_section_toe = t_toe - d_toe/1000;
q_at_critical = q_toe_factored + (q_max_factored - q_toe_factored) * ...
                critical_section_toe / B;
Vf_toe = (q_toe_factored + q_at_critical) / 2 * critical_section_toe - ...
         alpha_D * gamma_concrete * t_base * critical_section_toe;

fprintf('  Factored moment at face of stem = %.2f kN?m/m\n', Mf_toe);
fprintf('  Factored shear at d from face = %.2f kN/m\n\n', Vf_toe);

%% TOE FLEXURAL DESIGN

fprintf('TOE FLEXURAL DESIGN:\n');

% Effective depth
fprintf('  Effective depth (d) = %.1f mm\n', d_toe);

% Calculate required reinforcement
Mf_toe_Nmm = Mf_toe * 1e6;
omega_toe = 1 - sqrt(1 - 2 * Mf_toe_Nmm / (phi_c * alpha1 * fc_prime * b * d_toe^2));
rho_required_toe = omega_toe * alpha1 * phi_c * fc_prime / (phi_s * fy);
rho_toe = max(rho_required_toe, rho_min);

As_required_toe = rho_toe * b * d_toe;

fprintf('  Required steel area = %.1f mm?/m\n', As_required_toe);

% Select reinforcement
best_design = [];
min_excess = inf;

for i = 1:length(bar_sizes)
    for j = 1:length(spacing_options)
        As_provided = bar_areas(i) * 1000 / spacing_options(j);
        if As_provided >= As_required_toe
            excess = (As_provided - As_required_toe) / As_required_toe;
            if excess < min_excess && spacing_options(j) <= 300
                min_excess = excess;
                best_design = [bar_sizes(i), spacing_options(j), As_provided];
            end
        end
    end
end

if ~isempty(best_design)
    bar_dia_toe = best_design(1);
    spacing_toe = best_design(2);
    As_provided_toe = best_design(3);
    fprintf('  PROVIDED: %dM @ %d mm c/c (bottom)\n', bar_dia_toe, spacing_toe);
    fprintf('  Provided steel area = %.1f mm?/m\n', As_provided_toe);
else
    fprintf('  ERROR: Cannot find suitable reinforcement!\n');
    bar_dia_toe = 15;
    spacing_toe = 200;
    As_provided_toe = 200 * 1000 / 200;
end

% Design capacity
a = As_provided_toe * fy / (alpha1 * fc_prime * b);
Mr_toe = phi_s * As_provided_toe * fy * (d_toe - a/2) / 1e6;

fprintf('  Moment resistance (Mr) = %.2f kN?m/m\n', Mr_toe);
fprintf('  Utilization ratio = %.2f%%\n\n', Mf_toe/Mr_toe * 100);

%% TOE SHEAR DESIGN

fprintf('TOE SHEAR DESIGN:\n');
Vc_toe = phi_c * lambda * beta * sqrt(fc_prime) * b * d_toe / 1000;

fprintf('  Factored shear (Vf) = %.2f kN/m\n', Vf_toe);
fprintf('  Concrete shear resistance (Vc) = %.2f kN/m\n', Vc_toe);

if Vf_toe <= Vc_toe
    fprintf('  Shear reinforcement: NOT REQUIRED [OK] ?\n');
    fprintf('  Utilization ratio = %.2f%%\n\n', Vf_toe/Vc_toe * 100);
else
    fprintf('  Shear reinforcement: REQUIRED [WARNING]\n\n');
end

%% PART 3: BASE SLAB DESIGN (HEEL)

fprintf('='*70); fprintf('\n');
fprintf('HEEL SLAB DESIGN\n');
fprintf('='*70); fprintf('\n\n');

fprintf('HEEL SLAB:\n');
fprintf('  Heel length = %.2f m\n', t_heel);
fprintf('  Base thickness = %.3f m\n', t_base);

% Pressure distribution under heel
q_heel_soil = q_max_factored;  % At stem face
q_heel_end = q_max_factored - (q_max_factored - q_toe_factored) * (t_toe + t_stem) / B;

fprintf('  Soil pressure at stem (factored) = %.2f kPa\n', q_heel_soil);
fprintf('  Soil pressure at heel end (factored) = %.2f kPa\n', q_heel_end);

% Calculate moment at face of stem (heel acts as cantilever)
% Weight of soil on heel
W_soil_on_heel_factored = alpha_D * gamma_soil * h_stem * t_heel;

% Moment from soil weight (creates positive moment)
M_soil_weight = W_soil_on_heel_factored * t_heel / 2;

% Moment from soil pressure (creates negative moment)
M_soil_pressure = -(q_heel_soil + q_heel_end) / 2 * t_heel^2 / 2;

% Moment from self-weight
M_self_weight = -alpha_D * gamma_concrete * t_base * t_heel^2 / 2;

Mf_heel = M_soil_weight + M_soil_pressure + M_self_weight;

% Shear at distance d from face of stem
d_heel = t_base * 1000 - cover * 1000 - 10;
critical_section_heel = t_heel - d_heel/1000;

Vf_heel = W_soil_on_heel_factored - (q_heel_soil + q_heel_end) / 2 * critical_section_heel - ...
          alpha_D * gamma_concrete * t_base * critical_section_heel;

fprintf('  Factored moment at face of stem = %.2f kN?m/m\n', abs(Mf_heel));
fprintf('  Factored shear at d from face = %.2f kN/m\n\n', abs(Vf_heel));

%% HEEL FLEXURAL DESIGN

fprintf('HEEL FLEXURAL DESIGN:\n');

% Effective depth
fprintf('  Effective depth (d) = %.1f mm\n', d_heel);

% Calculate required reinforcement
Mf_heel_Nmm = abs(Mf_heel) * 1e6;
omega_heel = 1 - sqrt(1 - 2 * Mf_heel_Nmm / (phi_c * alpha1 * fc_prime * b * d_heel^2));
rho_required_heel = omega_heel * alpha1 * phi_c * fc_prime / (phi_s * fy);
rho_heel = max(rho_required_heel, rho_min);

As_required_heel = rho_heel * b * d_heel;

fprintf('  Required steel area = %.1f mm?/m\n', As_required_heel);

% Select reinforcement
best_design = [];
min_excess = inf;

for i = 1:length(bar_sizes)
    for j = 1:length(spacing_options)
        As_provided = bar_areas(i) * 1000 / spacing_options(j);
        if As_provided >= As_required_heel
            excess = (As_provided - As_required_heel) / As_required_heel;
            if excess < min_excess && spacing_options(j) <= 300
                min_excess = excess;
                best_design = [bar_sizes(i), spacing_options(j), As_provided];
            end
        end
    end
end

if ~isempty(best_design)
    bar_dia_heel = best_design(1);
    spacing_heel = best_design(2);
    As_provided_heel = best_design(3);
    fprintf('  PROVIDED: %dM @ %d mm c/c (top)\n', bar_dia_heel, spacing_heel);
    fprintf('  Provided steel area = %.1f mm?/m\n', As_provided_heel);
else
    fprintf('  ERROR: Cannot find suitable reinforcement!\n');
    bar_dia_heel = 15;
    spacing_heel = 200;
    As_provided_heel = 200 * 1000 / 200;
end

% Design capacity
a = As_provided_heel * fy / (alpha1 * fc_prime * b);
Mr_heel = phi_s * As_provided_heel * fy * (d_heel - a/2) / 1e6;

fprintf('  Moment resistance (Mr) = %.2f kN?m/m\n', Mr_heel);
fprintf('  Utilization ratio = %.2f%%\n\n', abs(Mf_heel)/Mr_heel * 100);

%% HEEL SHEAR DESIGN

fprintf('HEEL SHEAR DESIGN:\n');
Vc_heel = phi_c * lambda * beta * sqrt(fc_prime) * b * d_heel / 1000;

fprintf('  Factored shear (Vf) = %.2f kN/m\n', abs(Vf_heel));
fprintf('  Concrete shear resistance (Vc) = %.2f kN/m\n', Vc_heel);

if abs(Vf_heel) <= Vc_heel
    fprintf('  Shear reinforcement: NOT REQUIRED [OK] ?\n');
    fprintf('  Utilization ratio = %.2f%%\n\n', abs(Vf_heel)/Vc_heel * 100);
else
    fprintf('  Shear reinforcement: REQUIRED [WARNING]\n\n');
end

%% DEVELOPMENT LENGTH (CSA A23.3 Cl. 12)

fprintf('='*70); fprintf('\n');
fprintf('DEVELOPMENT LENGTH CHECK (CSA A23.3 Cl. 12)\n');
fprintf('='*70); fprintf('\n\n');

% Basic development length (CSA A23.3 Cl. 12.2.3)
k1 = 1.0;  % Bar location factor (other than top bars)
k2 = 1.0;  % Concrete density factor (normal density)
k3 = 1.0;  % Bar size factor
k4 = 0.8;  % Bar surface factor (deformed bars)

% Development length calculation
ld_basic = @(db) k1 * k2 * k3 * k4 * fy * db / (1.0 * sqrt(fc_prime) * 1.15);

% Stem vertical bars
ld_stem = ld_basic(bar_dia_stem);
L_available_stem = t_base * 1000 - 2 * cover * 1000;  % Available length in base

fprintf('STEM VERTICAL BARS (%dM):\n', bar_dia_stem);
fprintf('  Required development length = %.0f mm\n', ld_stem);
fprintf('  Available length in base = %.0f mm\n', L_available_stem);
if L_available_stem >= ld_stem
    fprintf('  Development length: ADEQUATE [OK] ?\n\n');
else
    fprintf('  Development length: INADEQUATE [WARNING] - Hooks required\n\n');
end

% Toe bars
ld_toe = ld_basic(bar_dia_toe);
L_available_toe = (t_toe - cover) * 1000;

fprintf('TOE BARS (%dM):\n', bar_dia_toe);
fprintf('  Required development length = %.0f mm\n', ld_toe);
fprintf('  Available length = %.0f mm\n', L_available_toe);
if L_available_toe >= ld_toe
    fprintf('  Development length: ADEQUATE [OK] ?\n\n');
else
    fprintf('  Development length: INADEQUATE [WARNING] - Hooks required\n\n');
end

% Heel bars
ld_heel = ld_basic(bar_dia_heel);
L_available_heel = (t_heel - cover) * 1000;

fprintf('HEEL BARS (%dM):\n', bar_dia_heel);
fprintf('  Required development length = %.0f mm\n', ld_heel);
fprintf('  Available length = %.0f mm\n', L_available_heel);
if L_available_heel >= ld_heel
    fprintf('  Development length: ADEQUATE [OK] ?\n\n');
else
    fprintf('  Development length: INADEQUATE [WARNING] - Hooks required\n\n');
end

%% CRACK CONTROL (CSA A23.3 Cl. 10.6)

fprintf('='*70); fprintf('\n');
fprintf('CRACK CONTROL CHECK (CSA A23.3 Cl. 10.6)\n');
fprintf('='*70); fprintf('\n\n');

% Service load moments (unfactored with load factors = 1.0)
% For exposure class N or S, z ? 30000 N/mm for interior exposure

% Stem
Ms_stem = M_soil / alpha_D + M_surcharge / alpha_L;  % Service moment
fs_stem = Ms_stem * 1e6 / (As_provided_stem * (d_stem - a/3));  % Steel stress at service
dc_stem = cover * 1000 + bar_dia_stem / 2;  % Distance from tension face to center of bar
A_eff_stem = 2 * dc_stem * 1000;  % Effective tension area per bar (for 1m width)

% Number of bars in 1m
n_bars_stem = 1000 / spacing_stem;
A_per_bar_stem = A_eff_stem / n_bars_stem;

z_stem = fs_stem * (dc_stem * A_per_bar_stem)^(1/3);  % Crack control parameter

fprintf('STEM:\n');
fprintf('  Service moment = %.2f kN?m/m\n', Ms_stem);
fprintf('  Steel stress at service (fs) = %.1f MPa\n', fs_stem);
fprintf('  Crack control parameter (z) = %.0f N/mm\n', z_stem);
if z_stem <= 30000
    fprintf('  Crack control: ADEQUATE [OK] ?\n\n');
else
    fprintf('  Crack control: INADEQUATE [WARNING]\n\n');
end

%% TEMPERATURE AND SHRINKAGE REINFORCEMENT

fprintf('='*70); fprintf('\n');
fprintf('TEMPERATURE AND SHRINKAGE REINFORCEMENT\n');
fprintf('='*70); fprintf('\n\n');

% CSA A23.3 Cl. 7.8 - Minimum 0.002 of gross cross-sectional area
As_temp_stem = 0.002 * t_stem * 1000 * 1000;  % mm?/m
As_temp_base = 0.002 * t_base * 1000 * 1000;  % mm?/m

fprintf('STEM (horizontal bars):\n');
fprintf('  Minimum area required = %.1f mm?/m\n', As_temp_stem);
fprintf('  Provide: 10M @ 300 mm c/c (Area = %.1f mm?/m)\n\n', 100*1000/300);

fprintf('BASE SLAB (transverse bars):\n');
fprintf('  Minimum area required = %.1f mm?/m\n', As_temp_base);
fprintf('  Provide: 10M @ 300 mm c/c (Area = %.1f mm?/m)\n\n', 100*1000/300);

%% DESIGN SUMMARY

fprintf('='*70); fprintf('\n');
fprintf('REINFORCEMENT SUMMARY\n');
fprintf('='*70); fprintf('\n\n');

fprintf('STEM:\n');
fprintf('  Vertical (tension face):    %dM @ %d mm c/c\n', bar_dia_stem, spacing_stem);
fprintf('  Horizontal:                 10M @ 300 mm c/c\n\n');

fprintf('BASE SLAB:\n');
fprintf('  Toe (bottom):               %dM @ %d mm c/c\n', bar_dia_toe, spacing_toe);
fprintf('  Heel (top):                 %dM @ %d mm c/c\n', bar_dia_heel, spacing_heel);
fprintf('  Temperature (both faces):   10M @ 300 mm c/c\n\n');

fprintf('CONCRETE:\n');
fprintf('  Specified strength (f''c):   %.1f MPa\n', fc_prime);
fprintf('  Cover (earth face):         %.0f mm\n', cover*1000);
fprintf('  Cover (formed face):        50 mm (typical)\n\n');

%% VISUALIZATION

figure('Position', [100, 100, 1400, 800]);

% Subplot 1: Reinforcement Layout - Stem
subplot(2,3,1);
hold on; axis equal;

% Draw stem
rectangle('Position', [0, 0, t_stem*1000, h_stem*1000], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 2);

% Draw vertical reinforcement
bar_positions = cover*1000:spacing_stem:(h_stem*1000 - cover*1000);
for i = 1:length(bar_positions)
    plot(t_stem*1000 - cover*1000, bar_positions(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

% Draw horizontal reinforcement
horiz_bar_pos = cover*1000:300:(t_stem*1000 - cover*1000);
for i = 1:length(horiz_bar_pos)
    plot(horiz_bar_pos(i), h_stem*1000/2, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
end

xlabel('Width (mm)'); ylabel('Height (mm)');
title(sprintf('STEM REINFORCEMENT\nVertical: %dM @ %dmm\nHorizontal: 10M @ 300mm', bar_dia_stem, spacing_stem));
grid on;

% Subplot 2: Reinforcement Layout - Base
subplot(2,3,2);
hold on; axis equal;

% Draw base
rectangle('Position', [0, 0, B*1000, t_base*1000], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 2);

% Mark stem location
rectangle('Position', [t_toe*1000, t_base*1000, t_stem*1000, 10], 'FaceColor', [0.5 0.5 0.5]);

% Draw toe reinforcement (bottom)
toe_bar_pos = cover*1000:spacing_toe:t_toe*1000;
for i = 1:length(toe_bar_pos)
    plot(toe_bar_pos(i), cover*1000, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

% Draw heel reinforcement (top)
heel_start = (t_toe + t_stem)*1000;
heel_bar_pos = heel_start:spacing_heel:B*1000-cover*1000;
for i = 1:length(heel_bar_pos)
    plot(heel_bar_pos(i), t_base*1000-cover*1000, 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

xlabel('Length (mm)'); ylabel('Thickness (mm)');
title(sprintf('BASE SLAB REINFORCEMENT\nToe (bottom): %dM @ %dmm\nHeel (top): %dM @ %dmm', ...
      bar_dia_toe, spacing_toe, bar_dia_heel, spacing_heel));
grid on;
text(t_toe*1000/2, -50, 'TOE', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text((heel_start + B*1000)/2, -50, 'HEEL', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Subplot 3: Moment Diagram - Stem
subplot(2,3,4);
heights = linspace(0, h_stem, 50);
moments = zeros(size(heights));
for i = 1:length(heights)
    h = heights(i);
    moments(i) = (1/6) * alpha_D * Ka * gamma_soil * h^3 + 0.5 * alpha_L * Ka * q * h^2;
end
plot(moments, heights, 'b-', 'LineWidth', 2);
hold on;
plot([Mf_stem Mf_stem], [0 h_stem], 'r--', 'LineWidth', 1.5);
plot([Mr_stem Mr_stem], [0 h_stem], 'g--', 'LineWidth', 1.5);
xlabel('Moment (kN?m/m)'); ylabel('Height (m)');
title('STEM MOMENT DIAGRAM');
legend('Factored Moment', 'Mf at base', 'Mr (Resistance)', 'Location', 'best');
grid on;

% Subplot 4: Shear Diagram - Stem
subplot(2,3,5);
shears = zeros(size(heights));
for i = 1:length(heights)
    h = heights(i);
    shears(i) = 0.5 * alpha_D * Ka * gamma_soil * h^2 + alpha_L * Ka * q * h;
end
plot(shears, heights, 'b-', 'LineWidth', 2);
hold on;
plot([Vf_stem Vf_stem], [0 h_stem], 'r--', 'LineWidth', 1.5);
plot([Vc Vc], [0 h_stem], 'g--', 'LineWidth', 1.5);
xlabel('Shear (kN/m)'); ylabel('Height (m)');
title('STEM SHEAR DIAGRAM');
legend('Factored Shear', 'Vf at base', 'Vc (Resistance)', 'Location', 'best');
grid on;

% Subplot 5: Full Section View
subplot(2,3,[3,6]);
hold on; axis equal;

% Draw wall
wall_x = [t_toe, t_toe, t_toe+t_stem, t_toe+t_stem, B, B, 0, 0, t_toe];
wall_y = [0, t_base, t_base, H, H, 0, 0, t_base, 0];
fill(wall_x, wall_y, [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 2);

% Draw backfill
soil_x = [t_toe+t_stem, t_toe+t_stem, B+1, B+1, B, B, t_toe+t_stem];
soil_y = [t_base, H, H, 0, 0, t_base, t_base];
fill(soil_x, soil_y, [0.8 0.6 0.4], 'EdgeColor', 'k', 'LineWidth', 1);

% Add dimension lines
plot([0, B], [-0.3, -0.3], 'k-', 'LineWidth', 1);
text(B/2, -0.6, sprintf('B = %.2f m', B), 'HorizontalAlignment', 'center', 'FontSize', 10);
plot([t_toe, t_toe], [0, -0.2], 'k-', 'LineWidth', 1);
text(t_toe/2, -0.3, sprintf('%.2f m', t_toe), 'HorizontalAlignment', 'center', 'FontSize', 9);

% Height dimension
plot([-0.3, -0.3], [0, H], 'k-', 'LineWidth', 1);
text(-0.6, H/2, sprintf('H = %.2f m', H), 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', 10);

% Add reinforcement indicators
text(t_toe + t_stem/2, H/2, sprintf('%dM@%d', bar_dia_stem, spacing_stem), ...
     'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(t_toe/2, t_base/2, sprintf('%dM@%d', bar_dia_toe, spacing_toe), ...
     'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text((t_toe+t_stem+B)/2, t_base/2, sprintf('%dM@%d', bar_dia_heel, spacing_heel), ...
     'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

xlabel('Distance (m)'); ylabel('Height (m)');
title('COMPLETE SECTION WITH REINFORCEMENT');
xlim([-1, B+1.5]); ylim([-1, H+0.5]);
grid on;

sgtitle('RETAINING WALL STRUCTURAL DESIGN - CSA A23.3-19', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('='*70); fprintf('\n');
fprintf('STRUCTURAL DESIGN COMPLETE\n');
fprintf('='*70); fprintf('\n');
