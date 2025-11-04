%% WATER_RETAINING_WALL_ANALYSIS
% Stability assessment for a water-retaining wall section in accordance with
% CSA practice. The script evaluates overturning, sliding, and bearing-capacity
% criteria and produces diagnostic visualisations.
%
% All forces are computed per metre length of wall. Dimensions are provided in
% millimetres in the problem statement and converted to metres below. Several
% secondary geometric properties (e.g., stem thickness, base slab thickness)
% were not specified; reasonable default assumptions are introduced as
% parameters so they can be adjusted easily.

clear; clc; close all;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS
% --------------------------------------------------------------------------

% Geometry (converted to metres)
stem_height               = 3050e-3;   % m
toe_width                 = 610e-3;    % m (stem centreline to toe edge)
heel_width                = 1270e-3;   % m (stem centreline to heel edge)
stem_thickness            = 300e-3;    % m (assumed constant thickness)
base_slab_thickness       = 550e-3;    % m (assumed)
base_shear_key_width      = 280e-3;    % m
base_shear_key_depth      = 370e-3;    % m
gutter_horizontal_width   = 800e-3;    % m
gutter_horizontal_thk     = 150e-3;    % m
gutter_vertical_height    = 290e-3;    % m
gutter_vertical_thk       = 150e-3;    % m
gutter_offset_from_top    = 400e-3;    % m (distance below stem crest)
roof_slab_thickness       = 225e-3;    % m
roof_slab_tributary_span  = 2.75;      % m (width carried by wall)

% Loading
water_depth_heel          = 1600e-3;   % m (measured from heel top)

% Material properties
gamma_concrete            = 24.0;      % kN/m^3
gamma_water               = 9.81;      % kN/m^3
gamma_soil                = 19.0;      % kN/m^3
phi_soil                  = deg2rad(32);  % radians
delta_interface           = deg2rad(15);  % radians

% Design criteria (CSA practice)
FS_overturning_req        = 1.50;
FS_sliding_req            = 1.50;
FS_bearing_req            = 3.00;

%% -------------------------------------------------------------------------
%  DERIVED GEOMETRY AND COMPONENT WEIGHTS
% --------------------------------------------------------------------------

base_width        = toe_width + heel_width;   % total footing width
stem_centerline   = toe_width;                % distance from toe to stem CL

components = struct('Name', {}, 'Weight', {}, 'LeverArm', {}, 'Moment', {});

% Base slab (excluding shear key)
volume_base_slab  = base_width * base_slab_thickness;  % m^3 per metre length
weight_base_slab  = gamma_concrete * volume_base_slab;  % kN/m
x_base_slab       = base_width / 2;
components(end+1) = makeComponent("Base slab", weight_base_slab, x_base_slab);

% Shear key (additional volume below base)
volume_shear_key  = base_shear_key_width * base_shear_key_depth;
weight_shear_key  = gamma_concrete * volume_shear_key;
x_shear_key       = stem_centerline;
components(end+1) = makeComponent("Shear key", weight_shear_key, x_shear_key);

% Stem (assumed prismatic)
volume_stem       = stem_thickness * stem_height;
weight_stem       = gamma_concrete * volume_stem;
x_stem            = stem_centerline;
components(end+1) = makeComponent("Stem", weight_stem, x_stem);

% Concrete gutter: horizontal slab
volume_gutter_h   = gutter_horizontal_width * gutter_horizontal_thk;
weight_gutter_h   = gamma_concrete * volume_gutter_h;
x_gutter_h        = stem_centerline + gutter_horizontal_width / 2;
components(end+1) = makeComponent("Gutter (horizontal)", weight_gutter_h, x_gutter_h);

% Concrete gutter: vertical curb
volume_gutter_v   = gutter_vertical_height * gutter_vertical_thk;
weight_gutter_v   = gamma_concrete * volume_gutter_v;
x_gutter_v        = stem_centerline + gutter_horizontal_width - gutter_vertical_thk / 2;
components(end+1) = makeComponent("Gutter (vertical)", weight_gutter_v, x_gutter_v);

% Roof slab tributary load
volume_roof_load  = roof_slab_thickness * roof_slab_tributary_span;
weight_roof_load  = gamma_concrete * volume_roof_load;
x_roof            = stem_centerline;
components(end+1) = makeComponent("Roof slab load", weight_roof_load, x_roof);

% Summary arrays
weights      = [components.Weight];
lever_arms   = [components.LeverArm];
moments      = [components.Moment];

W_total      = sum(weights);
M_resisting  = sum(moments);

%% -------------------------------------------------------------------------
%  HYDROSTATIC ACTIONS AND UPLIFT
% --------------------------------------------------------------------------

% Hydrostatic thrust on the heel face of the stem (triangular distribution)
hydro_thrust       = 0.5 * gamma_water * water_depth_heel^2;  % kN/m
hydro_resultant_z  = base_slab_thickness + water_depth_heel / 3; % m (from toe pivot)
M_overturning      = hydro_thrust * hydro_resultant_z;

% Uplift pressure under the footing (triangular: zero at toe, max at heel)
uplift_pressure_max = gamma_water * water_depth_heel;        % kPa = kN/m^2
uplift_force        = 0.5 * uplift_pressure_max * base_width; % kN/m
uplift_lever_arm    = 2 * base_width / 3;                     % centroid from toe (m)
M_uplift            = uplift_force * uplift_lever_arm;        % kN·m/m (destabilising)

% Net vertical action for sliding and bearing
V_net = W_total - uplift_force;                               % kN/m (downward positive)

%% -------------------------------------------------------------------------
%  STABILITY CHECKS
% --------------------------------------------------------------------------

% Overturning factor of safety
M_resisting_net = M_resisting - M_uplift;                     % kN·m/m
FS_overturning  = M_resisting_net / M_overturning;

% Sliding factor of safety
friction_resistance = V_net * tan(delta_interface);
Kp                  = tan(pi/4 + phi_soil/2)^2;
passive_resistance  = 0.5 * Kp * gamma_soil * base_shear_key_depth^2 * base_shear_key_width;
H_driving           = hydro_thrust;
FS_sliding          = (friction_resistance + passive_resistance) / H_driving;

% Resultant location and bearing pressures
M_net_toe      = M_resisting - M_uplift - M_overturning;
x_resultant    = M_net_toe / V_net;                     % m from toe
eccentricity   = base_width/2 - x_resultant;            % positive = toward toe
q_avg          = V_net / base_width;                    % kPa (kN/m^2)
q_max          = q_avg * (1 + 6 * eccentricity / base_width);
q_min          = q_avg * (1 - 6 * eccentricity / base_width);

% Bearing capacity (Terzaghi, strip footing, c=0, Df≈0)
Nq     = exp(pi * tan(phi_soil)) * tan(pi/4 + phi_soil/2)^2;
Ngamma = 2 * (Nq + 1) * tan(phi_soil);
q_ult  = 0.5 * gamma_soil * base_width * Ngamma;         % kPa
q_allow = q_ult / FS_bearing_req;
FS_bearing = q_ult / q_max;

%% -------------------------------------------------------------------------
%  OUTPUT SUMMARY
% --------------------------------------------------------------------------

components_table = table( ...
    string({components.Name})', ...
    weights', ...
    lever_arms', ...
    moments', ...
    'VariableNames', {'Component', 'Weight_kN_per_m', 'LeverArm_m', 'MomentAboutToe_kNm_per_m'});

fprintf('\nCSA Water-Retaining Wall Stability Check\n');
fprintf('---------------------------------------\n');
fprintf('Base width                : %.3f m\n', base_width);
fprintf('Assumed base thickness    : %.0f mm\n', base_slab_thickness * 1e3);
fprintf('Assumed stem thickness    : %.0f mm\n', stem_thickness * 1e3);
fprintf('Water depth (heel)        : %.2f m\n', water_depth_heel);
fprintf('Total vertical weight     : %.2f kN/m\n', W_total);
fprintf('Net vertical load (W-U)   : %.2f kN/m\n', V_net);
fprintf('Hydrostatic thrust        : %.2f kN/m\n', hydro_thrust);
fprintf('Uplift force (net heel)   : %.2f kN/m\n\n', uplift_force);

disp(components_table);

% Factor of safety summary
metrics     = ["Overturning"; "Sliding"; "Bearing"];
FS_values   = [FS_overturning; FS_sliding; FS_bearing];
FS_required = [FS_overturning_req; FS_sliding_req; FS_bearing_req];
status      = arrayfun(@(fs,req) passFail(fs >= req), FS_values, FS_required, 'UniformOutput', false);

results_table = table(metrics, FS_values, FS_required, string(status), ...
    'VariableNames', {'Check', 'FS_Computed', 'FS_Required', 'Status'});

fprintf('Stability Factors of Safety\n');
fprintf('---------------------------\n');
disp(results_table);

fprintf('Bearing pressures (toe critical):\n');
fprintf('  q_max = %.2f kPa\n', q_max);
fprintf('  q_min = %.2f kPa\n', q_min);
fprintf('  Eccentricity = %.3f m (limit = %.3f m)\n', eccentricity, base_width/6);
fprintf('  Allowable bearing (FS=%.1f) = %.2f kPa\n', FS_bearing_req, q_allow);

if FS_overturning < FS_overturning_req
    fprintf('WARNING: Overturning FS below CSA target (%.2f < %.2f).\n', FS_overturning, FS_overturning_req);
end
if FS_sliding < FS_sliding_req
    fprintf('WARNING: Sliding FS below CSA target (%.2f < %.2f).\n', FS_sliding, FS_sliding_req);
end
if FS_bearing < FS_bearing_req
    fprintf('WARNING: Bearing capacity FS below CSA target (%.2f < %.2f).\n', FS_bearing, FS_bearing_req);
end
if q_min < 0
    fprintf('WARNING: Heel in tension (q_{min} < 0). Consider increasing base width or weight.\n');
end

%% -------------------------------------------------------------------------
%  VISUALISATIONS
% --------------------------------------------------------------------------

pressure_scale = 0.015; % horizontal scale for hydrostatic pressure diagram (m per kPa)

figure('Name', 'Wall Geometry and Stability Summary', 'Color', 'w');
tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Geometry plot -----------------------------------------------------------
nexttile(tl, 1);
hold on; axis equal;

% Base slab polygon
base_poly = [0, 0; base_width, 0; base_width, base_slab_thickness; 0, base_slab_thickness];
patch(base_poly(:,1), base_poly(:,2), [0.85 0.85 0.88], 'EdgeColor', 'k');

% Shear key polygon
key_left  = stem_centerline - base_shear_key_width/2;
key_right = stem_centerline + base_shear_key_width/2;
key_poly  = [key_left, 0; key_right, 0; key_right, -base_shear_key_depth; key_left, -base_shear_key_depth];
patch(key_poly(:,1), key_poly(:,2), [0.75 0.75 0.78], 'EdgeColor', 'k');

% Stem polygon
stem_left  = stem_centerline - stem_thickness/2;
stem_right = stem_centerline + stem_thickness/2;
stem_poly  = [stem_left, base_slab_thickness; stem_right, base_slab_thickness; ...
              stem_right, base_slab_thickness + stem_height; stem_left, base_slab_thickness + stem_height];
patch(stem_poly(:,1), stem_poly(:,2), [0.90 0.90 0.94], 'EdgeColor', 'k');

% Gutter horizontal slab
gutter_top    = base_slab_thickness + stem_height - gutter_offset_from_top;
gutter_bottom = gutter_top - gutter_horizontal_thk;
gutter_left   = stem_centerline;
gutter_right  = stem_centerline + gutter_horizontal_width;
gutter_h_poly = [gutter_left, gutter_bottom; gutter_right, gutter_bottom; ...
                  gutter_right, gutter_top; gutter_left, gutter_top];
patch(gutter_h_poly(:,1), gutter_h_poly(:,2), [0.80 0.88 0.90], 'EdgeColor', 'k');

% Gutter vertical curb
curb_left   = gutter_right - gutter_vertical_thk;
curb_right  = gutter_right;
curb_top    = gutter_top + gutter_vertical_height;
curb_poly   = [curb_left, gutter_top; curb_right, gutter_top; curb_right, curb_top; curb_left, curb_top];
patch(curb_poly(:,1), curb_poly(:,2), [0.78 0.86 0.88], 'EdgeColor', 'k');

% Water body
water_poly = [base_width, base_slab_thickness; base_width, base_slab_thickness + water_depth_heel; ...
              base_width + 0.25, base_slab_thickness + water_depth_heel; base_width + 0.25, base_slab_thickness];
patch(water_poly(:,1), water_poly(:,2), [0.70 0.80 0.95], 'FaceAlpha', 0.6, 'EdgeColor', [0.4 0.4 0.6]);

% Hydrostatic pressure diagram (triangle)
pressure_at_base = gamma_water * water_depth_heel; % kPa
pressure_poly = [stem_right, base_slab_thickness; ...
                 stem_right + pressure_at_base * pressure_scale, base_slab_thickness; ...
                 stem_right, base_slab_thickness + water_depth_heel];
patch(pressure_poly(:,1), pressure_poly(:,2), [0.40 0.55 0.90], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot([stem_right, stem_right], [base_slab_thickness, base_slab_thickness + water_depth_heel], 'b', 'LineWidth', 1.2);

% Resultant location marker
plot([x_resultant x_resultant], [0, base_slab_thickness], 'r--', 'LineWidth', 1.2);

text(base_width*0.02, base_slab_thickness + stem_height + 0.15, 'Wall Elevation', 'FontWeight', 'bold');
text(stem_right + pressure_at_base * pressure_scale + 0.05, base_slab_thickness + water_depth_heel/2, ...
    sprintf('p_{max}=%.1f kPa', pressure_at_base), 'Color', [0.1 0.2 0.6]);
text(x_resultant, -0.08, 'Resultant', 'Color', 'r', 'HorizontalAlignment', 'center');

ylim([-base_shear_key_depth-0.1, base_slab_thickness + stem_height + 0.4]);
xlim([-0.15, base_width + 0.4]);
xlabel('Horizontal distance (m)');
ylabel('Elevation (m)');
title('Geometry & Hydrostatic Actions');
grid on;

% Factors of safety plot --------------------------------------------------
nexttile(tl, 2);
hold on; grid on;

bar_data = [FS_values FS_required];
hb = bar(1:3, bar_data, 'grouped');
hb(1).FaceColor = [0.55 0.75 0.55];
hb(2).FaceColor = [0.80 0.60 0.60];

for k = 1:numel(hb)
    x_points = hb(k).XEndPoints;
    y_vals   = hb(k).YData;
    labels   = strings(size(y_vals));
    switch k
        case 1
            labels = compose('%.2f', y_vals);
            text_color = [0 0 0];
            font_weight = 'bold';
        case 2
            labels = compose('Req %.2f', y_vals);
            text_color = [0.6 0.1 0.1];
            font_weight = 'normal';
    end
    text(x_points, y_vals + 0.05, labels, ...
        'HorizontalAlignment', 'center', 'Color', text_color, 'FontWeight', font_weight);
end

set(gca, 'XTick', 1:3, 'XTickLabel', metrics);
xlim([0.5, 3.5]);
ylabel('Factor of Safety');
title('Factor of Safety Comparison');

legend({'Computed FS', 'CSA Target'}, 'Location', 'northwest');

%% -------------------------------------------------------------------------
%  LOCAL FUNCTIONS
% --------------------------------------------------------------------------

function component = makeComponent(name, weight, leverArm)
% makeComponent Creates a struct describing a gravity load component.
    component.Name    = name;
    component.Weight  = weight;
    component.LeverArm = leverArm;
    component.Moment  = weight * leverArm;
end

function status = passFail(condition)
% passFail Returns PASS/FAIL style labels for tabulated output.
    if condition
        status = "PASS";
    else
        status = "CHECK";
    end
end

