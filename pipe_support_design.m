%% PIPE SUPPORT DESIGN CALCULATION
% Design of a single post pipe support with curved plate connection
% Author: Auto-generated design tool
% Date: November 6, 2025

clear all;
clc;
fprintf('===============================================\n');
fprintf('   PIPE SUPPORT DESIGN CALCULATION\n');
fprintf('===============================================\n\n');

%% INPUT PARAMETERS

% Pipe Properties
pipe_OD = 12.75;           % 12" pipe outside diameter (inches)
pipe_schedule = 10;         % Schedule 10
pipe_wall_thick = 0.165;    % Wall thickness for 12" SCH 10 (inches)
pipe_ID = pipe_OD - 2*pipe_wall_thick;  % Inside diameter
pipe_material_density = 7930;  % SS 304L density (kg/m³)
water_density = 1000;       % Water density (kg/m³)
pipe_length_supported = 1000;  % Assumed span length (mm) - can be adjusted

% Post Properties (HSS 3x3x0.25)
post_size = 3.0;            % 3 inch square tube (inches)
post_wall_thick = 0.25;     % Wall thickness (inches)
post_height = 1000;         % Height from floor (mm)
post_material = 'A500 Gr.B'; % Steel grade
post_Fy = 317;             % Yield strength (MPa) for HSS A500 Gr.B
post_Fu = 400;             % Ultimate strength (MPa)
post_E = 200000;           % Elastic modulus (MPa)

% Top Plate Properties
plate_length = 230;         % mm
plate_height = 100;         % mm
plate_thickness = 6;        % mm
plate_Fy = 250;            % Yield strength (MPa) for plate steel
plate_Fu = 400;            % Ultimate strength (MPa)
plate_radius = pipe_OD * 25.4 / 2;  % Curved to fit pipe (mm)

% Top Bolts (Pipe Flange Connection)
top_bolt_dia = 7/8;         % 7/8 inch diameter bolts
top_bolt_grade = 'A325';    % Bolt grade
top_num_bolts = 2;          % 2 bolts (one each side)
top_bolt_Fu = 827;          % Ultimate strength (MPa) for A325

% Weld (Post to Plate)
weld_size = 6;              % 6mm fillet weld
weld_strength = 0.6 * 485;  % Weld electrode strength (MPa), E48XX
weld_length = 2 * (post_size * 25.4);  % Weld around post perimeter (mm)

% Base Plate Properties
base_plate_size = 8.0;      % 8 inch x 8 inch (inches)
base_plate_thick = 0.25;    % 0.25 inch (inches)
base_plate_Fy = 250;        % Yield strength (MPa)

% Anchor Bolts
anchor_dia = 0.5;           % 0.5 inch Hilti drop-in anchors
anchor_num = 4;             % 4 anchor bolts
anchor_circle_dia = 6.0;    % Bolt circle diameter (inches)
anchor_Fu = 400;            % Ultimate strength (MPa) for anchors
anchor_allowable_tension = 6.7;  % kN (typical for 1/2" Hilti drop-in)
anchor_allowable_shear = 7.8;    % kN (typical for 1/2" Hilti drop-in)

% Concrete Properties
concrete_fc = 25;           % Concrete strength (MPa)

% Safety Factors
SF_dead = 1.25;             % Dead load factor
SF_live = 1.5;              % Live load factor (if applicable)
load_factor = SF_dead;      % Using dead load only

%% UNIT CONVERSIONS
mm_to_m = 0.001;
inch_to_mm = 25.4;
inch_to_m = 0.0254;
N_to_kN = 0.001;
g = 9.81;  % Gravity (m/s²)

%% SECTION PROPERTIES CALCULATION

fprintf('1. SECTION PROPERTIES\n');
fprintf('   -------------------\n');

% Post Section Properties (HSS 3x3x0.25)
post_outer = post_size * inch_to_mm;  % mm
post_inner = post_outer - 2 * post_wall_thick * inch_to_mm;  % mm
post_area = post_outer^2 - post_inner^2;  % mm²
post_I = (post_outer^4 - post_inner^4) / 12;  % mm⁴ (moment of inertia)
post_S = post_I / (post_outer/2);  % mm³ (section modulus)
post_r = sqrt(post_I / post_area);  % mm (radius of gyration)

fprintf('   Post (HSS 3x3x0.25):\n');
fprintf('   - Cross-sectional area: %.2f mm²\n', post_area);
fprintf('   - Moment of inertia: %.2f x 10⁴ mm⁴\n', post_I/1e4);
fprintf('   - Section modulus: %.2f x 10³ mm³\n', post_S/1e3);
fprintf('   - Radius of gyration: %.2f mm\n\n', post_r);

%% LOAD CALCULATIONS

fprintf('2. LOAD CALCULATIONS\n');
fprintf('   -----------------\n');

% Pipe weight (per unit length)
pipe_volume_per_m = pi * ((pipe_OD*inch_to_mm/2)^2 - (pipe_ID*inch_to_mm/2)^2) * 1e-6;  % m³/m
pipe_weight_per_m = pipe_volume_per_m * pipe_material_density * g;  % N/m
pipe_weight = pipe_weight_per_m * pipe_length_supported * mm_to_m;  % N

fprintf('   Pipe (12" SCH 10 SS304L):\n');
fprintf('   - Outside diameter: %.2f mm\n', pipe_OD*inch_to_mm);
fprintf('   - Wall thickness: %.2f mm\n', pipe_wall_thick*inch_to_mm);
fprintf('   - Pipe weight/meter: %.2f N/m\n', pipe_weight_per_m);
fprintf('   - Total pipe weight: %.2f N (%.2f kg)\n', pipe_weight, pipe_weight/g);

% Water weight
water_volume_per_m = pi * (pipe_ID*inch_to_mm/2)^2 * 1e-6;  % m³/m
water_weight_per_m = water_volume_per_m * water_density * g;  % N/m
water_weight = water_weight_per_m * pipe_length_supported * mm_to_m;  % N

fprintf('   - Water weight/meter: %.2f N/m\n', water_weight_per_m);
fprintf('   - Total water weight: %.2f N (%.2f kg)\n', water_weight, water_weight/g);

% Support component weights
post_weight = post_area * 1e-6 * post_height * mm_to_m * 7850 * g;  % N
plate_volume = plate_length * plate_height * plate_thickness * 1e-9;  % m³
plate_weight = plate_volume * 7850 * g;  % N
base_plate_weight = (base_plate_size*inch_to_m)^2 * base_plate_thick*inch_to_m * 7850 * g;  % N

fprintf('   - Post weight: %.2f N (%.2f kg)\n', post_weight, post_weight/g);
fprintf('   - Top plate weight: %.2f N (%.2f kg)\n', plate_weight, plate_weight/g);
fprintf('   - Base plate weight: %.2f N (%.2f kg)\n\n', base_plate_weight, base_plate_weight/g);

% Total loads
P_dead = pipe_weight + water_weight + post_weight + plate_weight + base_plate_weight;  % N
P_factored = load_factor * P_dead;  % N

fprintf('   TOTAL LOADS:\n');
fprintf('   - Dead load: %.2f N (%.2f kg)\n', P_dead, P_dead/g);
fprintf('   - Factored load (SF=%.2f): %.2f N (%.2f kg)\n\n', load_factor, P_factored, P_factored/g);

%% ASSUME LATERAL LOADS (Wind/Seismic - 5% of vertical load as example)
H_lateral = 0.05 * P_dead;  % Lateral load (N)
H_factored = 1.5 * H_lateral;  % Factored lateral load

fprintf('   Assumed lateral load (5%% of vertical): %.2f N\n', H_lateral);
fprintf('   Factored lateral load: %.2f N\n\n', H_factored);

%% POST STRESS ANALYSIS

fprintf('3. POST STRUCTURAL ANALYSIS\n');
fprintf('   ------------------------\n');

% Axial compression
P_compression = P_factored;  % N
stress_axial = P_compression / post_area;  % MPa

fprintf('   Axial Compression:\n');
fprintf('   - Axial force: %.2f kN\n', P_compression * N_to_kN);
fprintf('   - Axial stress: %.2f MPa\n', stress_axial);
fprintf('   - Allowable stress: %.2f MPa (0.85*Fy)\n', 0.85*post_Fy);
fprintf('   - Utilization ratio: %.2f%%\n', (stress_axial/(0.85*post_Fy))*100);

% Check slenderness ratio
KL_r = post_height / post_r;  % Slenderness ratio (K=1 for pinned-fixed)
KL_r_limit = 200;  % Typical limit for compression members

fprintf('   - Slenderness ratio (KL/r): %.2f\n', KL_r);
fprintf('   - Limit: %.2f\n', KL_r_limit);
if KL_r < KL_r_limit
    fprintf('   - Status: OK ✓\n\n');
else
    fprintf('   - Status: FAIL ✗\n\n');
end

% Bending moment due to lateral load
M_lateral = H_factored * post_height;  % N·mm
stress_bending = M_lateral / post_S;  % MPa

fprintf('   Bending (Due to Lateral Load):\n');
fprintf('   - Bending moment: %.2f kN·m\n', M_lateral * 1e-6);
fprintf('   - Bending stress: %.2f MPa\n', stress_bending);
fprintf('   - Allowable stress: %.2f MPa (0.66*Fy)\n', 0.66*post_Fy);
fprintf('   - Utilization ratio: %.2f%%\n', (stress_bending/(0.66*post_Fy))*100);

% Combined stress (interaction equation)
stress_ratio = stress_axial/(0.85*post_Fy) + stress_bending/(0.66*post_Fy);
fprintf('   - Combined stress ratio: %.3f\n', stress_ratio);
if stress_ratio <= 1.0
    fprintf('   - Status: OK ✓\n\n');
else
    fprintf('   - Status: FAIL ✗\n\n');
end

% Deflection
E_steel = post_E;  % MPa
delta = (H_factored * post_height^3) / (3 * E_steel * post_I);  % mm
delta_limit = post_height / 200;  % Deflection limit (H/200)

fprintf('   Deflection:\n');
fprintf('   - Lateral deflection: %.2f mm\n', delta);
fprintf('   - Limit (H/200): %.2f mm\n', delta_limit);
if delta < delta_limit
    fprintf('   - Status: OK ✓\n\n');
else
    fprintf('   - Status: FAIL ✗\n\n');
end

%% WELD STRESS CHECK (POST TO PLATE)

fprintf('4. WELD STRESS CHECK (Post to Plate)\n');
fprintf('   ----------------------------------\n');

% Effective throat thickness
weld_throat = 0.707 * weld_size;  % mm
weld_effective_area = weld_throat * weld_length;  % mm²

% Shear stress in weld (from axial load)
weld_stress_axial = P_factored / weld_effective_area;  % MPa

% Stress from moment (if any eccentricity)
eccentricity = plate_height / 2;  % mm (approximate)
M_weld = H_factored * eccentricity;  % N·mm
% Section modulus of weld group (approximate)
weld_S_approx = weld_throat * weld_length^2 / 6;  % mm³
weld_stress_moment = M_weld / weld_S_approx;  % MPa

% Combined weld stress
weld_stress_total = sqrt(weld_stress_axial^2 + weld_stress_moment^2);  % MPa
weld_allowable = 0.3 * weld_strength;  % Allowable weld stress

fprintf('   Weld size: %d mm fillet\n', weld_size);
fprintf('   Throat thickness: %.2f mm\n', weld_throat);
fprintf('   Effective weld length: %.2f mm\n', weld_length);
fprintf('   Weld stress (axial): %.2f MPa\n', weld_stress_axial);
fprintf('   Weld stress (moment): %.2f MPa\n', weld_stress_moment);
fprintf('   Combined weld stress: %.2f MPa\n', weld_stress_total);
fprintf('   Allowable weld stress: %.2f MPa\n', weld_allowable);
fprintf('   Utilization ratio: %.2f%%\n', (weld_stress_total/weld_allowable)*100);
if weld_stress_total < weld_allowable
    fprintf('   Status: OK ✓\n\n');
else
    fprintf('   Status: FAIL ✗\n\n');
end

%% TOP BOLT CHECK (PIPE FLANGE CONNECTION)

fprintf('5. TOP BOLT CHECK (Pipe Flange Connection)\n');
fprintf('   ---------------------------------------\n');

% Bolt properties
top_bolt_area = pi * (top_bolt_dia * inch_to_mm)^2 / 4;  % mm² per bolt
top_bolt_total_area = top_num_bolts * top_bolt_area;  % mm²

% Shear stress (vertical load)
bolt_shear_stress = P_factored / top_bolt_total_area;  % MPa
bolt_allowable_shear = 0.4 * top_bolt_Fu;  % MPa

% Tension due to moment (if lateral load creates tension)
bolt_spacing = plate_length - 50;  % mm (approximate, 25mm edge distance each side)
M_bolt = H_factored * post_height;  % N·mm
bolt_tension_force = M_bolt / bolt_spacing;  % N (per bolt on tension side)
bolt_tension_stress = bolt_tension_force / top_bolt_area;  % MPa
bolt_allowable_tension = 0.75 * top_bolt_Fu;  % MPa

fprintf('   Bolt diameter: %.3f inch (%.2f mm)\n', top_bolt_dia, top_bolt_dia*inch_to_mm);
fprintf('   Number of bolts: %d\n', top_num_bolts);
fprintf('   Bolt area per bolt: %.2f mm²\n', top_bolt_area);
fprintf('\n   Shear Check:\n');
fprintf('   - Shear stress: %.2f MPa\n', bolt_shear_stress);
fprintf('   - Allowable shear stress: %.2f MPa\n', bolt_allowable_shear);
fprintf('   - Utilization ratio: %.2f%%\n', (bolt_shear_stress/bolt_allowable_shear)*100);
if bolt_shear_stress < bolt_allowable_shear
    fprintf('   - Status: OK ✓\n');
else
    fprintf('   - Status: FAIL ✗\n');
end

fprintf('\n   Tension Check (due to lateral moment):\n');
fprintf('   - Tension force per bolt: %.2f N\n', bolt_tension_force);
fprintf('   - Tension stress: %.2f MPa\n', bolt_tension_stress);
fprintf('   - Allowable tension stress: %.2f MPa\n', bolt_allowable_tension);
fprintf('   - Utilization ratio: %.2f%%\n', (bolt_tension_stress/bolt_allowable_tension)*100);
if bolt_tension_stress < bolt_allowable_tension
    fprintf('   - Status: OK ✓\n\n');
else
    fprintf('   - Status: FAIL ✗\n\n');
end

%% TOP PLATE STRESS CHECK

fprintf('6. TOP PLATE STRESS CHECK\n');
fprintf('   -----------------------\n');

% Bearing stress under pipe
plate_contact_area = plate_length * plate_thickness;  % mm² (approximate)
plate_bearing_stress = P_factored / plate_contact_area;  % MPa
plate_bearing_allowable = 0.9 * plate_Fy;  % MPa

fprintf('   Plate dimensions: %d x %d x %d mm\n', plate_length, plate_height, plate_thickness);
fprintf('   Bearing stress: %.2f MPa\n', plate_bearing_stress);
fprintf('   Allowable bearing stress: %.2f MPa\n', plate_bearing_allowable);
fprintf('   Utilization ratio: %.2f%%\n', (plate_bearing_stress/plate_bearing_allowable)*100);
if plate_bearing_stress < plate_bearing_allowable
    fprintf('   Status: OK ✓\n\n');
else
    fprintf('   Status: FAIL ✗\n\n');
end

%% BASE PLATE AND ANCHOR BOLT CHECK

fprintf('7. BASE PLATE AND ANCHOR BOLT CHECK\n');
fprintf('   ---------------------------------\n');

% Base plate bearing on concrete
base_plate_area = (base_plate_size * inch_to_mm)^2;  % mm²
concrete_bearing_stress = P_factored / base_plate_area;  % MPa
concrete_bearing_allowable = 0.35 * concrete_fc;  % MPa (ACI 318)

fprintf('   Base plate: %.1f x %.1f x %.2f inch\n', base_plate_size, base_plate_size, base_plate_thick);
fprintf('   Base plate area: %.2f mm²\n', base_plate_area);
fprintf('   Concrete bearing stress: %.2f MPa\n', concrete_bearing_stress);
fprintf('   Allowable bearing stress: %.2f MPa (0.35*f\'c)\n', concrete_bearing_allowable);
fprintf('   Utilization ratio: %.2f%%\n', (concrete_bearing_stress/concrete_bearing_allowable)*100);
if concrete_bearing_stress < concrete_bearing_allowable
    fprintf('   Status: OK ✓\n\n');
else
    fprintf('   Status: FAIL ✗\n\n');
end

% Anchor bolt check
anchor_bolt_area = pi * (anchor_dia * inch_to_mm)^2 / 4;  % mm² per bolt
anchor_total_area = anchor_num * anchor_bolt_area;  % mm²

% Shear force per anchor
anchor_shear_per_bolt = H_factored / anchor_num * N_to_kN;  % kN

% Tension force due to overturning moment
moment_arm = (base_plate_size * inch_to_mm) / 2;  % mm
M_overturning = H_factored * post_height;  % N·mm
% Tension in anchors (on tension side)
anchor_tension_total = M_overturning / moment_arm - P_factored;  % N
if anchor_tension_total < 0
    anchor_tension_total = 0;  % No tension if compression dominates
end
anchor_tension_per_bolt = anchor_tension_total / (anchor_num/2) * N_to_kN;  % kN (2 bolts on tension side)

fprintf('   Anchor Bolts:\n');
fprintf('   - Number: %d\n', anchor_num);
fprintf('   - Diameter: %.1f inch (%.2f mm)\n', anchor_dia, anchor_dia*inch_to_mm);
fprintf('   - Type: Hilti drop-in anchors\n\n');

fprintf('   Shear Check:\n');
fprintf('   - Shear per bolt: %.2f kN\n', anchor_shear_per_bolt);
fprintf('   - Allowable shear: %.2f kN\n', anchor_allowable_shear);
fprintf('   - Utilization ratio: %.2f%%\n', (anchor_shear_per_bolt/anchor_allowable_shear)*100);
if anchor_shear_per_bolt < anchor_allowable_shear
    fprintf('   - Status: OK ✓\n');
else
    fprintf('   - Status: FAIL ✗\n');
end

fprintf('\n   Tension Check:\n');
fprintf('   - Tension per bolt: %.2f kN\n', anchor_tension_per_bolt);
fprintf('   - Allowable tension: %.2f kN\n', anchor_allowable_tension);
fprintf('   - Utilization ratio: %.2f%%\n', (anchor_tension_per_bolt/anchor_allowable_tension)*100);
if anchor_tension_per_bolt < anchor_allowable_tension
    fprintf('   - Status: OK ✓\n');
else
    fprintf('   - Status: FAIL ✗\n');
end

% Combined interaction (shear-tension)
anchor_interaction = (anchor_shear_per_bolt/anchor_allowable_shear)^2 + ...
                     (anchor_tension_per_bolt/anchor_allowable_tension)^2;
fprintf('\n   Combined interaction: %.3f\n', anchor_interaction);
if anchor_interaction <= 1.0
    fprintf('   Status: OK ✓\n\n');
else
    fprintf('   Status: FAIL ✗\n\n');
end

% Base plate bending check
% Maximum moment in base plate
q = concrete_bearing_stress;  % MPa (uniform pressure)
cantilever = (base_plate_size*inch_to_mm - post_size*inch_to_mm) / 2;  % mm
M_plate = q * cantilever^2 / 2;  % N·mm/mm
base_plate_S = (base_plate_thick*inch_to_mm)^2 / 6;  % mm³/mm
plate_bending_stress = M_plate / base_plate_S;  % MPa
plate_bending_allowable = 0.75 * base_plate_Fy;  % MPa

fprintf('   Base Plate Bending:\n');
fprintf('   - Bending stress: %.2f MPa\n', plate_bending_stress);
fprintf('   - Allowable bending stress: %.2f MPa\n', plate_bending_allowable);
fprintf('   - Utilization ratio: %.2f%%\n', (plate_bending_stress/plate_bending_allowable)*100);
if plate_bending_stress < plate_bending_allowable
    fprintf('   - Status: OK ✓\n\n');
else
    fprintf('   - Status: FAIL ✗\n\n');
end

%% SUMMARY

fprintf('===============================================\n');
fprintf('   DESIGN SUMMARY\n');
fprintf('===============================================\n\n');

fprintf('Total Factored Load: %.2f kN (%.2f kg)\n', P_factored*N_to_kN, P_factored/g);
fprintf('Total Lateral Load: %.2f kN\n\n', H_factored*N_to_kN);

fprintf('COMPONENT CHECKS:\n');
fprintf('1. Post (HSS 3x3x0.25):\n');
fprintf('   - Axial stress: %.1f%% utilized\n', (stress_axial/(0.85*post_Fy))*100);
fprintf('   - Bending stress: %.1f%% utilized\n', (stress_bending/(0.66*post_Fy))*100);
fprintf('   - Combined: %.1f%% utilized\n', stress_ratio*100);
fprintf('   - Deflection: %.1f%% of limit\n\n', (delta/delta_limit)*100);

fprintf('2. Welds (6mm fillet):\n');
fprintf('   - Weld stress: %.1f%% utilized\n\n', (weld_stress_total/weld_allowable)*100);

fprintf('3. Top Bolts (7/8" dia):\n');
fprintf('   - Shear: %.1f%% utilized\n', (bolt_shear_stress/bolt_allowable_shear)*100);
fprintf('   - Tension: %.1f%% utilized\n\n', (bolt_tension_stress/bolt_allowable_tension)*100);

fprintf('4. Top Plate (230x100x6mm):\n');
fprintf('   - Bearing: %.1f%% utilized\n\n', (plate_bearing_stress/plate_bearing_allowable)*100);

fprintf('5. Base Plate (8x8x0.25"):\n');
fprintf('   - Concrete bearing: %.1f%% utilized\n', (concrete_bearing_stress/concrete_bearing_allowable)*100);
fprintf('   - Plate bending: %.1f%% utilized\n\n', (plate_bending_stress/plate_bending_allowable)*100);

fprintf('6. Anchor Bolts (1/2" Hilti):\n');
fprintf('   - Shear: %.1f%% utilized\n', (anchor_shear_per_bolt/anchor_allowable_shear)*100);
fprintf('   - Tension: %.1f%% utilized\n', (anchor_tension_per_bolt/anchor_allowable_tension)*100);
fprintf('   - Interaction: %.1f%% utilized\n\n', anchor_interaction*100);

fprintf('===============================================\n');
fprintf('   DESIGN CALCULATION COMPLETE\n');
fprintf('===============================================\n');

% Overall design status
all_checks_pass = true;
if stress_ratio > 1.0 || delta > delta_limit || ...
   weld_stress_total > weld_allowable || ...
   bolt_shear_stress > bolt_allowable_shear || ...
   bolt_tension_stress > bolt_allowable_tension || ...
   plate_bearing_stress > plate_bearing_allowable || ...
   concrete_bearing_stress > concrete_bearing_allowable || ...
   anchor_interaction > 1.0 || ...
   plate_bending_stress > plate_bending_allowable
    all_checks_pass = false;
end

fprintf('\nOVERALL DESIGN STATUS: ');
if all_checks_pass
    fprintf('ACCEPTABLE ✓\n');
else
    fprintf('REQUIRES REVISION ✗\n');
end

fprintf('\nNote: This calculation assumes standard conditions.\n');
fprintf('Verify all assumptions and code requirements before construction.\n');
fprintf('Consider dynamic loads, impact factors, and safety factors per local codes.\n');
