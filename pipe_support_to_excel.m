%% PIPE SUPPORT DESIGN CALCULATION - EXCEL EXPORT
% Exports all design calculations to a formatted Excel spreadsheet
% Author: Auto-generated design tool
% Date: November 6, 2025

clear all;
clc;
fprintf('===============================================\n');
fprintf('   PIPE SUPPORT DESIGN - EXCEL EXPORT\n');
fprintf('===============================================\n\n');

% Define output filename
excel_filename = 'Pipe_Support_Design_Calculation.xlsx';
fprintf('Creating Excel file: %s\n\n', excel_filename);

%% INPUT PARAMETERS (Same as original script)

% Pipe Properties
pipe_OD = 12.75;
pipe_schedule = 10;
pipe_wall_thick = 0.165;
pipe_ID = pipe_OD - 2*pipe_wall_thick;
pipe_material_density = 7930;
water_density = 1000;
pipe_length_supported = 1000;

% Post Properties (HSS 3x3x0.25)
post_size = 3.0;
post_wall_thick = 0.25;
post_height = 1000;
post_material = 'A500 Gr.B';
post_Fy = 317;
post_Fu = 400;
post_E = 200000;

% Top Plate Properties
plate_length = 230;
plate_height = 100;
plate_thickness = 6;
plate_Fy = 250;
plate_Fu = 400;
plate_radius = pipe_OD * 25.4 / 2;

% Top Bolts
top_bolt_dia = 7/8;
top_bolt_grade = 'A325';
top_num_bolts = 2;
top_bolt_Fu = 827;

% Weld
weld_size = 6;
weld_strength = 0.6 * 485;
weld_length = 2 * (post_size * 25.4);

% Base Plate
base_plate_size = 8.0;
base_plate_thick = 0.25;
base_plate_Fy = 250;

% Anchor Bolts
anchor_dia = 0.5;
anchor_num = 4;
anchor_circle_dia = 6.0;
anchor_Fu = 400;
anchor_allowable_tension = 6.7;
anchor_allowable_shear = 7.8;

% Concrete
concrete_fc = 25;

% Safety Factors
SF_dead = 1.25;
SF_live = 1.5;
load_factor = SF_dead;

%% UNIT CONVERSIONS
mm_to_m = 0.001;
inch_to_mm = 25.4;
inch_to_m = 0.0254;
N_to_kN = 0.001;
g = 9.81;

%% SECTION PROPERTIES CALCULATION

post_outer = post_size * inch_to_mm;
post_inner = post_outer - 2 * post_wall_thick * inch_to_mm;
post_area = post_outer^2 - post_inner^2;
post_I = (post_outer^4 - post_inner^4) / 12;
post_S = post_I / (post_outer/2);
post_r = sqrt(post_I / post_area);

%% LOAD CALCULATIONS

% Pipe weight
pipe_volume_per_m = pi * ((pipe_OD*inch_to_mm/2)^2 - (pipe_ID*inch_to_mm/2)^2) * 1e-6;
pipe_weight_per_m = pipe_volume_per_m * pipe_material_density * g;
pipe_weight = pipe_weight_per_m * pipe_length_supported * mm_to_m;

% Water weight
water_volume_per_m = pi * (pipe_ID*inch_to_mm/2)^2 * 1e-6;
water_weight_per_m = water_volume_per_m * water_density * g;
water_weight = water_weight_per_m * pipe_length_supported * mm_to_m;

% Support weights
post_weight = post_area * 1e-6 * post_height * mm_to_m * 7850 * g;
plate_volume = plate_length * plate_height * plate_thickness * 1e-9;
plate_weight = plate_volume * 7850 * g;
base_plate_weight = (base_plate_size*inch_to_m)^2 * base_plate_thick*inch_to_m * 7850 * g;

% Total loads
P_dead = pipe_weight + water_weight + post_weight + plate_weight + base_plate_weight;
P_factored = load_factor * P_dead;

% Lateral loads
H_lateral = 0.05 * P_dead;
H_factored = 1.5 * H_lateral;

%% POST STRESS ANALYSIS

% Axial compression
P_compression = P_factored;
stress_axial = P_compression / post_area;

% Slenderness ratio
KL_r = post_height / post_r;
KL_r_limit = 200;

% Bending
M_lateral = H_factored * post_height;
stress_bending = M_lateral / post_S;

% Combined stress
stress_ratio = stress_axial/(0.85*post_Fy) + stress_bending/(0.66*post_Fy);

% Deflection
E_steel = post_E;
delta = (H_factored * post_height^3) / (3 * E_steel * post_I);
delta_limit = post_height / 200;

%% WELD STRESS CHECK

weld_throat = 0.707 * weld_size;
weld_effective_area = weld_throat * weld_length;
weld_stress_axial = P_factored / weld_effective_area;

eccentricity = plate_height / 2;
M_weld = H_factored * eccentricity;
weld_S_approx = weld_throat * weld_length^2 / 6;
weld_stress_moment = M_weld / weld_S_approx;

weld_stress_total = sqrt(weld_stress_axial^2 + weld_stress_moment^2);
weld_allowable = 0.3 * weld_strength;

%% TOP BOLT CHECK

top_bolt_area = pi * (top_bolt_dia * inch_to_mm)^2 / 4;
top_bolt_total_area = top_num_bolts * top_bolt_area;

bolt_shear_stress = P_factored / top_bolt_total_area;
bolt_allowable_shear = 0.4 * top_bolt_Fu;

bolt_spacing = plate_length - 50;
M_bolt = H_factored * post_height;
bolt_tension_force = M_bolt / bolt_spacing;
bolt_tension_stress = bolt_tension_force / top_bolt_area;
bolt_allowable_tension = 0.75 * top_bolt_Fu;

%% TOP PLATE STRESS CHECK

plate_contact_area = plate_length * plate_thickness;
plate_bearing_stress = P_factored / plate_contact_area;
plate_bearing_allowable = 0.9 * plate_Fy;

%% BASE PLATE AND ANCHOR BOLT CHECK

base_plate_area = (base_plate_size * inch_to_mm)^2;
concrete_bearing_stress = P_factored / base_plate_area;
concrete_bearing_allowable = 0.35 * concrete_fc;

anchor_bolt_area = pi * (anchor_dia * inch_to_mm)^2 / 4;
anchor_total_area = anchor_num * anchor_bolt_area;

anchor_shear_per_bolt = H_factored / anchor_num * N_to_kN;

moment_arm = (base_plate_size * inch_to_mm) / 2;
M_overturning = H_factored * post_height;
anchor_tension_total = M_overturning / moment_arm - P_factored;
if anchor_tension_total < 0
    anchor_tension_total = 0;
end
anchor_tension_per_bolt = anchor_tension_total / (anchor_num/2) * N_to_kN;

anchor_interaction = (anchor_shear_per_bolt/anchor_allowable_shear)^2 + ...
                     (anchor_tension_per_bolt/anchor_allowable_tension)^2;

% Base plate bending
q = concrete_bearing_stress;
cantilever = (base_plate_size*inch_to_mm - post_size*inch_to_mm) / 2;
M_plate = q * cantilever^2 / 2;
base_plate_S = (base_plate_thick*inch_to_mm)^2 / 6;
plate_bending_stress = M_plate / base_plate_S;
plate_bending_allowable = 0.75 * base_plate_Fy;

%% DELETE EXISTING FILE IF IT EXISTS
if exist(excel_filename, 'file')
    delete(excel_filename);
    fprintf('Deleted existing file.\n');
end

%% SHEET 1: INPUT PARAMETERS
fprintf('Writing Sheet 1: Input Parameters...\n');

Sheet1_Data = {
    'PIPE SUPPORT DESIGN CALCULATION', '', '', '';
    'Project:', 'Water Pipe Support', '', '';
    'Date:', datestr(now, 'dd-mmm-yyyy'), '', '';
    '', '', '', '';
    'INPUT PARAMETERS', '', '', '';
    '==================', '', '', '';
    '', '', '', '';
    'PIPE PROPERTIES', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Pipe Nominal Size', 12, 'inch', '12" Pipe';
    'Outside Diameter', pipe_OD, 'inch', 'OD';
    'Schedule', pipe_schedule, '', 'SCH 10';
    'Wall Thickness', pipe_wall_thick, 'inch', 't_pipe';
    'Inside Diameter', pipe_ID, 'inch', 'ID';
    'Material', 'SS 304L', '', 'Stainless Steel';
    'Material Density', pipe_material_density, 'kg/m³', 'ρ_steel';
    'Water Density', water_density, 'kg/m³', 'ρ_water';
    'Supported Length', pipe_length_supported, 'mm', 'L_support';
    '', '', '', '';
    'POST PROPERTIES (HSS)', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Section Size', sprintf('%.0f x %.0f', post_size, post_size), 'inch', 'Square HSS';
    'Wall Thickness', post_wall_thick, 'inch', 't_wall';
    'Height', post_height, 'mm', 'H_post';
    'Material Grade', post_material, '', 'ASTM A500';
    'Yield Strength', post_Fy, 'MPa', 'F_y';
    'Ultimate Strength', post_Fu, 'MPa', 'F_u';
    'Elastic Modulus', post_E, 'MPa', 'E';
    '', '', '', '';
    'TOP PLATE PROPERTIES', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Length', plate_length, 'mm', 'L_plate';
    'Height', plate_height, 'mm', 'H_plate';
    'Thickness', plate_thickness, 'mm', 't_plate';
    'Yield Strength', plate_Fy, 'MPa', 'F_y';
    'Radius (curved)', plate_radius, 'mm', 'R_curve';
    '', '', '', '';
    'TOP BOLTS', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Bolt Diameter', top_bolt_dia, 'inch', 'd_bolt';
    'Bolt Grade', top_bolt_grade, '', 'ASTM A325';
    'Number of Bolts', top_num_bolts, '', 'n_bolts';
    'Ultimate Strength', top_bolt_Fu, 'MPa', 'F_u';
    '', '', '', '';
    'WELD', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Weld Size', weld_size, 'mm', 'Fillet weld';
    'Weld Strength', weld_strength, 'MPa', 'E48XX electrode';
    'Weld Length', weld_length, 'mm', 'Total length';
    '', '', '', '';
    'BASE PLATE', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Plate Size', sprintf('%.0f x %.0f', base_plate_size, base_plate_size), 'inch', 'Square plate';
    'Thickness', base_plate_thick, 'inch', 't_base';
    'Yield Strength', base_plate_Fy, 'MPa', 'F_y';
    '', '', '', '';
    'ANCHOR BOLTS', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Anchor Diameter', anchor_dia, 'inch', 'd_anchor';
    'Number of Anchors', anchor_num, '', 'n_anchors';
    'Type', 'Hilti Drop-in', '', '1/2" anchors';
    'Allowable Tension', anchor_allowable_tension, 'kN', 'T_allow';
    'Allowable Shear', anchor_allowable_shear, 'kN', 'V_allow';
    '', '', '', '';
    'CONCRETE', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Compressive Strength', concrete_fc, 'MPa', 'f''c';
    '', '', '', '';
    'LOAD FACTORS', '', '', '';
    'Parameter', 'Value', 'Unit', 'Description';
    'Dead Load Factor', SF_dead, '', 'γ_D';
    'Live Load Factor', SF_live, '', 'γ_L';
};

writecell(Sheet1_Data, excel_filename, 'Sheet', 'Input Parameters');

%% SHEET 2: SECTION PROPERTIES
fprintf('Writing Sheet 2: Section Properties...\n');

Sheet2_Data = {
    'SECTION PROPERTIES', '', '', '';
    '==================', '', '', '';
    '', '', '', '';
    'POST SECTION (HSS 3x3x0.25)', '', '', '';
    'Property', 'Symbol', 'Value', 'Unit';
    'Outer Dimension', 'b_o', post_outer, 'mm';
    'Inner Dimension', 'b_i', post_inner, 'mm';
    'Cross-sectional Area', 'A', post_area, 'mm²';
    'Moment of Inertia', 'I', post_I, 'mm⁴';
    'Section Modulus', 'S', post_S, 'mm³';
    'Radius of Gyration', 'r', post_r, 'mm';
};

writecell(Sheet2_Data, excel_filename, 'Sheet', 'Section Properties');

%% SHEET 3: LOAD CALCULATIONS
fprintf('Writing Sheet 3: Load Calculations...\n');

Sheet3_Data = {
    'LOAD CALCULATIONS', '', '', '';
    '==================', '', '', '';
    '', '', '', '';
    'VERTICAL LOADS', '', '', '';
    'Load Component', 'Unit Weight', 'Total Weight (N)', 'Total Weight (kg)';
    'Pipe (empty)', sprintf('%.2f N/m', pipe_weight_per_m), pipe_weight, pipe_weight/g;
    'Water (full)', sprintf('%.2f N/m', water_weight_per_m), water_weight, water_weight/g;
    'Post', '', post_weight, post_weight/g;
    'Top Plate', '', plate_weight, plate_weight/g;
    'Base Plate', '', base_plate_weight, base_plate_weight/g;
    '', '', '', '';
    'TOTAL DEAD LOAD', '', P_dead, P_dead/g;
    'FACTORED LOAD (γ=1.25)', '', P_factored, P_factored/g;
    '', '', '', '';
    '', '', '', '';
    'LATERAL LOADS', '', '', '';
    'Load Type', 'Factor', 'Unfactored (N)', 'Factored (N)';
    'Lateral Load (5% of vertical)', '5%', H_lateral, '';
    'Factored Lateral Load', '1.5', '', H_factored;
    '', '', '', '';
    'SUMMARY', '', '', '';
    'Description', '', 'Value', 'Unit';
    'Total Vertical Dead Load', '', P_dead*N_to_kN, 'kN';
    'Total Factored Vertical Load', '', P_factored*N_to_kN, 'kN';
    'Total Factored Lateral Load', '', H_factored*N_to_kN, 'kN';
};

writecell(Sheet3_Data, excel_filename, 'Sheet', 'Load Calculations');

%% SHEET 4: POST ANALYSIS
fprintf('Writing Sheet 4: Post Analysis...\n');

Sheet4_Data = {
    'POST STRUCTURAL ANALYSIS', '', '', '', '';
    '=========================', '', '', '', '';
    '', '', '', '', '';
    '1. AXIAL COMPRESSION', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Axial Force', 'P', P_compression*N_to_kN, 'kN', '';
    'Axial Stress', 'f_a', stress_axial, 'MPa', '';
    'Allowable Stress', '0.85F_y', 0.85*post_Fy, 'MPa', '';
    'Utilization Ratio', 'f_a/F_allow', (stress_axial/(0.85*post_Fy))*100, '%', ...
        sprintf('%s', ternary(stress_axial < 0.85*post_Fy, 'OK', 'FAIL'));
    '', '', '', '', '';
    '2. SLENDERNESS CHECK', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Effective Length Factor', 'K', 1.0, '', 'Pinned-Fixed';
    'Unbraced Length', 'L', post_height, 'mm', '';
    'Radius of Gyration', 'r', post_r, 'mm', '';
    'Slenderness Ratio', 'KL/r', KL_r, '', '';
    'Limit', '', KL_r_limit, '', '';
    'Check', '', '', '', sprintf('%s', ternary(KL_r < KL_r_limit, 'OK', 'FAIL'));
    '', '', '', '', '';
    '3. BENDING ANALYSIS', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Lateral Force', 'H', H_factored*N_to_kN, 'kN', '';
    'Bending Moment', 'M', M_lateral*1e-6, 'kN·m', '';
    'Bending Stress', 'f_b', stress_bending, 'MPa', '';
    'Allowable Stress', '0.66F_y', 0.66*post_Fy, 'MPa', '';
    'Utilization Ratio', 'f_b/F_allow', (stress_bending/(0.66*post_Fy))*100, '%', ...
        sprintf('%s', ternary(stress_bending < 0.66*post_Fy, 'OK', 'FAIL'));
    '', '', '', '', '';
    '4. COMBINED STRESS', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Interaction Equation', 'f_a/F_a + f_b/F_b', stress_ratio, '', '';
    'Limit', '', 1.0, '', '';
    'Check', '', '', '', sprintf('%s', ternary(stress_ratio <= 1.0, 'OK', 'FAIL'));
    '', '', '', '', '';
    '5. DEFLECTION', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Lateral Deflection', 'δ', delta, 'mm', '';
    'Limit', 'H/200', delta_limit, 'mm', '';
    'Ratio', 'δ/δ_limit', (delta/delta_limit)*100, '%', ...
        sprintf('%s', ternary(delta < delta_limit, 'OK', 'FAIL'));
};

writecell(Sheet4_Data, excel_filename, 'Sheet', 'Post Analysis');

%% SHEET 5: WELD CHECK
fprintf('Writing Sheet 5: Weld Check...\n');

Sheet5_Data = {
    'WELD STRESS CHECK (Post to Plate)', '', '', '', '';
    '===================================', '', '', '', '';
    '', '', '', '', '';
    'WELD PROPERTIES', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Notes';
    'Weld Size', 'w', weld_size, 'mm', 'Fillet weld';
    'Throat Thickness', 't_e', weld_throat, 'mm', '0.707 × w';
    'Weld Length', 'L_w', weld_length, 'mm', 'Perimeter';
    'Effective Area', 'A_e', weld_effective_area, 'mm²', 't_e × L_w';
    '', '', '', '', '';
    'STRESS ANALYSIS', '', '', '', '';
    'Load Type', 'Symbol', 'Stress (MPa)', 'Allowable (MPa)', 'Status';
    'Axial Stress', 'f_w,axial', weld_stress_axial, '', '';
    'Moment Stress', 'f_w,moment', weld_stress_moment, '', '';
    'Combined Stress', 'f_w,total', weld_stress_total, weld_allowable, ...
        sprintf('%s', ternary(weld_stress_total < weld_allowable, 'OK', 'FAIL'));
    '', '', '', '', '';
    'UTILIZATION', '', '', '', '';
    'Check', '', 'Value', 'Unit', 'Status';
    'Weld Stress Ratio', '', (weld_stress_total/weld_allowable)*100, '%', ...
        sprintf('%s', ternary(weld_stress_total < weld_allowable, 'PASS', 'FAIL'));
};

writecell(Sheet5_Data, excel_filename, 'Sheet', 'Weld Check');

%% SHEET 6: TOP BOLT CHECK
fprintf('Writing Sheet 6: Top Bolt Check...\n');

Sheet6_Data = {
    'TOP BOLT CHECK (Pipe Flange Connection)', '', '', '', '';
    '========================================', '', '', '', '';
    '', '', '', '', '';
    'BOLT PROPERTIES', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Notes';
    'Bolt Diameter', 'd', top_bolt_dia, 'inch', sprintf('%.2f mm', top_bolt_dia*inch_to_mm);
    'Bolt Grade', '', top_bolt_grade, '', 'ASTM A325';
    'Number of Bolts', 'n', top_num_bolts, '', '1 each side';
    'Bolt Area (each)', 'A_b', top_bolt_area, 'mm²', '';
    'Total Bolt Area', 'A_total', top_bolt_total_area, 'mm²', '';
    '', '', '', '', '';
    'SHEAR CHECK', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Vertical Load', 'P', P_factored*N_to_kN, 'kN', '';
    'Shear Stress', 'f_v', bolt_shear_stress, 'MPa', '';
    'Allowable Shear', '0.4F_u', bolt_allowable_shear, 'MPa', '';
    'Utilization', '', (bolt_shear_stress/bolt_allowable_shear)*100, '%', ...
        sprintf('%s', ternary(bolt_shear_stress < bolt_allowable_shear, 'OK', 'FAIL'));
    '', '', '', '', '';
    'TENSION CHECK', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Moment', 'M', M_bolt*1e-6, 'kN·m', '';
    'Bolt Spacing', 'e', bolt_spacing, 'mm', '';
    'Tension Force/bolt', 'T', bolt_tension_force, 'N', '';
    'Tension Stress', 'f_t', bolt_tension_stress, 'MPa', '';
    'Allowable Tension', '0.75F_u', bolt_allowable_tension, 'MPa', '';
    'Utilization', '', (bolt_tension_stress/bolt_allowable_tension)*100, '%', ...
        sprintf('%s', ternary(bolt_tension_stress < bolt_allowable_tension, 'OK', 'FAIL'));
};

writecell(Sheet6_Data, excel_filename, 'Sheet', 'Top Bolt Check');

%% SHEET 7: TOP PLATE CHECK
fprintf('Writing Sheet 7: Top Plate Check...\n');

Sheet7_Data = {
    'TOP PLATE STRESS CHECK', '', '', '', '';
    '======================', '', '', '', '';
    '', '', '', '', '';
    'PLATE PROPERTIES', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Notes';
    'Length', 'L', plate_length, 'mm', '';
    'Height', 'H', plate_height, 'mm', '';
    'Thickness', 't', plate_thickness, 'mm', '';
    'Yield Strength', 'F_y', plate_Fy, 'MPa', '';
    '', '', '', '', '';
    'BEARING STRESS CHECK', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Applied Load', 'P', P_factored*N_to_kN, 'kN', '';
    'Contact Area', 'A_c', plate_contact_area, 'mm²', 'Approximate';
    'Bearing Stress', 'f_p', plate_bearing_stress, 'MPa', '';
    'Allowable Bearing', '0.9F_y', plate_bearing_allowable, 'MPa', '';
    'Utilization', '', (plate_bearing_stress/plate_bearing_allowable)*100, '%', ...
        sprintf('%s', ternary(plate_bearing_stress < plate_bearing_allowable, 'OK', 'FAIL'));
};

writecell(Sheet7_Data, excel_filename, 'Sheet', 'Top Plate Check');

%% SHEET 8: BASE PLATE & ANCHORS
fprintf('Writing Sheet 8: Base Plate & Anchors...\n');

Sheet8_Data = {
    'BASE PLATE AND ANCHOR BOLT CHECK', '', '', '', '';
    '=================================', '', '', '', '';
    '', '', '', '', '';
    'CONCRETE BEARING', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Base Plate Size', '', sprintf('%.0f x %.0f', base_plate_size, base_plate_size), 'inch', '';
    'Base Plate Area', 'A_bp', base_plate_area, 'mm²', '';
    'Applied Load', 'P', P_factored*N_to_kN, 'kN', '';
    'Bearing Stress', 'f_p', concrete_bearing_stress, 'MPa', '';
    'Allowable (0.35f''c)', 'f_allow', concrete_bearing_allowable, 'MPa', '';
    'Utilization', '', (concrete_bearing_stress/concrete_bearing_allowable)*100, '%', ...
        sprintf('%s', ternary(concrete_bearing_stress < concrete_bearing_allowable, 'OK', 'FAIL'));
    '', '', '', '', '';
    'ANCHOR BOLT SHEAR', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Lateral Force', 'H', H_factored*N_to_kN, 'kN', '';
    'Number of Anchors', 'n', anchor_num, '', '';
    'Shear per Bolt', 'V_bolt', anchor_shear_per_bolt, 'kN', '';
    'Allowable Shear', 'V_allow', anchor_allowable_shear, 'kN', '';
    'Utilization', '', (anchor_shear_per_bolt/anchor_allowable_shear)*100, '%', ...
        sprintf('%s', ternary(anchor_shear_per_bolt < anchor_allowable_shear, 'OK', 'FAIL'));
    '', '', '', '', '';
    'ANCHOR BOLT TENSION', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Overturning Moment', 'M_ot', M_overturning*1e-6, 'kN·m', '';
    'Moment Arm', 'e', moment_arm, 'mm', '';
    'Tension per Bolt', 'T_bolt', anchor_tension_per_bolt, 'kN', '';
    'Allowable Tension', 'T_allow', anchor_allowable_tension, 'kN', '';
    'Utilization', '', (anchor_tension_per_bolt/anchor_allowable_tension)*100, '%', ...
        sprintf('%s', ternary(anchor_tension_per_bolt < anchor_allowable_tension, 'OK', 'FAIL'));
    '', '', '', '', '';
    'COMBINED INTERACTION', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Interaction Ratio', '(V/V_a)² + (T/T_a)²', anchor_interaction, '', '';
    'Limit', '', 1.0, '', '';
    'Check', '', '', '', sprintf('%s', ternary(anchor_interaction <= 1.0, 'OK', 'FAIL'));
    '', '', '', '', '';
    'BASE PLATE BENDING', '', '', '', '';
    'Parameter', 'Symbol', 'Value', 'Unit', 'Status';
    'Cantilever Length', 'c', cantilever, 'mm', '';
    'Bearing Pressure', 'q', concrete_bearing_stress, 'MPa', '';
    'Bending Stress', 'f_b', plate_bending_stress, 'MPa', '';
    'Allowable (0.75F_y)', 'f_allow', plate_bending_allowable, 'MPa', '';
    'Utilization', '', (plate_bending_stress/plate_bending_allowable)*100, '%', ...
        sprintf('%s', ternary(plate_bending_stress < plate_bending_allowable, 'OK', 'FAIL'));
};

writecell(Sheet8_Data, excel_filename, 'Sheet', 'Base Plate & Anchors');

%% SHEET 9: SUMMARY
fprintf('Writing Sheet 9: Summary...\n');

% Determine overall status
all_checks = [
    stress_axial < 0.85*post_Fy;
    KL_r < KL_r_limit;
    stress_bending < 0.66*post_Fy;
    stress_ratio <= 1.0;
    delta < delta_limit;
    weld_stress_total < weld_allowable;
    bolt_shear_stress < bolt_allowable_shear;
    bolt_tension_stress < bolt_allowable_tension;
    plate_bearing_stress < plate_bearing_allowable;
    concrete_bearing_stress < concrete_bearing_allowable;
    anchor_shear_per_bolt < anchor_allowable_shear;
    anchor_tension_per_bolt < anchor_allowable_tension;
    anchor_interaction <= 1.0;
    plate_bending_stress < plate_bending_allowable;
];

overall_status = all(all_checks);

Sheet9_Data = {
    'DESIGN SUMMARY', '', '', '';
    '==============', '', '', '';
    '', '', '', '';
    'PROJECT INFORMATION', '', '', '';
    'Description', '', 'Value', '';
    'Project Name', '', 'Water Pipe Support Design', '';
    'Calculation Date', '', datestr(now, 'dd-mmm-yyyy'), '';
    'Support Type', '', 'Single Post with Curved Top Plate', '';
    'Post Height', '', sprintf('%.0f mm', post_height), '';
    '', '', '', '';
    'LOADING SUMMARY', '', '', '';
    'Load Description', '', 'Value', 'Unit';
    'Total Dead Load', '', sprintf('%.2f', P_dead*N_to_kN), 'kN';
    'Factored Vertical Load', '', sprintf('%.2f', P_factored*N_to_kN), 'kN';
    'Factored Lateral Load', '', sprintf('%.2f', H_factored*N_to_kN), 'kN';
    '', '', '', '';
    'COMPONENT UTILIZATION', '', '', '';
    'Component', 'Check Type', 'Utilization (%)', 'Status';
    'Post', 'Axial Stress', sprintf('%.1f', (stress_axial/(0.85*post_Fy))*100), ...
        sprintf('%s', ternary(stress_axial < 0.85*post_Fy, 'OK', 'FAIL'));
    'Post', 'Bending Stress', sprintf('%.1f', (stress_bending/(0.66*post_Fy))*100), ...
        sprintf('%s', ternary(stress_bending < 0.66*post_Fy, 'OK', 'FAIL'));
    'Post', 'Combined Stress', sprintf('%.1f', stress_ratio*100), ...
        sprintf('%s', ternary(stress_ratio <= 1.0, 'OK', 'FAIL'));
    'Post', 'Deflection', sprintf('%.1f', (delta/delta_limit)*100), ...
        sprintf('%s', ternary(delta < delta_limit, 'OK', 'FAIL'));
    'Weld (6mm)', 'Combined Stress', sprintf('%.1f', (weld_stress_total/weld_allowable)*100), ...
        sprintf('%s', ternary(weld_stress_total < weld_allowable, 'OK', 'FAIL'));
    'Top Bolts (7/8")', 'Shear', sprintf('%.1f', (bolt_shear_stress/bolt_allowable_shear)*100), ...
        sprintf('%s', ternary(bolt_shear_stress < bolt_allowable_shear, 'OK', 'FAIL'));
    'Top Bolts (7/8")', 'Tension', sprintf('%.1f', (bolt_tension_stress/bolt_allowable_tension)*100), ...
        sprintf('%s', ternary(bolt_tension_stress < bolt_allowable_tension, 'OK', 'FAIL'));
    'Top Plate', 'Bearing', sprintf('%.1f', (plate_bearing_stress/plate_bearing_allowable)*100), ...
        sprintf('%s', ternary(plate_bearing_stress < plate_bearing_allowable, 'OK', 'FAIL'));
    'Base Plate', 'Concrete Bearing', sprintf('%.1f', (concrete_bearing_stress/concrete_bearing_allowable)*100), ...
        sprintf('%s', ternary(concrete_bearing_stress < concrete_bearing_allowable, 'OK', 'FAIL'));
    'Base Plate', 'Plate Bending', sprintf('%.1f', (plate_bending_stress/plate_bending_allowable)*100), ...
        sprintf('%s', ternary(plate_bending_stress < plate_bending_allowable, 'OK', 'FAIL'));
    'Anchors (1/2")', 'Shear', sprintf('%.1f', (anchor_shear_per_bolt/anchor_allowable_shear)*100), ...
        sprintf('%s', ternary(anchor_shear_per_bolt < anchor_allowable_shear, 'OK', 'FAIL'));
    'Anchors (1/2")', 'Tension', sprintf('%.1f', (anchor_tension_per_bolt/anchor_allowable_tension)*100), ...
        sprintf('%s', ternary(anchor_tension_per_bolt < anchor_allowable_tension, 'OK', 'FAIL'));
    'Anchors (1/2")', 'Interaction', sprintf('%.1f', anchor_interaction*100), ...
        sprintf('%s', ternary(anchor_interaction <= 1.0, 'OK', 'FAIL'));
    '', '', '', '';
    '', '', '', '';
    'OVERALL DESIGN STATUS', '', sprintf('%s', ternary(overall_status, 'ACCEPTABLE', 'REQUIRES REVISION')), '';
    '', '', '', '';
    'NOTES AND RECOMMENDATIONS', '', '', '';
    '1. All calculations based on AISC/CSA design standards', '', '', '';
    '2. Factored loads used for strength design', '', '', '';
    '3. Lateral load assumed as 5% of vertical load', '', '', '';
    '4. Consider dynamic effects and impact factors if applicable', '', '', '';
    '5. Verify anchor embedment depth and edge distances', '', '', '';
    '6. Professional engineer review required before construction', '', '', '';
    '7. Field verification of dimensions and conditions required', '', '', '';
};

writecell(Sheet9_Data, excel_filename, 'Sheet', 'Summary');

%% COMPLETE
fprintf('\n===============================================\n');
fprintf('   Excel file created successfully!\n');
fprintf('   Filename: %s\n', excel_filename);
fprintf('===============================================\n\n');

fprintf('The Excel file contains 9 sheets:\n');
fprintf('  1. Input Parameters\n');
fprintf('  2. Section Properties\n');
fprintf('  3. Load Calculations\n');
fprintf('  4. Post Analysis\n');
fprintf('  5. Weld Check\n');
fprintf('  6. Top Bolt Check\n');
fprintf('  7. Top Plate Check\n');
fprintf('  8. Base Plate & Anchors\n');
fprintf('  9. Summary\n\n');

fprintf('Overall Design Status: %s\n\n', ternary(overall_status, 'ACCEPTABLE ✓', 'REQUIRES REVISION ✗'));

%% HELPER FUNCTION
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
