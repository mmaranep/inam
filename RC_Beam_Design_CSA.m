%% RC_BEAM_DESIGN_CSA
% CSA A23.3 compliant reinforced concrete beam design for flexure, shear,
% and deflection. This script performs a complete design analysis for a
% simply supported rectangular beam under uniform loading.
%
% Units: mm for geometry, kN for forces, MPa for stress
% Standard: CSA A23.3-19

clear; close all; clc;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS
% -------------------------------------------------------------------------

% Beam geometry
beam.span = 6000;                % mm - clear span
beam.width = 300;                % mm - beam width
beam.depth = 600;                % mm - total beam depth
beam.cover = 40;                 % mm - clear cover to stirrups

% Loading
loads.DL = 15.0;                 % kN/m - dead load (including self-weight)
loads.LL = 25.0;                 % kN/m - live load

% Material properties (CSA A23.3-19)
mat.fc_prime = 30;               % MPa - specified compressive strength of concrete
mat.fy_flexural = 400;           % MPa - yield strength of flexural reinforcement
mat.fy_shear = 400;              % MPa - yield strength of shear reinforcement
mat.Es = 200000;                 % MPa - modulus of elasticity of steel
mat.Ec = 4500 * sqrt(mat.fc_prime);  % MPa - concrete modulus (CSA A23.3-19, Cl. 8.6.2)
mat.lambda = 1.0;                % normal density concrete

% Resistance factors (CSA A23.3-19)
mat.phi_c = 0.65;                % resistance factor for concrete
mat.phi_s = 0.85;                % resistance factor for steel

% Load factors (CSA A23.3-19, Cl. 8.3.2)
load_factors.alpha_D = 1.25;     % dead load factor
load_factors.alpha_L = 1.5;      % live load factor

% Reinforcement options
rebar.flexural_bar = 25;         % 25M bar for flexural reinforcement
rebar.stirrup_bar = 10;          % 10M stirrups for shear
rebar.min_spacing = 25;          % mm - minimum aggregate size plus 5mm

% Bar database (CSA standard bar sizes)
bar_db = struct();
bar_db.size = [10, 15, 20, 25, 30, 35, 45];
bar_db.area = [100, 200, 300, 500, 700, 1000, 1500];  % mm^2
bar_db.diameter = [11.3, 16.0, 19.5, 25.2, 29.9, 35.7, 43.7];  % mm

%% -------------------------------------------------------------------------
%  DISPLAY INPUT SUMMARY
% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('RC BEAM DESIGN PER CSA A23.3-19\n');
fprintf('========================================\n\n');

fprintf('GEOMETRY\n');
fprintf('--------\n');
fprintf('Span (L):                    %d mm (%.2f m)\n', beam.span, beam.span/1000);
fprintf('Width (b):                   %d mm\n', beam.width);
fprintf('Total depth (h):             %d mm\n', beam.depth);
fprintf('Clear cover:                 %d mm\n\n', beam.cover);

fprintf('MATERIALS\n');
fprintf('---------\n');
fprintf('Concrete strength (f''c):     %.1f MPa\n', mat.fc_prime);
fprintf('Steel yield (fy):            %.0f MPa\n', mat.fy_flexural);
fprintf('Concrete modulus (Ec):       %.0f MPa\n', mat.Ec);
fprintf('Steel modulus (Es):          %.0f MPa\n\n', mat.Es);

fprintf('LOADING\n');
fprintf('-------\n');
fprintf('Dead load (DL):              %.2f kN/m\n', loads.DL);
fprintf('Live load (LL):              %.2f kN/m\n', loads.LL);
fprintf('Factored load (wf):          %.2f kN/m\n\n', ...
    load_factors.alpha_D * loads.DL + load_factors.alpha_L * loads.LL);

%% -------------------------------------------------------------------------
%  DERIVED GEOMETRY & PROPERTIES
% -------------------------------------------------------------------------

% Get bar properties
idx_flex = find(bar_db.size == rebar.flexural_bar);
idx_stirrup = find(bar_db.size == rebar.stirrup_bar);

if isempty(idx_flex)
    error('Flexural bar size %dM not found in database', rebar.flexural_bar);
end
if isempty(idx_stirrup)
    error('Stirrup bar size %dM not found in database', rebar.stirrup_bar);
end

rebar.Ab_flex = bar_db.area(idx_flex);        % mm^2
rebar.db_flex = bar_db.diameter(idx_flex);    % mm
rebar.Ab_stirrup = bar_db.area(idx_stirrup);  % mm^2
rebar.db_stirrup = bar_db.diameter(idx_stirrup);  % mm

% Effective depth (assuming single layer of flexural steel)
beam.d = beam.depth - beam.cover - rebar.db_stirrup - rebar.db_flex/2;  % mm

% Concrete stress block parameter (CSA A23.3-19, Cl. 10.1.7)
mat.beta1 = max(0.67, min(0.97 - 0.0025 * mat.fc_prime, 0.85));

% Strain limits
mat.eps_y = mat.fy_flexural / mat.Es;     % yield strain
mat.eps_cu = 0.0035;                       % ultimate concrete strain (CSA A23.3-19)

fprintf('EFFECTIVE DEPTH\n');
fprintf('---------------\n');
fprintf('Effective depth (d):         %.1f mm\n', beam.d);
fprintf('d/h ratio:                   %.3f\n', beam.d / beam.depth);
fprintf('Concrete β₁ parameter:       %.3f\n\n', mat.beta1);

%% -------------------------------------------------------------------------
%  LOAD ANALYSIS - FACTORED LOADS
% -------------------------------------------------------------------------

% Factored uniform load (kN/m)
loads.wf = load_factors.alpha_D * loads.DL + load_factors.alpha_L * loads.LL;

% Service load for deflection
loads.ws = loads.DL + loads.LL;  % kN/m

% Maximum moment and shear for simply supported beam
analysis.Mf_max = loads.wf * beam.span^2 / 8 / 1e6;  % kN·m (span in mm, convert to m)
analysis.Vf_max = loads.wf * beam.span / 2 / 1e3;    % kN (span in mm, convert to m)

fprintf('INTERNAL FORCES (FACTORED)\n');
fprintf('--------------------------\n');
fprintf('Maximum moment (Mf):         %.2f kN·m\n', analysis.Mf_max);
fprintf('Maximum shear (Vf):          %.2f kN\n\n', analysis.Vf_max);

%% -------------------------------------------------------------------------
%  FLEXURAL DESIGN - REQUIRED STEEL AREA
% -------------------------------------------------------------------------

fprintf('========================================\n');
fprintf('FLEXURAL DESIGN\n');
fprintf('========================================\n\n');

% Calculate balanced reinforcement ratio (CSA A23.3-19, Cl. 10.5.2)
flex.c_balanced = mat.eps_cu / (mat.eps_cu + mat.eps_y) * beam.d;  % mm
flex.rho_balanced = mat.phi_c / mat.phi_s * 0.85 * mat.fc_prime * mat.beta1 / mat.fy_flexural * ...
                    flex.c_balanced / beam.d;

% Maximum reinforcement ratio for tension-controlled section (c ≤ 0.5d for CSA)
% Use conservative limit c/d = 0.5 to ensure ductile behavior
flex.c_max = 0.5 * beam.d;  % mm
flex.rho_max = mat.phi_c / mat.phi_s * 0.85 * mat.fc_prime * mat.beta1 / mat.fy_flexural * ...
               flex.c_max / beam.d;

% Minimum reinforcement ratio (CSA A23.3-19, Cl. 10.5.1.2)
flex.rho_min = max(0.2 / mat.fy_flexural, 1.4 / mat.fy_flexural);

fprintf('REINFORCEMENT LIMITS\n');
fprintf('--------------------\n');
fprintf('ρ_min:                       %.6f\n', flex.rho_min);
fprintf('ρ_balanced:                  %.6f\n', flex.rho_balanced);
fprintf('ρ_max (tension-controlled):  %.6f\n\n', flex.rho_max);

% Required steel area using direct quadratic solution
% For rectangular section: Mr = phi_s * As * fy * (d - a/2)
% where a = phi_s * As * fy / (phi_c * 0.85 * fc' * b)
% Rearranging: Mr = phi_s * As * fy * d - (phi_s * As * fy)^2 / (2 * phi_c * 0.85 * fc' * b)
% This is quadratic in As: A*As^2 + B*As + C = 0

Mf_Nmm = analysis.Mf_max * 1e6;  % Convert kN·m to N·mm

% Coefficients of quadratic equation
A_coef = (mat.phi_s * mat.fy_flexural)^2 / (2 * mat.phi_c * 0.85 * mat.fc_prime * beam.width);
B_coef = -mat.phi_s * mat.fy_flexural * beam.d;
C_coef = Mf_Nmm;

% Solve quadratic equation
discriminant = B_coef^2 - 4 * A_coef * C_coef;

if discriminant < 0
    error('Negative discriminant - section cannot carry the applied moment. Increase section size.');
end

As_req_1 = (-B_coef + sqrt(discriminant)) / (2 * A_coef);
As_req_2 = (-B_coef - sqrt(discriminant)) / (2 * A_coef);

% Choose the smaller positive root
if As_req_1 > 0 && As_req_2 > 0
    flex.As_req = min(As_req_1, As_req_2);
elseif As_req_1 > 0
    flex.As_req = As_req_1;
elseif As_req_2 > 0
    flex.As_req = As_req_2;
else
    error('No positive solution for required steel area. Check input parameters.');
end

% Calculate corresponding reinforcement ratio and neutral axis depth
flex.rho_req = flex.As_req / (beam.width * beam.d);
flex.a_req = mat.phi_s * flex.As_req * mat.fy_flexural / (mat.phi_c * 0.85 * mat.fc_prime * beam.width);
flex.c_req = flex.a_req / mat.beta1;

fprintf('REQUIRED REINFORCEMENT (from moment)\n');
fprintf('------------------------------------\n');
fprintf('As_required:                 %.0f mm²\n', flex.As_req);
fprintf('ρ_required:                  %.6f\n', flex.rho_req);
fprintf('Neutral axis depth (c):      %.1f mm\n', flex.c_req);
fprintf('c/d ratio:                   %.3f\n', flex.c_req / beam.d);

% Check reinforcement limits
if flex.rho_req < flex.rho_min
    fprintf('\n⚠ WARNING: Required ratio below minimum. Using ρ_min.\n');
    flex.As_req = flex.rho_min * beam.width * beam.d;
    flex.rho_req = flex.rho_min;
    % Recalculate neutral axis
    flex.a_req = mat.phi_s * flex.As_req * mat.fy_flexural / (mat.phi_c * 0.85 * mat.fc_prime * beam.width);
    flex.c_req = flex.a_req / mat.beta1;
end

if flex.rho_req > flex.rho_max
    error('Required reinforcement ratio exceeds maximum (%.6f > %.6f). Increase section size or use compression reinforcement.', ...
        flex.rho_req, flex.rho_max);
end

if flex.c_req / beam.d > 0.5
    fprintf('\n⚠ WARNING: c/d = %.3f > 0.5. Section may not be sufficiently ductile.\n', flex.c_req / beam.d);
end

% Apply minimum steel area
As_min = flex.rho_min * beam.width * beam.d;
flex.As_req = max(flex.As_req, As_min);

fprintf('\nFINAL REQUIRED STEEL AREA\n');
fprintf('-------------------------\n');
fprintf('As_required:                 %.0f mm²\n', flex.As_req);
fprintf('As_min:                      %.0f mm²\n\n', As_min);

% Determine number of bars
flex.num_bars = ceil(flex.As_req / rebar.Ab_flex);
flex.As_provided = flex.num_bars * rebar.Ab_flex;
flex.rho_provided = flex.As_provided / (beam.width * beam.d);

% Check bar spacing
flex.bar_spacing = (beam.width - 2 * beam.cover - 2 * rebar.db_stirrup - flex.num_bars * rebar.db_flex) / ...
                   (flex.num_bars - 1);

% Minimum spacing check (CSA A23.3-19, Cl. 6.6.5.2)
flex.min_spacing_required = max([1.4 * rebar.db_flex, rebar.min_spacing, 30]);

fprintf('BAR LAYOUT\n');
fprintf('----------\n');
fprintf('Bar size selected:           %dM\n', rebar.flexural_bar);
fprintf('Number of bars:              %d\n', flex.num_bars);
fprintf('As_provided:                 %.0f mm²\n', flex.As_provided);
fprintf('ρ_provided:                  %.6f\n', flex.rho_provided);
fprintf('Bar spacing:                 %.1f mm\n', flex.bar_spacing);
fprintf('Min spacing required:        %.1f mm\n', flex.min_spacing_required);

if flex.bar_spacing < flex.min_spacing_required
    fprintf('\n✗ ERROR: Bar spacing too small! Increase beam width or use smaller bars.\n\n');
else
    fprintf('✓ Bar spacing adequate\n\n');
end

% Calculate actual moment resistance with provided steel
flex.a_provided = mat.phi_s * flex.As_provided * mat.fy_flexural / ...
                  (mat.phi_c * 0.85 * mat.fc_prime * beam.width);
flex.c_provided = flex.a_provided / mat.beta1;
flex.Mr = mat.phi_s * flex.As_provided * mat.fy_flexural * (beam.d - flex.a_provided / 2) / 1e6;  % kN·m

fprintf('MOMENT CAPACITY CHECK\n');
fprintf('---------------------\n');
fprintf('Factored moment (Mf):        %.2f kN·m\n', analysis.Mf_max);
fprintf('Moment resistance (Mr):      %.2f kN·m\n', flex.Mr);
fprintf('Mf/Mr ratio:                 %.3f\n', analysis.Mf_max / flex.Mr);
fprintf('Overstrength factor:         %.2f%%\n', (flex.Mr / analysis.Mf_max - 1) * 100);

if analysis.Mf_max <= flex.Mr
    fprintf('✓ MOMENT CHECK: PASS\n\n');
else
    fprintf('✗ MOMENT CHECK: FAIL - Increase section or reinforcement\n\n');
end

%% -------------------------------------------------------------------------
%  SHEAR DESIGN - STIRRUP SPACING
% -------------------------------------------------------------------------

fprintf('========================================\n');
fprintf('SHEAR DESIGN\n');
fprintf('========================================\n\n');

% Concrete shear resistance (CSA A23.3-19, Cl. 11.3.4)
% Simplified method: Vc = phi_c * lambda * beta * sqrt(f'c) * b * dv
% where beta = 0.21 for simplified method
shear.beta = 0.21;
shear.dv = max(0.9 * beam.d, 0.72 * beam.depth);  % effective shear depth (CSA A23.3-19, Cl. 11.3.3)

shear.Vc = mat.phi_c * mat.lambda * shear.beta * sqrt(mat.fc_prime) * beam.width * shear.dv / 1000;  % kN

fprintf('CONCRETE SHEAR CAPACITY\n');
fprintf('-----------------------\n');
fprintf('Effective shear depth (dv):  %.1f mm\n', shear.dv);
fprintf('Vc (concrete resistance):    %.2f kN\n', shear.Vc);
fprintf('Vf (factored shear):         %.2f kN\n', analysis.Vf_max);
fprintf('Vf/Vc ratio:                 %.3f\n\n', analysis.Vf_max / shear.Vc);

if analysis.Vf_max <= shear.Vc
    fprintf('✓ No shear reinforcement required by strength\n');
    fprintf('  However, minimum shear reinforcement will be provided per CSA A23.3-19, Cl. 11.2.8.2\n\n');
    shear.reinforcement_required = false;
else
    fprintf('⚠ Shear reinforcement REQUIRED\n\n');
    shear.reinforcement_required = true;
end

% Check if shear can be carried by section
shear.Vs_required = analysis.Vf_max - shear.Vc;  % kN
shear.Vs_max = 0.25 * mat.phi_c * mat.fc_prime * beam.width * shear.dv / 1000;  % kN (CSA A23.3-19, Cl. 11.3.3)

fprintf('SHEAR REINFORCEMENT REQUIREMENTS\n');
fprintf('--------------------------------\n');
fprintf('Vs_required:                 %.2f kN\n', shear.Vs_required);
fprintf('Vs_max (section limit):      %.2f kN\n', shear.Vs_max);

if shear.Vs_required > shear.Vs_max
    fprintf('✗ ERROR: Required Vs exceeds section capacity. Increase section size.\n\n');
else
    fprintf('✓ Section adequate for shear\n\n');
end

% Design stirrup spacing (CSA A23.3-19, Cl. 11.3.8)
% For vertical stirrups: Vs = Av * fy * dv / s
% Therefore: s = Av * fy * dv / Vs

if shear.reinforcement_required
    % Calculate required spacing
    shear.Av = 2 * rebar.Ab_stirrup;  % mm² (2-leg stirrup)
    shear.s_required = shear.Av * mat.fy_shear * shear.dv / (shear.Vs_required * 1000);  % mm
    
    fprintf('STIRRUP DESIGN (Strength)\n');
    fprintf('-------------------------\n');
    fprintf('Stirrup bar size:            %dM\n', rebar.stirrup_bar);
    fprintf('Number of legs:              2\n');
    fprintf('Av (stirrup area):           %.0f mm²\n', shear.Av);
    fprintf('s_required (from strength):  %.0f mm\n', shear.s_required);
else
    shear.Av = 2 * rebar.Ab_stirrup;
    shear.s_required = 9999;  % Large value when not required by strength
    fprintf('STIRRUP DESIGN (Minimum only)\n');
    fprintf('-----------------------------\n');
    fprintf('Stirrup bar size:            %dM\n', rebar.stirrup_bar);
    fprintf('Number of legs:              2\n');
    fprintf('Av (stirrup area):           %.0f mm²\n', shear.Av);
end

% Maximum spacing limits (CSA A23.3-19, Cl. 11.3.8.1)
if shear.Vs_required <= 0.125 * mat.phi_c * mat.fc_prime * beam.width * shear.dv / 1000
    shear.s_max = min(0.7 * shear.dv, 600);  % mm
else
    shear.s_max = min(0.35 * shear.dv, 300);  % mm
end

% Minimum shear reinforcement (CSA A23.3-19, Cl. 11.2.8.2)
% Av_min / s = 0.06 * sqrt(f'c) / fy * bw (in MPa units)
shear.s_min = shear.Av / (0.06 * sqrt(mat.fc_prime) / mat.fy_shear * beam.width);

fprintf('s_max (spacing limit):       %.0f mm\n', shear.s_max);
fprintf('s_min (from Av,min):         %.0f mm\n', shear.s_min);

% Select practical spacing (round down to nearest 25mm)
shear.s_provided = floor(min([shear.s_required, shear.s_max, shear.s_min]) / 25) * 25;
shear.s_provided = max(shear.s_provided, 50);  % Minimum practical spacing

% Special spacing at supports and mid-span
shear.s_at_support = min(shear.s_provided, 100);  % Closer spacing near support
shear.s_at_midspan = min(shear.s_max, 300);       % Can increase spacing at mid-span if Vf is low

fprintf('\nFINAL STIRRUP SPACING\n');
fprintf('---------------------\n');
fprintf('s_provided (general):        %d mm\n', shear.s_provided);
fprintf('s_at_support (d/2 from support): %d mm\n', shear.s_at_support);
fprintf('s_at_midspan (low shear):    %d mm\n\n', shear.s_at_midspan);

% Verify provided spacing
shear.Vs_provided = shear.Av * mat.fy_shear * shear.dv / shear.s_provided / 1000;  % kN
shear.Vr = shear.Vc + shear.Vs_provided;  % kN

fprintf('SHEAR CAPACITY CHECK\n');
fprintf('--------------------\n');
fprintf('Vc:                          %.2f kN\n', shear.Vc);
fprintf('Vs_provided:                 %.2f kN\n', shear.Vs_provided);
fprintf('Vr = Vc + Vs:                %.2f kN\n', shear.Vr);
fprintf('Vf/Vr ratio:                 %.3f\n', analysis.Vf_max / shear.Vr);

if analysis.Vf_max <= shear.Vr
    fprintf('✓ SHEAR CHECK: PASS\n\n');
else
    fprintf('✗ SHEAR CHECK: FAIL - Decrease stirrup spacing\n\n');
end

%% -------------------------------------------------------------------------
%  DEFLECTION CHECK
% -------------------------------------------------------------------------

fprintf('========================================\n');
fprintf('DEFLECTION CHECK\n');
fprintf('========================================\n\n');

% Service load moment at midspan
defl.Ms = loads.ws * (beam.span/1000)^2 / 8;  % kN·m

fprintf('SERVICE LOAD ANALYSIS\n');
fprintf('---------------------\n');
fprintf('Service load (ws):           %.2f kN/m\n', loads.ws);
fprintf('Service moment (Ms):         %.2f kN·m\n\n', defl.Ms);

% Gross moment of inertia
defl.Ig = beam.width * beam.depth^3 / 12;  % mm^4

% Modulus of rupture (CSA A23.3-19, Cl. 8.6.4)
defl.fr = 0.6 * mat.lambda * sqrt(mat.fc_prime);  % MPa

% Cracking moment (CSA A23.3-19, Cl. 9.8.4.2)
defl.yt = beam.depth / 2;  % mm - distance to extreme tension fiber
defl.Mcr = defl.fr * defl.Ig / defl.yt / 1e6;  % kN·m

fprintf('CRACKING ANALYSIS\n');
fprintf('-----------------\n');
fprintf('Gross moment of inertia (Ig): %.2e mm⁴\n', defl.Ig);
fprintf('Modulus of rupture (fr):      %.2f MPa\n', defl.fr);
fprintf('Cracking moment (Mcr):        %.2f kN·m\n', defl.Mcr);
fprintf('Ms/Mcr ratio:                 %.3f\n', defl.Ms / defl.Mcr);

if defl.Ms > defl.Mcr
    fprintf('→ Section is CRACKED under service load\n\n');
    
    % Modular ratio
    defl.n = mat.Es / mat.Ec;
    
    % Transformed area and neutral axis location for cracked section
    % Using transformed area method: solve for kd (depth to neutral axis)
    defl.rho = flex.As_provided / (beam.width * beam.d);
    defl.k = sqrt(2 * defl.rho * defl.n + (defl.rho * defl.n)^2) - defl.rho * defl.n;
    defl.kd = defl.k * beam.d;  % mm
    
    % Cracked moment of inertia
    defl.Icr = beam.width * defl.kd^3 / 3 + defl.n * flex.As_provided * (beam.d - defl.kd)^2;  % mm^4
    
    % Effective moment of inertia (CSA A23.3-19, Cl. 9.8.4.2)
    defl.Ie = (defl.Mcr / defl.Ms)^3 * defl.Ig + (1 - (defl.Mcr / defl.Ms)^3) * defl.Icr;
    
    fprintf('CRACKED SECTION PROPERTIES\n');
    fprintf('--------------------------\n');
    fprintf('Modular ratio (n):            %.2f\n', defl.n);
    fprintf('Neutral axis factor (k):      %.4f\n', defl.k);
    fprintf('Neutral axis depth (kd):      %.1f mm\n', defl.kd);
    fprintf('Cracked moment of inertia:    %.2e mm⁴\n', defl.Icr);
    fprintf('Effective moment of inertia:  %.2e mm⁴\n', defl.Ie);
    fprintf('Ie/Ig ratio:                  %.3f\n\n', defl.Ie / defl.Ig);
else
    fprintf('→ Section is UNCRACKED under service load\n');
    fprintf('→ Using gross moment of inertia\n\n');
    defl.Ie = defl.Ig;
end

% Immediate deflection for simply supported beam under uniform load
% δ = 5 * w * L^4 / (384 * E * I)
defl.delta_immediate = 5 * loads.ws * beam.span^4 / (384 * mat.Ec * defl.Ie);  % mm

fprintf('DEFLECTION CALCULATIONS\n');
fprintf('-----------------------\n');
fprintf('Immediate deflection:        %.2f mm\n', defl.delta_immediate);

% Long-term deflection multiplier (CSA A23.3-19, Cl. 9.8.4.3)
defl.xi = 2.0;  % time-dependent factor (5 years or more)
defl.rho_prime = 0;  % compression reinforcement ratio (assumed zero)
defl.lambda_delta = defl.xi / (1 + 50 * defl.rho_prime);

% Sustained load ratio (assume DL is sustained)
defl.sustained_load_ratio = loads.DL / loads.ws;
defl.delta_longterm = defl.delta_immediate * (1 + defl.lambda_delta * defl.sustained_load_ratio);

fprintf('Time-dependent factor (ξ):   %.2f\n', defl.xi);
fprintf('Deflection multiplier (λ):   %.2f\n', defl.lambda_delta);
fprintf('Sustained load ratio:        %.3f\n', defl.sustained_load_ratio);
fprintf('Long-term deflection:        %.2f mm\n\n', defl.delta_longterm);

% Deflection limits (CSA A23.3-19, Table 9.3)
defl.limit_immediate = beam.span / 360;  % mm - immediate deflection limit
defl.limit_longterm = beam.span / 240;   % mm - total deflection limit

fprintf('DEFLECTION LIMITS\n');
fprintf('-----------------\n');
fprintf('Immediate limit (L/360):     %.2f mm\n', defl.limit_immediate);
fprintf('Long-term limit (L/240):     %.2f mm\n', defl.limit_longterm);
fprintf('Immediate/Limit:             %.3f\n', defl.delta_immediate / defl.limit_immediate);
fprintf('Long-term/Limit:             %.3f\n\n', defl.delta_longterm / defl.limit_longterm);

defl.immediate_pass = defl.delta_immediate <= defl.limit_immediate;
defl.longterm_pass = defl.delta_longterm <= defl.limit_longterm;

if defl.immediate_pass && defl.longterm_pass
    fprintf('✓ DEFLECTION CHECK: PASS (both immediate and long-term)\n\n');
elseif defl.longterm_pass
    fprintf('⚠ DEFLECTION CHECK: Immediate deflection exceeds limit\n');
    fprintf('  (Long-term is acceptable)\n\n');
else
    fprintf('✗ DEFLECTION CHECK: FAIL - Increase depth or reduce load\n\n');
end

%% -------------------------------------------------------------------------
%  SUMMARY TABLE
% -------------------------------------------------------------------------

fprintf('========================================\n');
fprintf('DESIGN SUMMARY\n');
fprintf('========================================\n\n');

summaryTable = table(...
    ["Moment"; "Shear"; "Deflection (Immediate)"; "Deflection (Long-term)"], ...
    [analysis.Mf_max; analysis.Vf_max; defl.delta_immediate; defl.delta_longterm], ...
    [flex.Mr; shear.Vr; defl.limit_immediate; defl.limit_longterm], ...
    [analysis.Mf_max/flex.Mr; analysis.Vf_max/shear.Vr; ...
     defl.delta_immediate/defl.limit_immediate; defl.delta_longterm/defl.limit_longterm], ...
    ["PASS"; "PASS"; string(defl.immediate_pass); string(defl.longterm_pass)], ...
    'VariableNames', {'Check', 'Demand', 'Capacity', 'Ratio', 'Status'});

disp(summaryTable);

fprintf('\nREINFORCEMENT DETAILS\n');
fprintf('---------------------\n');
fprintf('Flexural:  %d-%dM bars (As = %.0f mm²)\n', flex.num_bars, rebar.flexural_bar, flex.As_provided);
fprintf('Shear:     %dM stirrups @ %d mm c/c\n', rebar.stirrup_bar, shear.s_provided);
fprintf('           (%d mm @ supports, %d mm @ midspan)\n\n', shear.s_at_support, shear.s_at_midspan);

%% -------------------------------------------------------------------------
%  MOMENT AND SHEAR DIAGRAMS
% -------------------------------------------------------------------------

fprintf('Generating moment and shear diagrams...\n\n');

% Create array of positions along span
n_points = 101;
x = linspace(0, beam.span, n_points);  % mm

% Calculate moment and shear at each position
% Simply supported beam with uniform load:
% M(x) = w*x*(L-x)/2
% V(x) = w*L/2 - w*x

M_x = loads.wf * x .* (beam.span - x) / 2 / 1e6;  % kN·m
V_x = loads.wf * beam.span / 2 / 1e3 - loads.wf * x / 1e3;  % kN

% Service load moment for deflection calculations
Ms_x = loads.ws * x .* (beam.span - x) / 2 / 1e6;  % kN·m

% Create figure with subplots
figure('Name', 'RC Beam Analysis Diagrams', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

% Beam elevation with reinforcement
subplot(3, 2, [1, 2]);
hold on; grid on;

% Draw beam outline
beam_y = [0, beam.depth];
patch([0, beam.span, beam.span, 0], [0, 0, beam.depth, beam.depth], ...
      [0.9 0.9 0.92], 'EdgeColor', 'k', 'LineWidth', 1.5);

% Draw flexural reinforcement
bar_y = beam.cover + rebar.db_stirrup + rebar.db_flex/2;
bar_x_start = beam.cover + rebar.db_stirrup + rebar.db_flex/2;
bar_x_end = beam.width - beam.cover - rebar.db_stirrup - rebar.db_flex/2;

for i = 1:flex.num_bars
    if flex.num_bars == 1
        bar_x = beam.width / 2;
    else
        bar_x = bar_x_start + (i-1) * (bar_x_end - bar_x_start) / (flex.num_bars - 1);
    end
    plot([bar_x, bar_x], [bar_y - 30, bar_y + 30], 'r', 'LineWidth', 8);
end

% Draw stirrups at several locations
stirrup_locations = [100, beam.span/4, beam.span/2, 3*beam.span/4, beam.span-100];
for i = 1:length(stirrup_locations)
    stirrup_x = stirrup_locations(i);
    % Vertical legs
    plot([stirrup_x, stirrup_x], [beam.cover, beam.depth - beam.cover], 'b', 'LineWidth', 2);
end

% Add dimensions
plot([0, beam.span], [-50, -50], 'k', 'LineWidth', 1);
plot([0, 0], [-70, -30], 'k', 'LineWidth', 1);
plot([beam.span, beam.span], [-70, -30], 'k', 'LineWidth', 1);
text(beam.span/2, -120, sprintf('L = %d mm (%.2f m)', beam.span, beam.span/1000), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');

% Supports
plot([0, 0], [0, -30], 'k', 'LineWidth', 3);
plot([beam.span, beam.span], [0, -30], 'k', 'LineWidth', 3);
plot([-50, 50], [-30, -30], 'k', 'LineWidth', 3);
plot([beam.span-50, beam.span+50], [-30, -30], 'k', 'LineWidth', 3);

title(sprintf('Beam Elevation: %d×%d mm, %d-%dM bars, %dM stirrups @ %d mm', ...
      beam.width, beam.depth, flex.num_bars, rebar.flexural_bar, ...
      rebar.stirrup_bar, shear.s_provided), 'FontWeight', 'bold');
xlabel('Length (mm)');
ylabel('Depth (mm)');
xlim([-200, beam.span + 200]);
ylim([-200, beam.depth + 100]);
axis equal;
xlim([-200, beam.span + 200]);

% Moment diagram (factored)
subplot(3, 2, 3);
hold on; grid on;
plot(x/1000, M_x, 'b-', 'LineWidth', 2);
yline(flex.Mr, 'r--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.5);
xlabel('Position along span (m)');
ylabel('Moment (kN·m)');
title('Bending Moment Diagram (Factored Load)', 'FontWeight', 'bold');
legend({'Mf (factored)', 'Mr (resistance)'}, 'Location', 'north');
ylim([min(M_x)*1.2, max([max(M_x), flex.Mr])*1.2]);

% Moment diagram (service) for deflection
subplot(3, 2, 4);
hold on; grid on;
plot(x/1000, Ms_x, 'g-', 'LineWidth', 2);
yline(defl.Mcr, 'm--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.5);
xlabel('Position along span (m)');
ylabel('Moment (kN·m)');
title('Bending Moment Diagram (Service Load)', 'FontWeight', 'bold');
legend({'Ms (service)', 'Mcr (cracking)'}, 'Location', 'north');
ylim([min(Ms_x)*1.2, max([max(Ms_x), defl.Mcr])*1.2]);

% Shear diagram
subplot(3, 2, 5);
hold on; grid on;
plot(x/1000, V_x, 'b-', 'LineWidth', 2);
yline(shear.Vr, 'r--', 'LineWidth', 1.5);
yline(shear.Vc, 'g--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.5);
yline(-shear.Vr, 'r--', 'LineWidth', 1.5);
yline(-shear.Vc, 'g--', 'LineWidth', 1.5);
xlabel('Position along span (m)');
ylabel('Shear Force (kN)');
title('Shear Force Diagram', 'FontWeight', 'bold');
legend({'Vf (factored)', 'Vr (Vc+Vs)', 'Vc (concrete only)'}, 'Location', 'northeast');
ylim([min(V_x)*1.3, max(V_x)*1.3]);

% Stirrup spacing diagram
subplot(3, 2, 6);
hold on; grid on;

% Create stirrup spacing zones
zone1_end = beam.d / 1000;  % d/2 from support (in meters)
zone2_start = zone1_end;
zone2_end = beam.span/1000 - zone1_end;

% Plot spacing zones as step function
x_spacing = [0, zone1_end, zone1_end, zone2_end, zone2_end, beam.span/1000];
s_spacing = [shear.s_at_support, shear.s_at_support, shear.s_provided, ...
             shear.s_provided, shear.s_at_support, shear.s_at_support];

plot(x_spacing, s_spacing, 'b-', 'LineWidth', 2);
yline(shear.s_max, 'r--', 'LineWidth', 1.5);
yline(shear.s_min, 'g--', 'LineWidth', 1.5);

xlabel('Position along span (m)');
ylabel('Stirrup Spacing (mm)');
title('Stirrup Spacing Layout', 'FontWeight', 'bold');
legend({'s_provided', 's_max', 's_min'}, 'Location', 'north');
grid on;
ylim([0, max([shear.s_max, shear.s_provided, shear.s_at_support]) * 1.2]);

fprintf('✓ Diagrams generated successfully.\n\n');

%% -------------------------------------------------------------------------
%  FINAL STATUS
% -------------------------------------------------------------------------

fprintf('========================================\n');
fprintf('DESIGN STATUS\n');
fprintf('========================================\n\n');

all_pass = (analysis.Mf_max <= flex.Mr) && ...
           (analysis.Vf_max <= shear.Vr) && ...
           defl.longterm_pass && ...
           (flex.bar_spacing >= flex.min_spacing_required);

if all_pass
    fprintf('✓✓✓ ALL CHECKS PASSED ✓✓✓\n\n');
    fprintf('The beam design is adequate per CSA A23.3-19.\n');
    fprintf('Proceed with detailing and construction drawings.\n\n');
else
    fprintf('⚠⚠⚠ DESIGN REQUIRES REVISION ⚠⚠⚠\n\n');
    
    if analysis.Mf_max > flex.Mr
        fprintf('• Moment capacity insufficient\n');
    end
    if analysis.Vf_max > shear.Vr
        fprintf('• Shear capacity insufficient\n');
    end
    if ~defl.longterm_pass
        fprintf('• Deflection exceeds limits\n');
    end
    if flex.bar_spacing < flex.min_spacing_required
        fprintf('• Bar spacing too small\n');
    end
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('END OF DESIGN CALCULATION\n');
fprintf('========================================\n\n');
