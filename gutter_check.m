%% WATER_RETAINING_WALL_CSA
% CSA-based stability analysis for a water-retaining wall section.
% Checks overturning, sliding, and bearing capacity for the supplied geometry
% and hydrostatic loading. All forces are reported per metre length of wall.
%
% Coordinate system: toe edge at x = 0, positive x runs toward the heel.
% Several geometric values (e.g., stem and base thickness) are assumed and
% should be updated if project-specific data are available.

clear; close all; clc;

mm = 1e-3;

%% -------------------------------------------------------------------------
%  INPUT PARAMETERS
% -------------------------------------------------------------------------

% Geometry (metres unless noted)
geom.stemHeight = 3050 * mm;
geom.heelWidth  = 1300 * mm;
geom.toeWidth   = 610  * mm;
geom.stemThk    = 0.35;         % average constant thickness
geom.baseThk    = 0.4;         % footing thickness
geom.keyWidth   = 350  * mm;
geom.keyDepth   = 370  * mm;

geom.gutter.offsetBelowTop = 400 * mm;
geom.gutter.horizWidth     = 1000 * mm;
geom.gutter.horizThk       = 150 * mm;
geom.gutter.vertHeight     = 290 * mm;
geom.gutter.vertThk        = 150 * mm;

geom.roof.thk           = 250 * mm;
geom.roof.tributarySpan = 3.0; % m

% Optional heel floor slab weight (set thickness > 0 to include)
geom.floorSlabThk   = 0.0;      % m
geom.floorSlabWidth = geom.heelWidth;

% Hydrostatic loading
loads.waterDepth = 0.0 * mm;   % measured from heel slab surface

% Material properties
mat.gammaConcrete = 25.0;       % kN/m^3
mat.gammaWater    = 9.81;       % kN/m^3
mat.gammaSoil     = 19.0;       % kN/m^3
mat.phi           = deg2rad(31);
mat.delta         = deg2rad(15);

% CSA target factors of safety(temporary works)
criteria.overturning = 1.20;
criteria.sliding     = 1.20;
criteria.bearing     = 2.0;

ls.phi_bearing = 0.5;
ls.gamma_dead  = 1.25;
ls.gamma_hydro = 1.00;

%% -------------------------------------------------------------------------
%  DERIVED GEOMETRY
% -------------------------------------------------------------------------

geom.baseWidth        = geom.toeWidth + geom.heelWidth + geom.stemThk;
geom.stemCL           = geom.toeWidth + geom.stemThk/2;
geom.baseTopElevation = geom.baseThk;
geom.stemLeft         = geom.stemCL - geom.stemThk/2;
geom.stemRight        = geom.stemCL + geom.stemThk/2;

%% -------------------------------------------------------------------------
%  GRAVITY LOAD COMPONENTS
% -------------------------------------------------------------------------

components = struct('Name', {}, 'Weight', {}, 'LeverArm', {}, 'Moment', {}, 'Notes', {});

components = addComponent(components, "Base slab", ...
    mat.gammaConcrete * geom.baseWidth * geom.baseThk, geom.baseWidth / 2, ...
    "Footing self-weight");

components = addComponent(components, "Shear key", ...
    mat.gammaConcrete * geom.keyWidth * geom.keyDepth, geom.stemCL, ...
    "Key volume below footing");

components = addComponent(components, "Stem", ...
    mat.gammaConcrete * geom.stemThk * geom.stemHeight, geom.stemCL, ...
    "Rectangular stem");

components = addComponent(components, "Gutter (horizontal)", ...
    mat.gammaConcrete * geom.gutter.horizWidth * geom.gutter.horizThk, ...
    geom.stemCL - geom.gutter.horizWidth / 2, ...
    "Toe-side gutter slab");

components = addComponent(components, "Gutter (vertical)", ...
    mat.gammaConcrete * geom.gutter.vertHeight * geom.gutter.vertThk, ...
    geom.stemCL - geom.gutter.horizWidth + geom.gutter.vertThk / 2, ...
    "Toe-side curb");

components = addComponent(components, "Roof slab tributary load", ...
    mat.gammaConcrete * geom.roof.thk * geom.roof.tributarySpan, geom.stemCL, ...
    "Roof weight carried by stem");

if geom.floorSlabThk > 0
    components = addComponent(components, "Heel floor slab (optional)", ...
        mat.gammaConcrete * geom.floorSlabThk * geom.floorSlabWidth, ...
        geom.stemCL + geom.heelWidth / 2, ...
        "Adjust thickness if present");
end

weights     = [components.Weight]';
leverArms   = [components.LeverArm]';
moments     = [components.Moment]';

W_total     = sum(weights);
M_resisting = sum(moments);

%% -------------------------------------------------------------------------
%  HYDROSTATIC ACTIONS
% -------------------------------------------------------------------------

hydro.thrust = 0.5 * mat.gammaWater * loads.waterDepth^2;            % kN/m
hydro.arm    = geom.baseTopElevation + loads.waterDepth / 3;          % m
hydro.moment = hydro.thrust * hydro.arm;                              % kN·m/m

V_net = W_total;                                                       % kN/m (downwards positive)

%% -------------------------------------------------------------------------
%  STABILITY CHECKS (CSA)
% -------------------------------------------------------------------------

FS_overturning = M_resisting / hydro.moment;

friction_resistance = V_net * tan(mat.delta);
Kp = tan(pi/4 + mat.phi/2)^2;
passive_resistance = 0.5 * Kp * mat.gammaSoil * geom.keyDepth^2 * geom.keyWidth;

FS_sliding = (friction_resistance + passive_resistance) / hydro.thrust;

M_net_about_toe = M_resisting - hydro.moment;

if V_net <= 0
    warning('Net vertical load is non-positive. Sliding and bearing checks are invalid.');
    x_resultant  = NaN;
    eccentricity = NaN;
    q_avg        = NaN;
    q_toe        = NaN;
    q_heel       = NaN;
else
    x_resultant  = M_net_about_toe / V_net;                  % m from toe
    eccentricity = geom.baseWidth / 2 - x_resultant;         % +ve toward toe
    q_avg        = V_net / geom.baseWidth;                   % kPa
    q_toe        = q_avg * (1 + 6 * eccentricity / geom.baseWidth);
    q_heel       = q_avg * (1 - 6 * eccentricity / geom.baseWidth);
end

W_uls              = ls.gamma_dead * W_total;
M_resisting_uls    = ls.gamma_dead * M_resisting;
hydro_uls.thrust   = ls.gamma_hydro * hydro.thrust;
hydro_uls.moment   = ls.gamma_hydro * hydro.moment;
V_net_uls          = W_uls;
M_net_about_toe_uls = M_resisting_uls - hydro_uls.moment;

if V_net_uls <= 0
    x_resultant_uls  = NaN;
    eccentricity_uls = NaN;
    q_avg_uls        = NaN;
    q_toe_uls        = NaN;
    q_heel_uls       = NaN;
else
    x_resultant_uls  = M_net_about_toe_uls / V_net_uls;
    eccentricity_uls = geom.baseWidth / 2 - x_resultant_uls;
    q_avg_uls        = V_net_uls / geom.baseWidth;
    q_toe_uls        = q_avg_uls * (1 + 6 * eccentricity_uls / geom.baseWidth);
    q_heel_uls       = q_avg_uls * (1 - 6 * eccentricity_uls / geom.baseWidth);
end

q_pair = [q_toe, q_heel];
valid_q = q_pair(~isnan(q_pair));
if isempty(valid_q)
    q_max = NaN;
    q_min = NaN;
else
    q_max = max(valid_q);
    q_min = min(valid_q);
end

Nq     = exp(pi * tan(mat.phi)) * tan(pi/4 + mat.phi/2)^2;
Ngamma = 2 * (Nq + 1) * tan(mat.phi);

% ULS effective width method (CFEM): no-tension, use B_eff = B - 2|e_uls|
if isnan(eccentricity_uls)
    B_eff_uls = NaN;
    q_dem_uls = NaN;
    phi_q_ult = NaN;
else
    B_eff_uls = max(geom.baseWidth - 2 * abs(eccentricity_uls), 1e-6);
    q_ult_eff = 0.5 * mat.gammaSoil * B_eff_uls * Ngamma;   % kPa
    phi_q_ult = ls.phi_bearing * q_ult_eff;                 % kPa
    q_dem_uls = V_net_uls / B_eff_uls;                      % kPa
end

% Also report ULS linear pressures (may show tension if e_uls > B/6)
if isnan(q_toe_uls) || isnan(q_heel_uls)
    q_max_uls = NaN;
    q_min_uls = NaN;
else
    q_max_uls = max([q_toe_uls, q_heel_uls]);
    q_min_uls = min([q_toe_uls, q_heel_uls]);
end

if isnan(q_dem_uls) || q_dem_uls <= 0
    FS_bearing = NaN;
else
    FS_bearing = phi_q_ult / q_dem_uls;
end

if isnan(eccentricity)
    within_middle_third = false;
else
    within_middle_third = abs(eccentricity) <= geom.baseWidth / 6;
end

%% -------------------------------------------------------------------------
%  TABULATED RESULTS
% -------------------------------------------------------------------------

componentTable = table(string({components.Name})', weights, leverArms, moments, string({components.Notes})', ...
    'VariableNames', {'Component', 'Weight_kN_per_m', 'LeverArm_m', 'Moment_kNm_per_m', 'Notes'});

FS_values   = [FS_overturning; FS_sliding; FS_bearing];
FS_required = [criteria.overturning; criteria.sliding; 1.0];
statusStr   = repmat("PASS", size(FS_values));
statusStr(FS_values < FS_required) = "CHECK";

resultsTable = table(["Overturning"; "Sliding"; "Bearing"], FS_values, FS_required, statusStr, ...
    'VariableNames', {'Check', 'FS_Computed', 'FS_Required', 'Status'});

bearingTable = table(q_toe, q_heel, q_avg, q_min, q_max, q_toe_uls, q_heel_uls, q_max_uls, B_eff_uls, q_dem_uls, phi_q_ult, within_middle_third, ...
    'VariableNames', {'q_toe_SLS_kPa', 'q_heel_SLS_kPa', 'q_avg_SLS_kPa', 'q_min_SLS_kPa', 'q_max_SLS_kPa', ...
                      'q_toe_ULS_kPa', 'q_heel_ULS_kPa', 'q_max_ULS_kPa', 'B_eff_ULS_m', 'q_dem_ULS_kPa', 'phi_q_ult_kPa', 'WithinMiddleThird'});

%% -------------------------------------------------------------------------
%  CONSOLE OUTPUT
% -------------------------------------------------------------------------

fprintf('\nCSA Water-Retaining Wall Stability Check (per metre length)\n');
fprintf('-----------------------------------------------------------\n');
fprintf('Toe width                : %.3f m\n', geom.toeWidth);
fprintf('Heel width               : %.3f m\n', geom.heelWidth);
fprintf('Overall base width       : %.3f m\n', geom.baseWidth);
fprintf('Stem height              : %.3f m\n', geom.stemHeight);
fprintf('Water depth at heel      : %.3f m\n', loads.waterDepth);
fprintf('Assumed base thickness   : %.0f mm\n', geom.baseThk * 1e3);
fprintf('Assumed stem thickness   : %.0f mm\n\n', geom.stemThk * 1e3);

fprintf('Total dead load          : %.2f kN/m\n', W_total);
fprintf('Hydrostatic thrust       : %.2f kN/m\n', hydro.thrust);
fprintf('Net vertical load        : %.2f kN/m\n\n', V_net);

disp(componentTable);

fprintf('Stability Factors of Safety\n');
fprintf('---------------------------\n');
disp(resultsTable);

fprintf('Bearing summary (toe positive)\n');
fprintf('-------------------------------\n');
disp(bearingTable);

fprintf('Resultant from toe        : %.3f m\n', x_resultant);
fprintf('Eccentricity (toe +ve)    : %.3f m (limit %.3f m)\n', eccentricity, geom.baseWidth/6);

if FS_overturning < criteria.overturning
    fprintf('WARNING: Overturning FS below CSA minimum (%.2f < %.2f).\n', FS_overturning, criteria.overturning);
end
if FS_sliding < criteria.sliding
    fprintf('WARNING: Sliding FS below CSA minimum (%.2f < %.2f).\n', FS_sliding, criteria.sliding);
end
if FS_bearing < 1.0
    fprintf('WARNING: Bearing ULS check not satisfied (FS = %.2f < 1.00). Factored demand exceeds factored resistance.\n', FS_bearing);
end
if ~within_middle_third
    fprintf('WARNING: Resultant lies outside the middle third; heel tension likely.\n');
end
if ~isnan(q_min) && q_min < 0
    fprintf('WARNING: Heel contact pressure negative (q_{min} = %.2f kPa). Revise geometry/loading.\n', q_min);
end

%% -------------------------------------------------------------------------
%  VISUALISATIONS
% -------------------------------------------------------------------------

figure('Name', 'Wall Geometry & Hydrostatics', 'Color', 'w');
hold on; axis equal;

% Base slab
basePoly = [0, 0;
            geom.baseWidth, 0;
            geom.baseWidth, geom.baseTopElevation;
            0, geom.baseTopElevation];
patch(basePoly(:,1), basePoly(:,2), [0.85 0.85 0.88], 'EdgeColor', 'k');

% Shear key
keyLeft  = geom.stemCL - geom.keyWidth/2;
keyRight = geom.stemCL + geom.keyWidth/2;
keyPoly = [keyLeft, 0;
           keyRight, 0;
           keyRight, -geom.keyDepth;
           keyLeft, -geom.keyDepth];
patch(keyPoly(:,1), keyPoly(:,2), [0.75 0.75 0.80], 'EdgeColor', 'k');

% Stem
stemPoly = [geom.stemLeft,  geom.baseTopElevation;
            geom.stemRight, geom.baseTopElevation;
            geom.stemRight, geom.baseTopElevation + geom.stemHeight;
            geom.stemLeft,  geom.baseTopElevation + geom.stemHeight];
patch(stemPoly(:,1), stemPoly(:,2), [0.90 0.90 0.94], 'EdgeColor', 'k');

% Gutter horizontal slab (toe side)
gutterTop    = geom.baseTopElevation + geom.stemHeight - geom.gutter.offsetBelowTop;
gutterBottom = gutterTop - geom.gutter.horizThk;
gutterLeft   = geom.stemCL - geom.gutter.horizWidth;
gutterRight  = geom.stemCL;
gutterHPoly  = [gutterLeft, gutterBottom;
                gutterRight, gutterBottom;
                gutterRight, gutterTop;
                gutterLeft, gutterTop];
patch(gutterHPoly(:,1), gutterHPoly(:,2), [0.80 0.88 0.90], 'EdgeColor', 'k');

% Gutter vertical curb
curbBottom = gutterTop;
curbTop    = gutterTop + geom.gutter.vertHeight;
curbLeft   = gutterLeft;
curbRight  = gutterLeft + geom.gutter.vertThk;
curbPoly   = [curbLeft, curbBottom;
              curbRight, curbBottom;
              curbRight, curbTop;
              curbLeft, curbTop];
patch(curbPoly(:,1), curbPoly(:,2), [0.78 0.86 0.88], 'EdgeColor', 'k');

% Water retained on heel side
waterTop = geom.baseTopElevation + loads.waterDepth;
waterPoly = [geom.stemRight,           geom.baseTopElevation;
             geom.stemRight,           waterTop;
             geom.stemRight + 0.35,    waterTop;
             geom.stemRight + 0.35,    geom.baseTopElevation];
patch(waterPoly(:,1), waterPoly(:,2), [0.70 0.80 0.95], 'FaceAlpha', 0.6, 'EdgeColor', [0.3 0.4 0.7]);

% Hydrostatic pressure diagram
pressureScale = 0.012; % m per kPa for plotting
pressureAtBase = mat.gammaWater * loads.waterDepth;
pressurePoly = [geom.stemRight,                               geom.baseTopElevation;
                geom.stemRight + pressureAtBase * pressureScale, geom.baseTopElevation;
                geom.stemRight,                               waterTop];
patch(pressurePoly(:,1), pressurePoly(:,2), [0.40 0.55 0.90], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
plot([geom.stemRight, geom.stemRight], [geom.baseTopElevation, waterTop], 'b', 'LineWidth', 1.2);

% Resultant location marker
if ~isnan(x_resultant)
    plot([x_resultant, x_resultant], [-geom.keyDepth, geom.baseTopElevation], 'r--', 'LineWidth', 1.2);
    text(x_resultant, -geom.keyDepth - 0.05, sprintf('Resultant %.3f m from toe', x_resultant), ...
        'Color', 'r', 'HorizontalAlignment', 'center');
end

title('Wall Elevation & Hydrostatic Actions');
xlabel('Horizontal distance from toe (m)');
ylabel('Elevation (m)');
grid on;
xlim([gutterLeft - 0.4, geom.baseWidth + 0.5]);
ylim([-geom.keyDepth - 0.2, geom.baseTopElevation + geom.stemHeight + 0.4]);

% Stability summary plots -------------------------------------------------
figure('Name', 'CSA Stability Summary', 'Color', 'w');
tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Factors of safety bar chart
nexttile(tl, 1);
barData = [FS_values FS_required];
hb = bar(1:3, barData, 'grouped');
hb(1).FaceColor = [0.55 0.75 0.55];
hb(2).FaceColor = [0.85 0.60 0.60];
hold on; grid on;

for k = 1:numel(hb)
    xPts = hb(k).XEndPoints;
    yVals = hb(k).YData;
    for j = 1:numel(xPts)
        switch k
            case 1
                labelText = sprintf('%.2f', yVals(j));
                labelColor = [0 0 0];
                labelWeight = 'bold';
            case 2
                labelText = sprintf('Req %.2f', yVals(j));
                labelColor = [0.6 0.1 0.1];
                labelWeight = 'normal';
        end
        text(xPts(j), yVals(j) + 0.05, labelText, 'HorizontalAlignment', 'center', ...
            'Color', labelColor, 'FontWeight', labelWeight);
    end
end

set(gca, 'XTick', 1:3, 'XTickLabel', resultsTable.Check);
ylabel('Factor of Safety');
title('Factors of Safety vs CSA Targets');
legend({'Computed FS', 'CSA Target'}, 'Location', 'northwest');

% Bearing pressure profile
nexttile(tl, 2);
hold on; grid on;

if all(isnan(q_pair))
    q_plot = [0, 0];
else
    q_plot = [q_toe, q_heel];
end

plot([0, geom.baseWidth], q_plot, '-o', 'LineWidth', 2, 'MarkerFaceColor', [0.2 0.4 0.7]);
yline(0, 'k:');
yline(phi_q_ult, 'Color', [0.3 0.6 0.3], 'LineStyle', '--', 'LineWidth', 1.2);

if ~isnan(q_min) && q_min < 0
    text(geom.baseWidth * 0.8, min(q_plot) - 10, 'Heel tension (q_{min} < 0)', ...
        'Color', [0.6 0.1 0.1], 'FontWeight', 'bold');
end

xlabel('Distance from toe (m)');
ylabel('Bearing pressure (kPa)');
title('Base Contact Pressure');
xlim([0, geom.baseWidth + 0.4]);
ylim([min([0, q_plot, phi_q_ult], [], 'omitnan') - 20, max([q_plot, phi_q_ult], [], 'omitnan') + 20]);
legend({'Contact pressure (SLS)', 'Zero', 'Factored resistance (\phi q_{ult})'}, 'Location', 'best');

%% -------------------------------------------------------------------------
%  GUTTER STRUCTURAL ADEQUACY CHECK (CSA A23.3)
% -------------------------------------------------------------------------

fprintf('\n\nGutter Structural Adequacy Check per CSA A23.3\n');
fprintf('==============================================\n\n');

% Material properties for concrete and steel (CSA A23.3)
gutter.fc_prime = 30;              % MPa - specified compressive strength of concrete
gutter.fy = 400;                   % MPa - specified yield strength of reinforcement (Grade 400)
gutter.Es = 200000;                % MPa - modulus of elasticity of steel
gutter.Ec = 4500 * sqrt(gutter.fc_prime);  % MPa - modulus of elasticity of concrete (CSA A23.3-14, Cl. 8.6.2.2)

% Material resistance factors (CSA A23.3-14)
gutter.phi_c = 0.65;               % resistance factor for concrete
gutter.phi_s = 0.85;               % resistance factor for steel

% Geometric properties
gutter.span = geom.gutter.horizWidth;     % m - clear span of gutter
gutter.thickness = geom.gutter.horizThk;  % m - slab thickness
gutter.width = 1.0;                       % m - design width (per metre)

% Reinforcement details
gutter.barSize = 15;               % 15M bar
gutter.spacing = 0.300;            % m - 300mm spacing
gutter.cover = 0.040;              % m - 40mm clear cover (CSA A23.3-14, Table 10.6.1 for severe exposure)

% 15M bar properties (CSA A23.3)
gutter.Ab = 200;                   % mm^2 - area of 15M bar
gutter.db = 16;                    % mm - nominal diameter of 15M bar

% Calculate effective depth
gutter.d = gutter.thickness - gutter.cover - gutter.db/(2*1000);  % m

% Steel area per metre width
gutter.As = gutter.Ab * (gutter.width / gutter.spacing);  % mm^2/m

fprintf('Gutter Geometry:\n');
fprintf('  Span:                   %.3f m\n', gutter.span);
fprintf('  Thickness:              %.0f mm\n', gutter.thickness * 1000);
fprintf('  Effective depth (d):    %.0f mm\n', gutter.d * 1000);
fprintf('  Reinforcement:          15M @ %.0f mm\n', gutter.spacing * 1000);
fprintf('  As provided:            %.0f mm²/m\n', gutter.As);
fprintf('  Cover:                  %.0f mm\n\n', gutter.cover * 1000);

% Loading on gutter (simplified - self-weight + water weight)
gutter.DL = mat.gammaConcrete * gutter.thickness;  % kN/m² - dead load (self-weight)
gutter.waterDepthInGutter = 0.0;                 % m - assumed water depth in gutter (200mm)
gutter.LL = mat.gammaWater * gutter.waterDepthInGutter;  % kN/m² - live load (water)

% Load factors (CSA A23.3-14, Cl. 8.3.2)
gutter.alphaD = 1.25;              % dead load factor
gutter.alphaL = 1.5;               % live load factor

% Factored load per metre width
gutter.wf = (gutter.alphaD * gutter.DL + gutter.alphaL * gutter.LL) * gutter.width;  % kN/m

fprintf('Loading:\n');
fprintf('  Dead load (DL):         %.2f kN/m²\n', gutter.DL);
fprintf('  Live load (LL):         %.2f kN/m² (water in gutter)\n', gutter.LL);
fprintf('  Factored load (wf):     %.2f kN/m\n\n', gutter.wf);

%% --- MOMENT CHECK ---
fprintf('--- MOMENT CHECK (CSA A23.3-14, Cl. 10.1.7) ---\n');

% Maximum factored moment (simply supported or cantilever from stem)
% Assuming cantilever from stem support
gutter.Mf = gutter.wf * gutter.span^2 / 2;  % kN·m/m

fprintf('  Factored moment (Mf):   %.2f kN·m/m\n', gutter.Mf);

% Calculate balanced reinforcement ratio (CSA A23.3-14, Cl. 10.5.2)
gutter.beta1 = max(0.67, 0.97 - 0.0025 * gutter.fc_prime);  % stress block parameter
gutter.eps_y = gutter.fy / gutter.Es;                        % yield strain
gutter.eps_cu = 0.0035;                                       % ultimate concrete strain (CSA A23.3-14)
gutter.cb_balanced = gutter.eps_cu / (gutter.eps_cu + gutter.eps_y) * gutter.d;  % balanced neutral axis depth
gutter.rho_balanced = gutter.phi_c / gutter.phi_s * 0.85 * gutter.fc_prime * gutter.beta1 / gutter.fy * ...
                      gutter.cb_balanced / gutter.d;

% Actual reinforcement ratio
gutter.rho = gutter.As / (gutter.width * 1000 * gutter.d * 1000);  % convert to consistent units

fprintf('  Reinforcement ratio (ρ):         %.4f\n', gutter.rho);
fprintf('  Balanced ratio (ρ_balanced):     %.4f\n', gutter.rho_balanced);

if gutter.rho > gutter.rho_balanced
    fprintf('  WARNING: Over-reinforced section (compression failure)\n');
end

% Calculate factored moment resistance (CSA A23.3-14, Cl. 10.10)
% Assume rectangular section with tension reinforcement only
gutter.a = gutter.phi_s * gutter.As * gutter.fy / (gutter.phi_c * 0.85 * gutter.fc_prime * gutter.width * 1000);  % mm
gutter.c = gutter.a / gutter.beta1;  % mm - neutral axis depth

% Check if section is tension-controlled (c ≤ 0.375d per CSA)
if gutter.c <= 0.375 * gutter.d * 1000
    fprintf('  Section is tension-controlled (c/d = %.3f ≤ 0.375)\n', gutter.c / (gutter.d * 1000));
else
    fprintf('  WARNING: Section may not be tension-controlled (c/d = %.3f > 0.375)\n', gutter.c / (gutter.d * 1000));
end

gutter.Mr = gutter.phi_s * gutter.As * gutter.fy * (gutter.d * 1000 - gutter.a / 2) / 1e6;  % kN·m/m

fprintf('  Neutral axis depth (c):          %.1f mm\n', gutter.c);
fprintf('  Moment resistance (Mr):          %.2f kN·m/m\n', gutter.Mr);
fprintf('  Demand/Capacity (Mf/Mr):         %.3f\n', gutter.Mf / gutter.Mr);

if gutter.Mf <= gutter.Mr
    fprintf('  ✓ MOMENT CHECK: PASS (Mf ≤ Mr)\n\n');
else
    fprintf('  ✗ MOMENT CHECK: FAIL (Mf > Mr) - Increase reinforcement or thickness\n\n');
end

%% --- SHEAR CHECK ---
fprintf('--- SHEAR CHECK (CSA A23.3-14, Cl. 11.3) ---\n');

% Maximum factored shear (at support)
gutter.Vf = gutter.wf * gutter.span;  % kN/m

fprintf('  Factored shear (Vf):             %.2f kN/m\n', gutter.Vf);

% Concrete shear resistance (CSA A23.3-14, Cl. 11.3.4)
% For members without shear reinforcement
gutter.lambda = 1.0;                  % normal density concrete
gutter.beta = 0.21;                   % factor for simplified method
gutter.Vc = gutter.phi_c * gutter.lambda * gutter.beta * sqrt(gutter.fc_prime) * ...
            gutter.width * 1000 * gutter.d * 1000 / 1000;  % N -> kN

fprintf('  Concrete shear resistance (Vc): %.2f kN/m\n', gutter.Vc);
fprintf('  Demand/Capacity (Vf/Vc):         %.3f\n', gutter.Vf / gutter.Vc);

if gutter.Vf <= gutter.Vc
    fprintf('  ✓ SHEAR CHECK: PASS (Vf ≤ Vc)\n');
    fprintf('  No shear reinforcement required\n\n');
else
    fprintf('  ✗ SHEAR CHECK: FAIL (Vf > Vc)\n');
    fprintf('  Shear reinforcement required or increase thickness\n\n');
end

%% --- DEFLECTION CHECK ---
fprintf('--- DEFLECTION CHECK (CSA A23.3-14, Cl. 9.8) ---\n');

% Service load (unfactored)
gutter.ws = (gutter.DL + gutter.LL) * gutter.width;  % kN/m

% Calculate cracking moment (CSA A23.3-14, Cl. 9.8.4.2)
gutter.fr = 0.6 * gutter.lambda * sqrt(gutter.fc_prime);  % MPa - modulus of rupture
gutter.Ig = gutter.width * 1000 * (gutter.thickness * 1000)^3 / 12;  % mm^4 - gross moment of inertia
gutter.yt = gutter.thickness * 1000 / 2;  % mm - distance from neutral axis to extreme tension fibre
gutter.Mcr = gutter.fr * gutter.Ig / gutter.yt / 1e6;  % kN·m/m - cracking moment

fprintf('  Service load (ws):               %.2f kN/m\n', gutter.ws);
fprintf('  Cracking moment (Mcr):           %.2f kN·m/m\n', gutter.Mcr);

% Service moment (cantilever)
gutter.Ms = gutter.ws * gutter.span^2 / 2;  % kN·m/m

fprintf('  Service moment (Ms):             %.2f kN·m/m\n', gutter.Ms);

% Determine if section is cracked
if gutter.Ms > gutter.Mcr
    fprintf('  Section is CRACKED under service load\n');
    
    % Calculate cracked moment of inertia (transformed section)
    gutter.n = gutter.Es / gutter.Ec;  % modular ratio
    
    % Solve for neutral axis of cracked section (kd)
    % Quadratic equation: (b/2)(kd)^2 = n*As*(d - kd)
    gutter.k = sqrt(2 * gutter.rho * gutter.n + (gutter.rho * gutter.n)^2) - gutter.rho * gutter.n;
    gutter.kd = gutter.k * gutter.d * 1000;  % mm
    
    % Cracked moment of inertia
    gutter.Icr = gutter.width * 1000 * gutter.kd^3 / 3 + ...
                 gutter.n * gutter.As * (gutter.d * 1000 - gutter.kd)^2;  % mm^4
    
    % Effective moment of inertia (CSA A23.3-14, Cl. 9.8.4.2)
    gutter.Ie = (gutter.Mcr / gutter.Ms)^3 * gutter.Ig + ...
                (1 - (gutter.Mcr / gutter.Ms)^3) * gutter.Icr;
    
    fprintf('  Effective moment of inertia (Ie): %.2e mm^4\n', gutter.Ie);
else
    fprintf('  Section is UNCRACKED under service load\n');
    gutter.Ie = gutter.Ig;
    fprintf('  Using gross moment of inertia (Ig): %.2e mm^4\n', gutter.Ie);
end

% Calculate immediate deflection (cantilever)
gutter.delta_i = gutter.ws * (gutter.span * 1000)^4 / (8 * gutter.Ec * gutter.Ie);  % mm

fprintf('  Immediate deflection:            %.2f mm\n', gutter.delta_i);

% Long-term deflection (CSA A23.3-14, Cl. 9.8.4.3)
gutter.xi = 2.0;                   % time-dependent factor (5 years or more)
gutter.rho_prime = 0;              % compression reinforcement ratio (assumed zero)
gutter.lambda_delta = gutter.xi / (1 + 50 * gutter.rho_prime);  % deflection multiplier

gutter.delta_total = gutter.delta_i * (1 + gutter.lambda_delta);  % mm

fprintf('  Long-term multiplier (λ):        %.2f\n', gutter.lambda_delta);
fprintf('  Total deflection (δ_total):      %.2f mm\n', gutter.delta_total);

% Deflection limit (CSA A23.3-14, Table 9.3)
gutter.delta_limit = gutter.span * 1000 / 240;  % mm (for members supporting non-structural elements)

fprintf('  Deflection limit (L/240):        %.2f mm\n', gutter.delta_limit);
fprintf('  Demand/Limit:                    %.3f\n', gutter.delta_total / gutter.delta_limit);

if gutter.delta_total <= gutter.delta_limit
    fprintf('  ✓ DEFLECTION CHECK: PASS (δ ≤ L/240)\n\n');
else
    fprintf('  ✗ DEFLECTION CHECK: FAIL (δ > L/240) - Consider increasing thickness\n\n');
end

fprintf('--- END OF GUTTER STRUCTURAL ADEQUACY CHECK ---\n\n');

%% -------------------------------------------------------------------------
%  LOCAL FUNCTION(S)
% -------------------------------------------------------------------------

function components = addComponent(components, name, weight, leverArm, notes)
% addComponent helper to append gravity load data to the component list
    component.Name     = name;
    component.Weight   = weight;
    component.LeverArm = leverArm;
    component.Moment   = weight * leverArm;
    component.Notes    = notes;
    components = [components; component];
end

