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
geom.heelWidth  = 1270 * mm;
geom.toeWidth   = 610  * mm;
geom.stemThk    = 0.30;         % assumed constant thickness
geom.baseThk    = 0.55;         % assumed footing thickness
geom.keyWidth   = 280  * mm;
geom.keyDepth   = 370  * mm;

geom.gutter.offsetBelowTop = 400 * mm;
geom.gutter.horizWidth     = 800 * mm;
geom.gutter.horizThk       = 150 * mm;
geom.gutter.vertHeight     = 290 * mm;
geom.gutter.vertThk        = 150 * mm;

geom.roof.thk           = 225 * mm;
geom.roof.tributarySpan = 2.75; % m

% Optional heel floor slab weight (set thickness > 0 to include)
geom.floorSlabThk   = 0.0;      % m
geom.floorSlabWidth = geom.heelWidth;

% Hydrostatic loading
loads.waterDepth = 1600 * mm;   % measured from heel slab surface

% Material properties
mat.gammaConcrete = 24.0;       % kN/m^3
mat.gammaWater    = 9.81;       % kN/m^3
mat.gammaSoil     = 19.0;       % kN/m^3
mat.phi           = deg2rad(32);
mat.delta         = deg2rad(15);

% CSA target factors of safety
criteria.overturning = 1.50;
criteria.sliding     = 1.50;
criteria.bearing     = 3.00;

%% -------------------------------------------------------------------------
%  DERIVED GEOMETRY
% -------------------------------------------------------------------------

geom.baseWidth        = geom.toeWidth + geom.heelWidth;
geom.stemCL           = geom.toeWidth;
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
%  HYDROSTATIC & UPLIFT ACTIONS
% -------------------------------------------------------------------------

hydro.thrust = 0.5 * mat.gammaWater * loads.waterDepth^2;            % kN/m
hydro.arm    = geom.baseTopElevation + loads.waterDepth / 3;          % m
hydro.moment = hydro.thrust * hydro.arm;                              % kN·m/m

uplift.p_max = mat.gammaWater * loads.waterDepth;                     % kPa
uplift.force = 0.5 * uplift.p_max * geom.baseWidth;                   % kN/m
uplift.arm   = 2 * geom.baseWidth / 3;                                % m
uplift.moment = uplift.force * uplift.arm;                            % kN·m/m

V_net = W_total - uplift.force;                                       % kN/m (downwards positive)

%% -------------------------------------------------------------------------
%  STABILITY CHECKS (CSA)
% -------------------------------------------------------------------------

FS_overturning = (M_resisting - uplift.moment) / hydro.moment;

friction_resistance = V_net * tan(mat.delta);
Kp = tan(pi/4 + mat.phi/2)^2;
passive_resistance = 0.5 * Kp * mat.gammaSoil * geom.keyDepth^2 * geom.keyWidth;

FS_sliding = (friction_resistance + passive_resistance) / hydro.thrust;

M_net_about_toe = M_resisting - uplift.moment - hydro.moment;

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
q_ult  = 0.5 * mat.gammaSoil * geom.baseWidth * Ngamma;     % kPa
q_allow = q_ult / criteria.bearing;

if isnan(q_max) || q_max <= 0
    FS_bearing = NaN;
else
    FS_bearing = q_ult / q_max;
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
FS_required = [criteria.overturning; criteria.sliding; criteria.bearing];
statusStr   = repmat("PASS", size(FS_values));
statusStr(FS_values < FS_required) = "CHECK";

resultsTable = table(["Overturning"; "Sliding"; "Bearing"], FS_values, FS_required, statusStr, ...
    'VariableNames', {'Check', 'FS_Computed', 'FS_Required', 'Status'});

bearingTable = table(q_toe, q_heel, q_avg, q_min, q_max, q_allow, within_middle_third, ...
    'VariableNames', {'q_toe_kPa', 'q_heel_kPa', 'q_avg_kPa', 'q_min_kPa', 'q_max_kPa', ...
                      'q_allowable_kPa', 'WithinMiddleThird'});

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
fprintf('Uplift force             : %.2f kN/m\n', uplift.force);
fprintf('Net vertical load (W-U)  : %.2f kN/m\n\n', V_net);

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
if FS_bearing < criteria.bearing
    fprintf('WARNING: Bearing FS below CSA minimum (%.2f < %.2f).\n', FS_bearing, criteria.bearing);
end
if ~within_middle_third
    fprintf('WARNING: Resultant lies outside the middle third. Heel uplift is expected.\n');
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
yline(q_allow, 'Color', [0.3 0.6 0.3], 'LineStyle', '--', 'LineWidth', 1.2);

if ~isnan(q_min) && q_min < 0
    text(geom.baseWidth * 0.8, min(q_plot) - 10, 'Heel uplift (q_{min} < 0)', ...
        'Color', [0.6 0.1 0.1], 'FontWeight', 'bold');
end

xlabel('Distance from toe (m)');
ylabel('Bearing pressure (kPa)');
title('Base Contact Pressure');
xlim([0, geom.baseWidth + 0.4]);
ylim([min([0, q_plot, q_allow], [], 'omitnan') - 20, max([q_plot, q_allow], [], 'omitnan') + 20]);
legend({'Contact pressure', 'Zero', 'Allowable (FS=3.0)'}, 'Location', 'best');

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

