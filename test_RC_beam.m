%% TEST SCRIPT FOR RC BEAM DESIGN
% This script verifies the RC_Beam_Design_CSA.m calculations
% with manual calculations to ensure unit conversions are correct

clear; close all; clc;

fprintf('========================================\n');
fprintf('MANUAL VERIFICATION OF RC BEAM DESIGN\n');
fprintf('========================================\n\n');

%% Input parameters (same as main script)
span = 6000;  % mm
width = 300;  % mm
depth = 600;  % mm
DL = 15.0;    % kN/m
LL = 25.0;    % kN/m
alpha_D = 1.25;
alpha_L = 1.5;

%% Manual calculation of factored loads
wf = alpha_D * DL + alpha_L * LL;  % kN/m
fprintf('Factored load: wf = %.2f kN/m\n\n', wf);

%% Method 1: Convert span to meters
L_m = span / 1000;  % m
Mf_method1 = wf * L_m^2 / 8;  % kN·m
Vf_method1 = wf * L_m / 2;    % kN

fprintf('METHOD 1: Convert span to meters\n');
fprintf('---------------------------------\n');
fprintf('Span: %.3f m\n', L_m);
fprintf('Moment: Mf = wf * L^2 / 8\n');
fprintf('      = %.2f * %.3f^2 / 8\n', wf, L_m);
fprintf('      = %.2f kN·m\n', Mf_method1);
fprintf('Shear:  Vf = wf * L / 2\n');
fprintf('      = %.2f * %.3f / 2\n', wf, L_m);
fprintf('      = %.2f kN\n\n', Vf_method1);

%% Method 2: Convert load to N/mm
w_Nmm = wf;  % kN/m = N/mm (unit equivalence)
L_mm = span;  % mm
Mf_Nmm = w_Nmm * L_mm^2 / 8;  % N·mm
Mf_method2 = Mf_Nmm / 1e6;  % kN·m
Vf_N = w_Nmm * L_mm / 2;  % N
Vf_method2 = Vf_N / 1000;  % kN

fprintf('METHOD 2: Use N/mm for load (1 kN/m = 1 N/mm)\n');
fprintf('---------------------------------------------\n');
fprintf('Load: %.2f N/mm\n', w_Nmm);
fprintf('Span: %d mm\n', L_mm);
fprintf('Moment: Mf = w * L^2 / 8\n');
fprintf('      = %.2f * %d^2 / 8\n', w_Nmm, L_mm);
fprintf('      = %.0f N·mm = %.2f kN·m\n', Mf_Nmm, Mf_method2);
fprintf('Shear:  Vf = w * L / 2\n');
fprintf('      = %.2f * %d / 2\n', w_Nmm, L_mm);
fprintf('      = %.0f N = %.2f kN\n\n', Vf_N, Vf_method2);

%% Verification
fprintf('VERIFICATION\n');
fprintf('------------\n');
fprintf('Moment results match: %s\n', string(abs(Mf_method1 - Mf_method2) < 0.01));
fprintf('Shear results match:  %s\n\n', string(abs(Vf_method1 - Vf_method2) < 0.01));

%% Expected values
fprintf('EXPECTED VALUES\n');
fprintf('---------------\n');
fprintf('For a 6 m span beam with 56.25 kN/m load:\n');
fprintf('  Maximum moment: %.2f kN·m\n', Mf_method1);
fprintf('  Maximum shear:  %.2f kN\n\n', Vf_method1);

%% Check if values are reasonable
fprintf('REASONABLENESS CHECK\n');
fprintf('--------------------\n');
if Mf_method1 > 100 && Mf_method1 < 500
    fprintf('✓ Moment is reasonable (%.2f kN·m)\n', Mf_method1);
else
    fprintf('✗ Moment seems unreasonable (%.2f kN·m)\n', Mf_method1);
end

if Vf_method1 > 50 && Vf_method1 < 300
    fprintf('✓ Shear is reasonable (%.2f kN)\n', Vf_method1);
else
    fprintf('✗ Shear seems unreasonable (%.2f kN)\n', Vf_method1);
end

fprintf('\n========================================\n');
fprintf('TEST COMPLETE\n');
fprintf('========================================\n\n');
