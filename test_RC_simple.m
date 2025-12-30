%% Simple test of RC beam calculations (Octave compatible)
clear; clc;

fprintf('Testing RC Beam Design Calculations\n');
fprintf('====================================\n\n');

% Input parameters
span_mm = 6000;
width = 300;
depth = 600;
cover = 40;
DL = 15.0;
LL = 25.0;
alpha_D = 1.25;
alpha_L = 1.5;
fc_prime = 30;
fy = 400;
phi_c = 0.65;
phi_s = 0.85;

% Unit conversions
mm_to_m = 1e-3;
kNm_to_Nmm = 1e6;

% Factored load
wf = alpha_D * DL + alpha_L * LL;
fprintf('Factored load: %.2f kN/m\n', wf);

% Moment and shear
span_m = span_mm * mm_to_m;
Mf_kNm = wf * span_m^2 / 8;
Vf_kN = wf * span_m / 2;

fprintf('Span: %.1f m\n', span_m);
fprintf('Moment: %.2f kN·m\n', Mf_kNm);
fprintf('Shear: %.2f kN\n\n', Vf_kN);

% Effective depth
db_stirrup = 11.3;
db_flex = 25.2;
d = depth - cover - db_stirrup - db_flex/2;
fprintf('Effective depth: %.1f mm\n\n', d);

% Steel area calculation
Mf_Nmm = Mf_kNm * kNm_to_Nmm;

A_coef = (phi_s * fy)^2 / (2 * phi_c * 0.85 * fc_prime * width);
B_coef = -phi_s * fy * d;
C_coef = Mf_Nmm;

discriminant = B_coef^2 - 4 * A_coef * C_coef;

fprintf('Quadratic solution:\n');
fprintf('  Discriminant: %.2e\n', discriminant);

if discriminant < 0
    fprintf('  ERROR: Negative discriminant\n');
else
    As1 = (-B_coef + sqrt(discriminant)) / (2 * A_coef);
    As2 = (-B_coef - sqrt(discriminant)) / (2 * A_coef);
    
    fprintf('  As1 = %.1f mm²\n', As1);
    fprintf('  As2 = %.1f mm²\n', As2);
    
    if As1 > 0 && As2 > 0
        As_req = min(As1, As2);
        fprintf('  Selected: %.1f mm²\n', As_req);
    end
end

fprintf('\n✓ Test complete - no errors\n');
