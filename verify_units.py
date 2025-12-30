#!/usr/bin/env python3
"""
Verification script for RC beam design calculations
Tests unit conversions and validates formulas
"""

import math

print("="*60)
print("RC BEAM DESIGN UNIT VERIFICATION")
print("="*60)
print()

# Input parameters
span_mm = 6000  # mm
width_mm = 300  # mm
depth_mm = 600  # mm
DL = 15.0  # kN/m
LL = 25.0  # kN/m
alpha_D = 1.25
alpha_L = 1.5

# Material properties
fc_prime = 30  # MPa
fy = 400  # MPa
phi_c = 0.65
phi_s = 0.85
cover = 40  # mm

# Factored load
wf = alpha_D * DL + alpha_L * LL  # kN/m
ws = DL + LL  # kN/m (service)

print(f"INPUT PARAMETERS")
print(f"----------------")
print(f"Span:            {span_mm} mm = {span_mm/1000} m")
print(f"Width:           {width_mm} mm")
print(f"Depth:           {depth_mm} mm")
print(f"Dead load:       {DL} kN/m")
print(f"Live load:       {LL} kN/m")
print(f"Factored load:   {wf} kN/m")
print()

# Method 1: Convert span to meters
print("METHOD 1: Convert span to meters")
print("-"*60)
L_m = span_mm / 1000  # m
Mf_kNm = wf * L_m**2 / 8  # kN·m
Vf_kN = wf * L_m / 2  # kN

print(f"Span in meters:  {L_m} m")
print(f"Moment formula:  Mf = wf × L² / 8")
print(f"               = {wf} × {L_m}² / 8")
print(f"               = {wf} × {L_m**2} / 8")
print(f"               = {Mf_kNm:.2f} kN·m")
print()
print(f"Shear formula:   Vf = wf × L / 2")
print(f"               = {wf} × {L_m} / 2")
print(f"               = {Vf_kN:.2f} kN")
print()

# Method 2: Use N/mm (note: 1 kN/m = 1 N/mm)
print("METHOD 2: Use N/mm for load")
print("-"*60)
print("Note: 1 kN/m = 1 kN / 1000 mm = 0.001 kN/mm = 1 N/mm")
w_Nmm = wf  # N/mm (numerically equal to kN/m)
L_mm = span_mm  # mm
Mf_Nmm = w_Nmm * L_mm**2 / 8  # N·mm
Mf_kNm_method2 = Mf_Nmm / 1e6  # kN·m
Vf_N = w_Nmm * L_mm / 2  # N
Vf_kN_method2 = Vf_N / 1000  # kN

print(f"Load in N/mm:    {w_Nmm} N/mm")
print(f"Span in mm:      {L_mm} mm")
print(f"Moment:          {Mf_Nmm:,.0f} N·mm = {Mf_kNm_method2:.2f} kN·m")
print(f"Shear:           {Vf_N:,.0f} N = {Vf_kN_method2:.2f} kN")
print()

# Verification
print("VERIFICATION")
print("-"*60)
moment_match = abs(Mf_kNm - Mf_kNm_method2) < 0.01
shear_match = abs(Vf_kN - Vf_kN_method2) < 0.01
print(f"Moment results match:  {'✓ YES' if moment_match else '✗ NO'}")
print(f"  Method 1: {Mf_kNm:.2f} kN·m")
print(f"  Method 2: {Mf_kNm_method2:.2f} kN·m")
print(f"Shear results match:   {'✓ YES' if shear_match else '✗ NO'}")
print(f"  Method 1: {Vf_kN:.2f} kN")
print(f"  Method 2: {Vf_kN_method2:.2f} kN")
print()

# Reasonableness check
print("REASONABLENESS CHECK")
print("-"*60)
moment_ok = 100 < Mf_kNm < 500
shear_ok = 50 < Vf_kN < 300
print(f"Moment: {Mf_kNm:.2f} kN·m  {'✓ Reasonable' if moment_ok else '✗ Unreasonable'}")
print(f"Shear:  {Vf_kN:.2f} kN    {'✓ Reasonable' if shear_ok else '✗ Unreasonable'}")
print()

# Steel area calculation using quadratic formula
print("STEEL AREA CALCULATION")
print("-"*60)

# Assume effective depth
db_stirrup = 11.3  # mm (10M)
db_flex = 25.2  # mm (25M)
d = depth_mm - cover - db_stirrup - db_flex/2
print(f"Effective depth: d = {d:.1f} mm")
print()

# Convert moment to N·mm
Mf_Nmm = Mf_kNm * 1e6
print(f"Factored moment: {Mf_kNm:.2f} kN·m = {Mf_Nmm:,.0f} N·mm")
print()

# Quadratic coefficients: A*As² + B*As + C = 0
A = (phi_s * fy)**2 / (2 * phi_c * 0.85 * fc_prime * width_mm)
B = -phi_s * fy * d
C = Mf_Nmm

print(f"Quadratic equation: A·As² + B·As + C = 0")
print(f"  A = {A:.6e}")
print(f"  B = {B:.6e}")
print(f"  C = {C:.6e}")
print()

# Discriminant
discriminant = B**2 - 4*A*C
print(f"Discriminant: {discriminant:.6e}")

if discriminant < 0:
    print("✗ ERROR: Negative discriminant - section inadequate")
else:
    print("✓ Discriminant is positive - solutions exist")
    
    # Solve
    As1 = (-B + math.sqrt(discriminant)) / (2*A)
    As2 = (-B - math.sqrt(discriminant)) / (2*A)
    
    print(f"  As₁ = {As1:.1f} mm²")
    print(f"  As₂ = {As2:.1f} mm²")
    
    if As1 > 0 and As2 > 0:
        As_req = min(As1, As2)
        print(f"  Selected (smaller positive): {As_req:.1f} mm²")
    elif As1 > 0:
        As_req = As1
        print(f"  Selected (As₁, only positive): {As_req:.1f} mm²")
    elif As2 > 0:
        As_req = As2
        print(f"  Selected (As₂, only positive): {As_req:.1f} mm²")
    else:
        print("  ✗ ERROR: No positive solution")
        As_req = None
    
    if As_req:
        # Verify the solution
        a = phi_s * As_req * fy / (phi_c * 0.85 * fc_prime * width_mm)
        Mr_Nmm = phi_s * As_req * fy * (d - a/2)
        Mr_kNm = Mr_Nmm / 1e6
        
        print()
        print(f"VERIFICATION OF SOLUTION")
        print(f"  Stress block depth: a = {a:.2f} mm")
        print(f"  Moment resistance:  Mr = {Mr_kNm:.2f} kN·m")
        print(f"  Factored moment:    Mf = {Mf_kNm:.2f} kN·m")
        print(f"  Ratio Mf/Mr:        {Mf_kNm/Mr_kNm:.3f}")
        print(f"  Status: {'✓ PASS' if Mf_kNm <= Mr_kNm else '✗ FAIL'}")

print()
print("="*60)
print("VERIFICATION COMPLETE")
print("="*60)
