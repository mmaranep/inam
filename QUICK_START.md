# RC Beam Design CSA - Quick Start Guide

## What's Fixed

✅ **Unit conversion issues** - Moment now correctly calculated as ~250 kN·m (not 390,000!)
✅ **Iteration removed** - Direct quadratic formula (no convergence problems)
✅ **Validation added** - Comprehensive error handling
✅ **Production ready** - All checks passing

## Quick Test

### Run the main script:
```matlab
% In MATLAB or Octave
RC_Beam_Design_CSA
```

### Expected key outputs:
```
Factored load:              56.25 kN/m
Maximum moment (Mf):        253.12 kN·m  ← CORRECT (not 390,000!)
Maximum shear (Vf):         168.75 kN
Required steel:             1540 mm²
✓ MOMENT CHECK: PASS
✓ SHEAR CHECK: PASS
✓ DEFLECTION CHECK: PASS
```

## How It Works

### 1. Direct Quadratic Formula (No Iteration!)
```
A·As² + B·As + C = 0

where:
  A = (φs·fy)² / (2·φc·0.85·f'c·b)
  B = -φs·fy·d  
  C = Mf (factored moment in N·mm)

Solution: As = (-B ± √(B²-4AC)) / (2A)
```

### 2. Unit Conversions (Made Explicit)
```matlab
% Define unit converters
units.mm_to_m = 1e-3;
units.kNm_to_Nmm = 1e6;

% Use them consistently
Mf = loads.wf * (beam.span * units.mm_to_m)^2 / 8;  % kN·m
```

### 3. Key Formulas
```
Simply supported beam with uniform load w (kN/m), span L (m):
  M_max = w·L²/8  (kN·m)
  V_max = w·L/2   (kN)

Important: Always convert span to meters when using kN/m loads!
```

## Modify Parameters

Edit lines 31-62 in RC_Beam_Design_CSA.m:

```matlab
% Beam geometry
beam.span = 6000;      % mm - try 5000, 7000, 8000
beam.width = 300;      % mm - try 250, 350, 400
beam.depth = 600;      % mm - try 500, 700, 800

% Loading
loads.DL = 15.0;       % kN/m - adjust as needed
loads.LL = 25.0;       % kN/m - adjust as needed

% Materials
mat.fc_prime = 30;     % MPa - try 25, 35, 40
mat.fy_flexural = 400; % MPa - try 300, 500
```

## Verification Scripts

### Python verification:
```bash
python3 verify_units.py
```
Shows step-by-step calculation with both unit methods.

### Simple MATLAB/Octave test:
```matlab
test_RC_simple
```
Quick validation without plotting.

## What Changed

| Aspect | Before | After |
|--------|--------|-------|
| Method | Iteration | Direct formula |
| Moment | 390,000 kN·m? | 253.12 kN·m ✓ |
| Convergence | Issues | None |
| Units | Unclear | Explicit |
| Errors | Generic | Helpful |

## Common Issues

### "Negative discriminant" error
**Meaning:** Section too small to carry the moment
**Solutions:**
1. Increase beam depth
2. Increase beam width  
3. Increase concrete strength
4. Use compression reinforcement

### "Bar spacing too small" warning
**Meaning:** Bars are too close together
**Solutions:**
1. Increase beam width
2. Use smaller bar size
3. Use fewer bars (if strength allows)

### Unrealistic values
**Check:**
- Moment should be 100-500 kN·m for typical beams
- Shear should be 50-300 kN for typical beams
- Steel area should be 500-5000 mm² for typical beams

## Output Files

The script generates:
1. Console output with detailed calculations
2. Figure with 6 subplots:
   - Beam elevation with reinforcement
   - Factored moment diagram
   - Service moment diagram
   - Shear force diagram
   - Stirrup spacing layout
3. Summary table with pass/fail status

## Documentation

- `RC_BEAM_DESIGN_README.md` - Comprehensive guide
- `FIXES_SUMMARY.md` - What was fixed and why
- `VERIFICATION_RESULTS.md` - Detailed test results
- `QUICK_START.md` - This file

## Need Help?

1. Check if values are reasonable (moment ~250 kN·m for default case)
2. Review error messages (they suggest solutions)
3. Run verification scripts to compare calculations
4. Check input parameters are in correct units (mm, kN/m, MPa)

## Success Indicators

✓ Moment between 100-500 kN·m (for typical beams)
✓ Shear between 50-300 kN (for typical beams)
✓ Steel area between 500-5000 mm² (for typical beams)
✓ Discriminant is positive
✓ All checks show PASS
✓ No error messages
✓ Diagrams generate successfully

**If you see these, the script is working correctly!**
