# Pipe Support Design - Quick Start Guide

## 🚀 Quick Start (3 Easy Steps)

### Step 1: Open MATLAB
Navigate to the workspace folder:
```matlab
cd /path/to/workspace
```

### Step 2: Run the Analysis
Choose one of these options:

#### Option A: Complete Analysis (Recommended)
```matlab
run_complete_analysis
```
This runs both the console calculation AND creates the Excel file.

#### Option B: Console Only
```matlab
pipe_support_design
```
View results in MATLAB console only.

#### Option C: Excel Export Only
```matlab
pipe_support_to_excel
```
Creates Excel file without console output.

### Step 3: Review Results
Open `Pipe_Support_Design_Calculation.xlsx` and check the **Summary** sheet.

---

## 📋 What Gets Calculated?

### Pipe System
- **Pipe**: 12" Schedule 10 Stainless Steel 304L
- **Contents**: Full of water
- **Weight**: ~375 kg total (pipe + water)

### Support Structure
- **Post**: HSS 3×3×0.25, Height: 1000mm
- **Top Plate**: 230×100×6mm curved plate
- **Top Bolts**: 2× 7/8" bolts
- **Weld**: 6mm fillet weld
- **Base Plate**: 8×8×0.25 inch
- **Anchors**: 4× Hilti 1/2" drop-in anchors

### Design Checks Performed
✓ Post axial stress  
✓ Post bending stress  
✓ Post deflection  
✓ Slenderness ratio  
✓ Weld stress  
✓ Bolt shear & tension  
✓ Plate bearing  
✓ Concrete bearing  
✓ Anchor bolts  
✓ Base plate bending  

---

## 🎯 Understanding the Results

### In the Console
Look for lines like:
```
Utilization ratio: 45.2%
Status: OK ✓
```

- **< 100%** = Good (component is adequate)
- **≥ 100%** = Problem (component needs revision)

### In the Excel File
Go to **Sheet 9: Summary** and find:
```
OVERALL DESIGN STATUS: ACCEPTABLE ✓
```

If you see **REQUIRES REVISION**, check individual component utilizations.

---

## ⚙️ Customizing Parameters

To modify the design, edit these files:

### For Console Output:
Edit `pipe_support_design.m` - Lines 15-60

### For Excel Export:
Edit `pipe_support_to_excel.m` - Lines 15-60

### Common Changes:

```matlab
% Change post height
post_height = 1200;  % mm (was 1000)

% Change pipe span
pipe_length_supported = 1500;  % mm (was 1000)

% Change pipe size
pipe_OD = 10.75;  % For 10" pipe (was 12.75)
pipe_wall_thick = 0.165;

% Change post size
post_size = 4.0;  % 4" square (was 3")

% Change lateral load assumption
H_lateral = 0.10 * P_dead;  % 10% of vertical (was 5%)
```

After editing, save and re-run the script.

---

## 📊 Excel File Structure

| Sheet # | Name | Contents |
|---------|------|----------|
| 1 | Input Parameters | All design inputs |
| 2 | Section Properties | Geometric properties |
| 3 | Load Calculations | Weights and forces |
| 4 | Post Analysis | Post stress & deflection |
| 5 | Weld Check | Weld stress analysis |
| 6 | Top Bolt Check | Bolt capacity |
| 7 | Top Plate Check | Plate bearing |
| 8 | Base Plate & Anchors | Foundation design |
| 9 | **Summary** | **Overall results** ⭐ |

👉 **Always check Sheet 9 first!**

---

## ⚠️ Important Notes

### Load Assumptions
- **Dead load only** (no live load, no seismic, no wind)
- **Lateral load**: 5% of vertical (adjustable)
- **Impact factor**: None included
- **Dynamic effects**: Not considered

### Design Standards
- Based on **AISC/CSA** steel design
- **ACI 318** for concrete
- Conservative assumptions used
- Professional review required

### Before Construction
✅ Have a licensed engineer review the design  
✅ Verify local building code requirements  
✅ Check seismic and wind load requirements  
✅ Confirm anchor embedment depths  
✅ Verify concrete strength in field  
✅ Check edge distances for anchors  

---

## 🔧 Troubleshooting

### Problem: "File already open" error
**Solution**: Close the Excel file and re-run the script

### Problem: Calculation shows "FAIL"
**Solution**: 
1. Increase component size (post, plate, bolts)
2. Reduce span length
3. Check if loads are realistic

### Problem: MATLAB can't find script
**Solution**: 
```matlab
cd /path/to/workspace
ls  % Verify files are present
```

### Problem: Excel file not created
**Solution**: Check MATLAB version (needs R2019a+)

---

## 📁 File Organization

```
workspace/
├── pipe_support_design.m           ← Console calculation
├── pipe_support_to_excel.m         ← Excel export
├── run_complete_analysis.m         ← Run both (recommended)
├── README_pipe_support.md          ← Detailed calculation guide
├── README_excel_export.md          ← Excel export guide
├── QUICK_START.md                  ← This file
└── Pipe_Support_Design_Calculation.xlsx  ← Output (generated)
```

---

## 💡 Tips

1. **Always run `run_complete_analysis`** to get both console and Excel output
2. **Check the Summary sheet first** in Excel for quick status
3. **Save different versions** when trying different parameters
4. **Add notes** in the Excel file for future reference
5. **Print the Summary sheet** for field use

---

## 📞 Need Help?

1. Review the detailed README files
2. Check inline comments in the MATLAB code
3. Verify input parameters are reasonable
4. Consult with a licensed structural engineer

---

## ✅ Quick Checklist

Before finalizing your design:

- [ ] Ran the complete analysis
- [ ] Checked Excel Summary sheet
- [ ] All components show "OK"
- [ ] Overall status is "ACCEPTABLE"
- [ ] Verified load assumptions
- [ ] Confirmed material grades
- [ ] Checked local code requirements
- [ ] Had engineer review the design
- [ ] Verified field conditions match assumptions

---

**Version**: 1.0  
**Last Updated**: November 6, 2025  

**Ready to start? Run this command:**
```matlab
run_complete_analysis
```
