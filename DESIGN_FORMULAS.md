# Retaining Wall Design Formulas - Quick Reference

## 1. EARTH PRESSURE THEORY

### Active Earth Pressure Coefficient (Coulomb's Theory)

```
Ka = cos?(?) / [cos(?) ? (1 + ?[sin(?+?)?sin(?)/cos(?)])?]
```

Where:
- ? = soil internal friction angle
- ? = wall-soil friction angle (typically 2?/3)

### Lateral Earth Pressure Forces

```
Pa(soil) = ? ? Ka ? ?soil ? H?           [kN/m]
Pa(surcharge) = Ka ? q ? H               [kN/m]
Ph = Pa(total)                            [kN/m, horizontal]
Pv = Pa(total) ? tan(?)                   [kN/m, vertical]
```

Where:
- ?soil = unit weight of soil
- H = wall height
- q = surcharge pressure

### Pressure Distribution

```
p(z) = Ka ? ?soil ? z + Ka ? q
```

At depth z below surface.

---

## 2. STABILITY ANALYSIS

### Overturning Stability

```
FS(overturning) = ?MR / ?MO ? 2.0

MR = Resisting moments about toe
   = W1?x1 + W2?x2 + ... + Wn?xn

MO = Overturning moment about toe
   = Ph ? H/3
```

### Sliding Stability

```
FS(sliding) = FR / FD ? 1.5

FR = W ? tan(?f) + c ? B + Pp
FD = Ph

where:
  Pp = ? ? Kp ? ? ? D? (passive resistance, often neglected)
  Kp = tan?(45? + ?f/2)
```

### Bearing Capacity

**Resultant Location:**
```
x? = ?MR / ?W
e = B/2 - x?
```

**Bearing Pressure:**

If |e| ? B/6 (within middle third):
```
qmax = (W/B) ? (1 + 6e/B)
qmin = (W/B) ? (1 - 6e/B)
```

If |e| > B/6 (outside middle third):
```
L' = 3(B/2 - e)  [effective length]
qmax = 2W/L'
qmin = 0
```

**Ultimate Bearing Capacity (Terzaghi):**
```
qult = c?Nc + ??B'?N?/2 + ??Df?Nq

where:
  B' = B - 2e (effective width)
  Nq = e^(??tan ?) ? tan?(45? + ?/2)
  Nc = (Nq - 1) / tan(?)
  N? = 2(Nq - 1) ? tan(?)
```

**Factor of Safety:**
```
FS(bearing) = qult / qmax ? 3.0
```

---

## 3. STRUCTURAL DESIGN (CSA A23.3-19)

### Load Factors (Ultimate Limit State)

```
Factored Load = ?D?D + ?L?L + ?E?E

where:
  ?D = 1.25 (dead load)
  ?L = 1.5 (live load)
  ?E = 1.0 (earth pressure, used with ?D)
```

### Material Properties

**Concrete Stress Block Parameters:**
```
?1 = 0.85 - 0.0015?f'c  (but ? 0.67)
?1 = 0.97 - 0.0025?f'c  (but ? 0.67)
```

**Resistance Factors:**
```
?c = 0.65  (concrete)
?s = 0.85  (steel)
```

### Flexural Design

**Required Steel Ratio:**
```
Mf = ?c ? ?1 ? f'c ? b ? d? ? ? ? (1 - 0.59?)

? = 1 - ?(1 - 2Mf / (?c ? ?1 ? f'c ? b ? d?))

? = ? ? ?1 ? ?c ? f'c / (?s ? fy)

As = ? ? b ? d
```

**Minimum Reinforcement (CSA Cl. 10.5):**
```
?min = max(0.002, 1.4/fy)
As,min = ?min ? b ? h
```

**Maximum Reinforcement (Balanced Condition):**
```
?c = 0.0035  (concrete ultimate strain)
?y = fy/Es   (steel yield strain)

cb = ?c / (?c + ?y) ? d
ab = ?1 ? cb

?max = ?1 ? ?c ? f'c ? ab / (?s ? fy ? d)
```

**Moment Resistance:**
```
a = As ? fy / (?1 ? f'c ? b)
Mr = ?s ? As ? fy ? (d - a/2)
```

### Shear Design (CSA Cl. 11.3)

**Concrete Shear Resistance (Simplified Method):**
```
Vc = ?c ? ? ? ? ? ?(f'c) ? b ? d

where:
  ? = 1.0 (normal density concrete)
  ? = 0.21 (simplified method, conservative)
```

**Shear Check:**
```
If Vf ? Vc:  No shear reinforcement required
If Vf > Vc:  Provide stirrups
```

**Required Shear Reinforcement:**
```
Vs = Vf/?s - Vc
Av/s = Vs / (fy ? dv)

where dv = max(0.9d, 0.72h)
```

### Development Length (CSA Cl. 12.2)

**Basic Development Length:**
```
ld = k1 ? k2 ? k3 ? k4 ? fy ? db / (dcs ? ?(f'c))

where:
  k1 = 1.0 (bar location factor)
  k2 = 1.0 (normal density concrete)
  k3 = 1.0 (bar size factor)
  k4 = 0.8 (deformed bars)
  dcs = 1.0 (simplified)
  
Simplified: ld ? 0.45 ? fy ? db / ?(f'c)
```

### Crack Control (CSA Cl. 10.6)

**Crack Width Parameter:**
```
z = fs ? ??(dc ? A)  ? 30,000 N/mm

where:
  fs = steel stress at service loads
  dc = cover to center of bar
  A = effective tension area per bar
```

### Temperature & Shrinkage Reinforcement (CSA Cl. 7.8)

```
As,temp = 0.002 ? Ag

For 1-meter width:
As,temp = 0.002 ? t ? 1000  [mm?/m]

Typical: 10M @ 300 mm c/c
```

---

## 4. DESIGN MOMENTS AND SHEARS

### Stem (Cantilever from Base)

**Moment at Base:**
```
Mstem = (1/6) ? ?soil ? Ka ? H? + (1/2) ? q ? Ka ? H?
```

**Shear at Base:**
```
Vstem = (1/2) ? ?soil ? Ka ? H? + q ? Ka ? H
```

### Toe Slab

**Moment at Face of Stem:**
```
Mtoe = ?????? [q(x) ? x] dx - Wself ? L?toe/2

where q(x) varies linearly from toe to stem
```

**Shear at Distance d from Face:**
```
Vtoe = ??????? q(x) dx - Wself ? Lcrit
```

### Heel Slab

**Moment at Face of Stem:**
```
Mheel = Wsoil ? L/2 - ??????? [q(x) ? x] dx - Wself ? L?/2

(Note: Soil weight creates positive moment, 
       pressure creates negative moment)
```

**Shear at Distance d from Face:**
```
Vheel = Wsoil - ??????? q(x) dx - Wself ? Lcrit
```

---

## 5. TYPICAL DESIGN VALUES

### Soil Properties

| Soil Type | ? (kN/m?) | ? (deg) | c (kPa) | Ka |
|-----------|-----------|---------|---------|-----|
| Loose sand | 16-18 | 28-32 | 0 | 0.30-0.35 |
| Dense sand | 18-20 | 32-38 | 0 | 0.25-0.30 |
| Sandy gravel | 18-20 | 34-40 | 0 | 0.23-0.27 |
| Silty sand | 17-20 | 26-32 | 0-5 | 0.30-0.36 |
| Stiff clay | 18-20 | 20-28 | 20-50 | 0.36-0.45 |

### Concrete Properties

| f'c (MPa) | Ec (MPa) | ?1 | ?1 |
|-----------|----------|----|----|
| 25 | 24,800 | 0.813 | 0.908 |
| 30 | 26,800 | 0.805 | 0.895 |
| 35 | 28,600 | 0.798 | 0.883 |
| 40 | 30,100 | 0.790 | 0.870 |

### Steel Properties

| Grade | fy (MPa) | Es (MPa) | ?y |
|-------|----------|----------|----|
| 400 | 400 | 200,000 | 0.002 |
| 500 | 500 | 200,000 | 0.0025 |

---

## 6. DESIGN AIDS

### Preliminary Sizing Rules

```
B ? 0.5H to 0.7H
tstem,base ? H/16 to H/12
tbase ? B/8 to B/6
Ltoe ? B/3 to B/2
Lheel ? B/2 to 2B/3
```

### Bar Spacing Limits (CSA)

**Maximum spacing:**
```
smax = min(3h, 500 mm)  for walls
smax = min(2h, 300 mm)  for slabs in negative moment regions
```

**Minimum spacing:**
```
smin = max(db, 25 mm, 1.4?dagg)
```

### Standard Bar Spacings (mm)

```
100, 125, 150, 175, 200, 250, 300
```

Use these for easier construction and detailing.

---

## 7. CONSTRUCTION DETAILS

### Concrete Cover (CSA A23.1)

| Exposure | Cover (mm) |
|----------|------------|
| Earth face | 75 |
| Formed face (protected) | 50 |
| Formed face (exposed) | 60 |

### Hook and Bend Details

**Standard 90? Hook:**
```
Inside radius ? 3db
Extension ? 12db or 150 mm
```

**Standard 180? Hook:**
```
Inside radius ? 3db
Extension ? 4db or 65 mm
```

### Lap Splices

```
Llap = 1.0 ? ld  (tension)
Llap = 0.83 ? ld  (compression)

Minimum lap = 300 mm
```

---

## 8. SIGN CONVENTIONS

**Forces:**
- Horizontal forces: positive to right
- Vertical forces: positive downward
- Active pressure: positive (pushes wall)

**Moments:**
- Clockwise: positive
- Overturning: positive
- Resisting: positive

**Distances:**
- From toe: positive to heel
- From base: positive upward

---

## 9. UNIT CONVERSIONS

```
1 MPa = 1000 kPa = 1 N/mm?
1 kN/m? = 0.102 tf/m?
1 m = 1000 mm
1 kN?m = 1,000,000 N?mm
```

---

## 10. QUALITY CHECKS

Before finalizing design, verify:

? All safety factors met
? Reinforcement within limits (?min < ? < ?max)
? Spacing within code limits
? Development lengths adequate
? Crack control satisfied
? Bearing pressure < allowable
? Resultant within middle third
? Drainage provisions included
? Construction joints planned
? Seismic requirements met (if applicable)

---

## REFERENCES

1. CSA A23.3-19: Design of Concrete Structures
2. Foundation Engineering Handbook (2nd Ed.)
3. CRSI Design Handbook
4. Reinforced Concrete Mechanics and Design (Wight & MacGregor)

---

**Note:** All formulas implement SI units (kN, m, MPa, kPa) unless otherwise stated.

**Version:** 1.0 - October 2025
