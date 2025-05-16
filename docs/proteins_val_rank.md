# **Scoring Criteria for Metric Calculation**

| **Metric**            | **Weight** | **Criteria**                                              | **Range**            | **Better When**            |
|-----------------------|-----------|-----------------------------------------------------------|----------------------|----------------------------|
| **Binding Affinity**  | **0.30**   | Dissociation Constant (Kd)                                 | e-3 to e-12           | Closer to **e-12**          |
|                       |           |                                                         | > e-12 → **0**      |                            |
|                       |           |                                                         | < e-3 → **0**       |                            |
| **GMQE Score**        | **0.20**   | Score                                                   | (0.6 - 1)           | Closer to **1**            |
| **Glycosylation Sites** | **0.10**  | N-Glycosylation Sites Count                            | (1 - 4)             | Closer to **4**            |
| **Aggregation**       | **0.10**   | Aggregation Propensity (Low, Medium, High)            | Low                 | **Low is better**          |
|                       |           | Aggregation-Prone Regions                             | (0 - 8)             | Closer to **0**            |
| **ProtParam**         | **0.10**   | GRAVY Score                                           | (-0.5, -1.5)        | Closer to **-0.5**         |
|                       |           | Solubility                                            | Soluble/Not Soluble | **Soluble is better**      |
| **Immunogenicity**    | **0.05**   | Immunogenic Score                                     | (0 - 1)             | Closer to **1**            |
| **Conservancy**       | **0.05**   | Conservancy Score                                     | (0 - 1)             | Closer to **1**            |
| **Stability**         | **0.05**   | Melting Temperature (Tm)                              | ≤ 65°C → **0** & > 90°C → **0**      |                            |
|                       |           |                                                       | 65°C to 90°C        | Closer to **90°C**         |
| **Epitope**          | **0.03**   | Epitope Score                                        | (0 - 1)             | Closer to **1**            |
| **Developability**    | **0.02**   | Developability Score                                 | (0 - 1)             | Closer to **1**            |

## 1. Binding Affinity (Weight: 0.30)

- Dissociation Constant (Kd): Threshold between e-3 and e-12
  - Closer to e-12 is better → Higher score.
  - Anything above e-12 → Score = 0.
  - Anything below e-3 → Score = 0.

## 2. GMQE Score (Weight: 0.20)

- Score Range: 0.6 - 1 → Higher is better (Closer to 1 gets a higher score).

## 3. Glycosylation Sites (Weight: 0.10)

- N-Glycosylation Sites Count: Range (1, 4) → Closer to 4 is better (Higher score).

## 4. Aggregation (Weight: 0.10)

- Aggregation Propensity: Categories → Low, Medium, High
  - Low is better → Higher score.
- Aggregation-Prone Regions: Range (0, 8) → Closer to 0 is better (Higher score).

## 5. ProtParam (Weight: 0.10)

- GRAVY Score: Range (-0.5, -1.5) → Closer to -0.5 is better (Higher score).
- Solubility: Soluble is better → Higher score.

## 6. Immunogenicity (Weight: 0.05)

- Immunogenic Score: Range (0-1) → Higher is better (Closer to 1 gets a higher score).

## 7. Conservancy (Weight: 0.05)

- Conservancy Score: Range (0-1) → Higher is better (Closer to 1 gets a higher score).

## 8. Stability (Weight: 0.05)

- Melting Temperature (Tm):
  - Tm ≤ 65°C → Score = 0.
  - 65°C to 90°C → Higher is better (Closer to 90°C gets a higher score).

## 9. Epitope (Weight: 0.03)

- Epitope Score: Range (0-1) → Higher is better (Closer to 1 gets a higher score).

## 10. Developability (Weight: 0.02)

- Developability Score: Range (0-1) → Higher is better (Closer to 1 gets a higher score)
