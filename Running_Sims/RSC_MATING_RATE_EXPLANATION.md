# RSC to Mating Rate Translation Explanation

## Answer to PI's Question:

**RSC is NOT directly equal to the number of mates.** Instead, RSC is a scaling factor that affects the Poisson lambda parameter for mating rate.

### Current Implementation:

1. **Base mean mating rate**: 1.75 mates per female
2. **RSC scaling**: 
   - `rsc_factor = (female_rsc_phenotype / mean_pop_rsc)`
   - If female RSC equals population mean → factor = 1.0
   - Higher RSC → factor > 1.0
   - Lower RSC → factor < 1.0

3. **Lambda calculation**:
   - `lambda_mates = base_mean_mates * clamp(rsc_factor, 0.5, 2.0)`
   - Range: 1.75 × 0.5 = 0.875 to 1.75 × 2.0 = 3.5
   - Final lambda clamped to [0.8, 3.0]

4. **Number of mates**:
   - Drawn from Poisson(lambda)
   - Then clamped to [1, 5] mates

### What RSC = 1 Means:

If we interpret "RSC = 1" as meaning the RSC phenotype equals 1.0:
- This would give `rsc_factor ≈ 1.0 / mean_pop_rsc`
- If mean_pop_rsc ≈ 2.46 (from our simulation), then factor ≈ 0.41
- Lambda ≈ 1.75 × 0.41 ≈ 0.72
- Mean mates ≈ 0.72 (but clamped to minimum 1)

**However**, if "RSC = 1" means RSC equals the population mean (factor = 1.0):
- Lambda = 1.75 × 1.0 = 1.75
- Mean mates ≈ 1.75 (between 1 and 2)

### Conclusion:

**RSC = 1 does NOT mean monogamy (1 mate).** 

With the current implementation:
- RSC at population mean → ~1.75 mates on average (between 1 and 2)
- Lower RSC → fewer mates (closer to 1)
- Higher RSC → more mates (up to ~3-4 on average, max 5)

### Recommendation for PI:

To make RSC more interpretable, we could modify the code so:
- RSC = 1.0 → exactly 1 mate (monogamy)
- RSC = 2.0 → exactly 2 mates
- etc.

This would require changing the lambda calculation to be more directly proportional to RSC.

