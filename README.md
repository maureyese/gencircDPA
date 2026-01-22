# Genetic Circuits Simulation for DPAHelix Project

This repository corresponds to one of the five computational biology analyses of the DPAHelix project in the GOGEC Competition of 2026.

> **Objective**: Simulate the behavior of genetic circuits in engineered _B. subtilis_.

We developed a mathematical model of the SpoVF genetic circuit in engineered Bacillus subtilis to simulate the critical aspects from DNA transcription to DPA (dipicolinic acid) product synthesis. We analyzed transcriptional kinetics, translation dynamics, enzyme catalysis, and product accumulation across four constitutive promoters (p43, pylb, pyddf, pmsm) under varying metabolic conditions.

We worked with notebooks from Google Colab (https://colab.research.google.com/). This repository contains the finalized code and analysis.

## Key Resources

### Databases and Tools:
- **BioNumbers**: https://bionumbers.hms.harvard.edu/search.aspx
- **Salis Lab Promoter Calculator**: https://salislab.net/software/predict_promoter_calculator
- **Salis Lab RBS Calculator**: https://salislab.net/software/design_rbs_calculator
- **Benchling**: https://www.benchling.com/

### Python Libraries:
- `numpy`: Numerical computations
- `scipy.integrate`: ODE solvers (odeint, cumulative_trapezoid)
- `matplotlib`: Visualization

---

## Genetic Circuit Construction

Gene circuit sequences were stored in Benchling:

- pylb Promoter Gene Circuit: https://benchling.com/s/seq-yeFdvlVxDx5KuyEJT9Uk?m=slm-84pa474qitnLYuwiyLHR
- p43 Promoter Gene Circuit: https://benchling.com/s/seq-TvMQXjQG884q6chNU3Tc?m=slm-PkeQAjnd7UKXslDK0l3d
- pMSM Promoter Gene Circuit: https://benchling.com/s/seq-aRuC8I2qCFpeLJLDoCan?m=slm-3hCWN2WaaUbeBCutQDCu
- pYDDF Promoter Gene Circuit: https://benchling.com/s/seq-l8sWvjkjjitUSVCmJsSA?m=slm-3yy5lv6UMV8gGi3qGa4K

Specifications:

- Promoters were retrieved from NCBI: https://www.ncbi.nlm.nih.gov/
- RBS were designed according to CDS in Salis Lab Software: https://salislab.net/software/design_promoter_calculator
- CDS were retrieved from NCBI: https://www.ncbi.nlm.nih.gov/
- Terminator were retrieved from iGEM Parts: https://parts.igem.org/Main_Page

## Mathematical Modeling

We implemented a multi-stage model that couples transcription, translation, and enzymatic catalysis to predict DPA production in engineered B. subtilis.

### Activity 1: Promoter Transcription Rate Prediction

**Objective**: Determine the mRNA production rate for four constitutive promoters based on biological transcription kinetics.

#### Core Equation - mRNA Dynamics:
$$\frac{dm}{dt} = \underbrace{\alpha \cdot \text{Tx}_{\text{rate}} \cdot c}_{\text{Production Rate}} \;-\; \underbrace{d_1 \cdot m}_{\text{Degradation Rate}}$$

where:
- $m$: mRNA concentration (M)
- $\alpha$: Molarity per molecule per cell (M molecule⁻¹)
- $\text{Tx}_{\text{rate}}$: Transcription completion rate (transcripts s⁻¹)
- $c$: Gene copy number (dimensionless)
- $d_1$: mRNA degradation rate (s⁻¹)

#### Key Parameters:
- **Transcription initiation rate**: 20 transcripts min⁻¹ gene⁻¹
- **RNAP elongation rate**: 75 nt s⁻¹ in B. subtilis
- **Cell volume**: 0.9 µm³
- **mRNA half-life**: ~5 minutes
- **Gene length (SpoVF)**: 1437 nucleotides

#### Promoter Fold Changes (relative to p43):
| Promoter | Fold Change | Type |
|----------|-------------|------|
| p43 | 1.0× | Constitutive |
| pylb | 8.2× | Stationary-phase |
| pyddf | 2.4× | Log-phase |
| pmsm | 2.0× | Unknown |

#### Steady-State Solution:
$$m^* = \frac{k_{\text{on}} \cdot c}{d_1}$$

where $k_{\text{on}}$ is the molar transcription rate derived from:
$$k_{\text{on}} = \alpha \cdot \text{Tx}_{\text{rate}} \cdot c$$

---

### Activity 2: mRNA Concentration Prediction

**Objective**: Simulate temporal dynamics of mRNA accumulation for each promoter over 1 hour.

#### Analytical Solution:
$$m(t) = \frac{k_{\text{on}} \cdot c}{d_1} \left(1 - e^{-d_1 t}\right)$$

We solved this ODE using `scipy.integrate.odeint` for each promoter, capturing the approach to steady-state. Given the mRNA half-life of ~5 minutes, steady-state is reached within 20–30 minutes.

---

### Activity 3: Protein Concentration Prediction

**Objective**: Simulate coupled transcription-translation dynamics and protein accumulation over 24 hours.

#### Coupled ODE System - mRNA and Protein Dynamics:

**Equation (1) - mRNA dynamics:**
$$\frac{dm}{dt} = k_{\text{on}} \cdot c \;-\; d_1 \cdot m$$

**Equation (2) - Protein dynamics:**
$$\frac{dp}{dt} = \underbrace{k_{\text{trans}} \cdot m}_{\text{Translation Rate}} \;-\; \underbrace{d_2 \cdot p}_{\text{Protein Degradation}}$$

where:
- $p$: Protein concentration (M)
- $k_{\text{trans}}$: Translation rate constant (s⁻¹)
- $d_2$: Protein degradation rate (s⁻¹)

#### Key Parameters:
- **Translation initiation rate**: 5 initiations min⁻¹ mRNA⁻¹ → $k_{\text{trans}}$ = 0.083 s⁻¹
- **Ribosome elongation rate**: 33.3 amino acids s⁻¹
- **Protein half-life**: ~20 hours
- **RBS strength**: 7,500 AU (designed using Salis Lab RBS Calculator)

#### Steady-State Protein Concentration:
$$p^* = \frac{k_{\text{trans}} \cdot m^*}{d_2}$$

We solved the coupled system using `scipy.integrate.odeint` for all four promoters, capturing protein accumulation dynamics.

---

### Activity 4: Catalytic Velocity Prediction

**Objective**: Predict the enzymatic activity of SpoVF using Michaelis-Menten kinetics under different HTPA substrate concentrations.

#### Michaelis-Menten Equation:
$$v = \frac{V_{\text{max}} \cdot [S]}{K_m + [S]}$$

where:
- $v$: Reaction velocity (mM s⁻¹)
- $V_{\text{max}} = k_{\text{cat}} \cdot [E]$: Maximum reaction velocity
- $[S]$: Substrate concentration (HTPA, mM)
- $K_m$: Michaelis constant (mM)
- $[E]$: Enzyme concentration (SpoVF protein, mM)

#### Enzyme Kinetic Parameters:
- **$K_m$ (SpoVF)**: 0.776 mM
- **$k_{\text{cat}}$ (SpoVF)**: 0.301 s⁻¹

#### Substrate Conditions:
- **Glucose (Lower HTPA)**: ~0.0014 mM s⁻¹ (derived from 5 mmol gDW⁻¹ h⁻¹)
- **Glutamate (Higher HTPA)**: ~0.0019 mM s⁻¹ (derived from 6.7 mmol gDW⁻¹ h⁻¹)

We calculated instantaneous catalytic velocity over 4 hours by combining protein concentration profiles from Activity 3 with substrate availability.

---

### Activity 5: Product Accumulation Prediction (DPA Mass)

**Objective**: Translate enzymatic activity into practical bioprocess output by predicting DPA production over 24 hours under metabolic flux-based conditions.

#### Integration of Reaction Velocity:
$$\text{total\_mM\_internal} = \int_0^{T} v(t) \, dt$$

where $v(t)$ is calculated from Michaelis-Menten kinetics at each time point.

#### Scaling from Cellular Volume to Bioreactor Volume:
$$m_{\text{bioreactor}} \; [\text{mmol}] = m_{\text{internal}} \; [\text{mmol}] \times \left( C_{\text{biomass}} \; [\text{g L}^{-1}] \times \rho_{\text{cells}} \; [\text{cells g}^{-1}] \times V_{\text{cell}} \; [\text{L cell}^{-1}] \right)$$

#### Conversion to Mass:
$$C_{\text{dpa, g L}^{-1}} = \frac{m_{\text{bioreactor, mmol}} \cdot M_{\text{dpa, g mol}^{-1}}}{1000 \; \text{mmol mol}^{-1}}$$

#### Key Parameters:
- **Cell volume**: 0.9 × 10⁻¹⁵ L
- **Cell density**: 1.0 × 10¹² cells gDW⁻¹
- **DPA molecular weight**: 167.12 g mol⁻¹
- **Biomass concentrations tested**: 5, 15, 25 gDW L⁻¹
- **Production time**: 24 hours

We performed DPA production simulations using numerical integration (`scipy.integrate.cumulative_trapezoid`) and evaluated production across different biomass densities and substrate conditions.

---

## Outputs

The simulation generates four visualization outputs:

1. **01_mrna_concentrations.svg**: mRNA accumulation over 1 hour for all four promoters
2. **02_protein_concentrations.svg**: Protein accumulation over 24 hours for all four promoters
3. **03_catalytic_velocity.svg**: SpoVF catalytic activity under different HTPA conditions over 4 hours
4. **04_dpa_production_under_metabolic_flux.svg**: DPA production over 24 hours at 5 gDW/L biomass
5. **05_dpa_production_biomass_comparison_under_metabolic_flux.svg**: DPA production comparison at 5, 15, and 25 gDW/L biomass

---

## References

1. Yu et al. (2015). Promoter activity measurements. Nature Scientific Reports. https://doi.org/10.1038/srep18405
2. LaFleur et al. (2022). Salis Lab Promoter Calculator. https://doi.org/10.1038/s41467-022-32829-5
3. Salis et al. (2009). RBS Calculator for synthetic biology. Nature Biotechnology. https://www.nature.com/articles/nbt.1568
4. McClintock et al. (2025). SpoVF enzyme characterization and kinetic parameters.
5. BioNumbers: https://bionumbers.hms.harvard.edu