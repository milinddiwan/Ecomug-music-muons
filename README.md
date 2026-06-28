# muonsim — Cosmic Muon Simulation

Fortran 90 package that combines the **ECoMUG** surface-muon flux parametrization with the **MUSIC** Monte Carlo transport engine to simulate cosmic muons reaching underground or underwater depths.

Reference: D. Pagano et al., *NIM A* **1014** (2021) 165732

---

## Overview

The package has two main executables and a test/validation program:

| Executable | Source | Purpose |
|---|---|---|
| `muon_prop_200m` | `muon_prop_200m.f90` | Generate sea-level muons, propagate through 200 m of water, write survivors |
| `muon_range` | `muon_range.f90` | Compute range statistics for muons in rock at fixed initial energies |
| `test_ecomug` | `test_ecomug.f90` | Validate ECoMUG flux functions and Metropolis generators |

---

## Dependencies

- **gfortran** (GCC ≥ 9) or any Fortran 90 compiler
- Three MUSIC data files (included in repo):
  - `music-eloss-sr.dat`
  - `music-double-diff-rock.dat`
  - `music-cross-sections-sr.dat`
- Python 3 with `numpy` and `matplotlib` (optional, for plotting)

---

## Building

```bash
make          # builds muon_range, test_ecomug, muon_prop_200m
make clean    # remove objects, modules, and executables
```

The default compiler is `gfortran -O2`. Override with:

```bash
make FC=ifort FFLAGS="-O3 -xHost"
```

---

## Executables

### `muon_prop_200m`

Generates cosmic muons at sea level on a horizontal surface using the ECoMUG flux, propagates each through 200 m of water with MUSIC full 3-D Monte Carlo (bremsstrahlung, pair production, nuclear inelastic, ionisation, multiple scattering), and writes surviving muons to an ASCII file.

A pre-cut skips muons whose total energy is below the ionisation-only energy loss along the full slant path, before the expensive MUSIC call.

**Run:**

```bash
./muon_prop_200m
```

**Output — `muons_200m.dat`** (10 columns, one muon per row):

| Column | Symbol | Units | Description |
|---|---|---|---|
| 1 | `pmu0` | GeV/c | Initial momentum |
| 2 | `cth0` | — | Initial cos(zenith) |
| 3 | `phi0` | rad | Initial azimuthal angle |
| 4 | `emu_f` | GeV | Final total energy |
| 5 | `x_f` | cm | Final x position |
| 6 | `y_f` | cm | Final y position |
| 7 | `z_f` | cm | Final z (depth) |
| 8 | `cx_f` | — | Final direction cosine x |
| 9 | `cy_f` | — | Final direction cosine y |
| 10 | `cz_f` | — | Final direction cosine z (>0 = downward) |

**Key parameters** (edit `muon_prop_200m.f90` to change):

| Parameter | Default | Meaning |
|---|---|---|
| `N_TARGET` | 10000 | Stop after this many surviving muons |
| `NWARMUP` | 50000 | Metropolis burn-in steps |
| `DEPTH_M` | 200.0 m | Vertical depth of target surface |
| `CTH_MIN_PROP` | 0.02 | Drop near-horizontal tracks |
| `MEDIUM` | `'water'` | Propagation medium (`'water'`, `'seawater'`, `'rock'`) |

---

### `muon_range`

Propagates N_MU muons per energy bin through a medium until they stop (or exhaust `DEPTH_MAX`), reporting mean range, range spread, lateral displacement, and mean final angle.

**Run:**

```bash
./muon_range
```

**Output** (printed to stdout):

```
  E_mu     <R>        <R>      sigma_R  f_stop  <r_lat>  <theta_f>
  [GeV]    [cm]      [mwe]     [mwe]            [cm]     [deg]
```

- `<R>` — mean longitudinal range (depth of stopping point)
- `mwe` — meters water equivalent (g/cm²)
- `f_stop` — fraction of muons that stopped (should be 1.0; warn if not)
- `<r_lat>` — mean lateral displacement from beam axis
- `<theta_f>` — mean final zenith angle

**Key parameters** (edit `muon_range.f90` to change):

| Parameter | Default | Meaning |
|---|---|---|
| `N_MU` | 1000 | Muons per energy bin |
| `E_BINS` | 5–10000 GeV (8 bins) | Initial energies |
| `DEPTH_MAX` | 5×10⁶ cm | Maximum propagation depth |
| `MEDIUM` | `'rock'` | Medium (`'water'`, `'seawater'`, `'rock'`) |

---

### `test_ecomug`

Validates the ECoMUG flux functions and Metropolis generators with known reference values. Run after any modification to `ecomug.f90`.

```bash
./test_ecomug
```

---

## Source Modules

| File | Role |
|---|---|
| `ecomug.f90` | ECoMUG flux module — flux functions and Metropolis generators |
| `music-sr.f90` | MUSIC transport engine — `initialize_music`, `muon_transport` |
| `randn.f90` | Normal random number generator (Box-Muller) |
| `corset.f90` | Cholesky factorisation of a covariance matrix |
| `corgen.f90` | Correlated Gaussian sample generation |

### `ecomug_mod` public API

```fortran
use ecomug_mod
```

| Procedure | Signature | Description |
|---|---|---|
| `dalt_flux` | `(pmu, x) → real(8)` | Sea-level differential flux [m⁻² sr⁻¹ s⁻¹ (GeV/c)⁻¹] |
| `dmu_shallow` | `(pmu, x, depth_m) → real(8)` | Flux attenuated to vertical depth `depth_m` [m, water eq.] |
| `gen_horiz` | `(depth_m, nwarmup, n, pmu, cth, phi)` | Batch generator for a horizontal surface |
| `gen_vert` | `(depth_m, nwarmup, n, pmu, cth, phi_rel)` | Batch generator for a vertical surface |
| `gen_horiz_init` | `(depth_m, nwarmup, s_p, s_x, s_f)` | Warm up the horizontal Metropolis chain |
| `gen_horiz_step` | `(depth_m, s_p, s_x, s_f, pmu, cth, phi)` | Advance the chain by one thinned sample |

All generators use Metropolis-Hastings sampling with a multiplicative log₁₀(p) proposal (corrected Jacobian) and `NTHIN = 10` thinning steps per recorded sample. Momentum range: 1–2000 GeV/c.

### `music-sr.f90` interface

```fortran
call initialize_music(medium)          ! 'water', 'seawater', or 'rock'
call muon_transport(x, y, z, cx, cy, cz, emu, depth_max, ttime)
```

`muon_transport` propagates one muon in place. On return, `emu = 0` means the muon stopped; otherwise `emu` is the surviving total energy [GeV].

---

## Analysis and Plotting

### Python

```bash
python3 plot_muons.py               # reads muons_200m.dat
python3 plot_muons.py my_output.dat # read a different file
```

Produces `muon_distributions.svg` with four panels: initial momentum, final energy, initial zenith angle, and lateral displacement at depth.

Requires: `numpy`, `matplotlib`.

### Mathematica

- `read_muons.nb` — reads and analyses `muons_200m.dat`
- `ecomug.nb` — derivation and cross-checks of the ECoMUG parametrization

---

## File List

```
Makefile
ecomug.f90              ECoMUG flux module
ecomug.nb               Mathematica: ECoMUG derivation
music-sr.f90            MUSIC transport engine
music-eloss-sr.dat      MUSIC data: energy loss tables
music-double-diff-rock.dat
music-cross-sections-sr.dat
muon_prop_200m.f90      Main simulation (200 m water)
muon_range.f90          Range-vs-energy study
test_ecomug.f90         ECoMUG unit tests
randn.f90               Normal RNG
corset.f90              Cholesky factorisation
corgen.f90              Correlated Gaussian generator
plot_muons.py           Python: distribution plots
read_muons.nb           Mathematica: output analysis
```
