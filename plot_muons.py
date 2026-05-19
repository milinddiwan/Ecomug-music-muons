#!/usr/bin/env python3
"""
plot_muons.py — distributions of muons surviving 200 m of water.

Reads muons_200m.dat produced by muon_prop_200m and plots:
  - Initial momentum distribution
  - Final energy distribution
  - Initial zenith angle distribution
  - Lateral displacement at 200 m depth

Usage:
    python3 plot_muons.py [muons_200m.dat]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('svg')   # vector backend; works without a display
import matplotlib.pyplot as plt

DATA_FILE = sys.argv[1] if len(sys.argv) > 1 else 'muons_200m.dat'

# ---- Load data -------------------------------------------------------
try:
    data = np.loadtxt(DATA_FILE, comments='#')
except FileNotFoundError:
    sys.exit(f'Error: {DATA_FILE} not found — run muon_prop_200m first.')

if data.ndim == 1:
    data = data.reshape(1, -1)

pmu0  = data[:, 0]   # initial momentum [GeV/c]
cth0  = data[:, 1]   # initial cos(zenith)
phi0  = data[:, 2]   # initial azimuthal [rad]
emu_f = data[:, 3]   # final total energy [GeV]
x_f   = data[:, 4]   # final x [cm]
y_f   = data[:, 5]   # final y [cm]
cx_f  = data[:, 7]
cy_f  = data[:, 8]
cz_f  = data[:, 9]

th0_deg = np.degrees(np.arccos(np.clip(cth0, -1, 1)))   # initial zenith [deg]
thf_deg = np.degrees(np.arccos(np.clip(cz_f, -1, 1)))   # final zenith [deg]
r_lat   = np.sqrt(x_f**2 + y_f**2) * 1e-2               # lateral displacement [m]

n = len(pmu0)
print(f'Loaded {n:,} muons from {DATA_FILE}')
print(f'  pmu0  : {pmu0.min():.1f} – {pmu0.max():.1f} GeV/c  (mean {pmu0.mean():.1f})')
print(f'  emu_f : {emu_f.min():.2f} – {emu_f.max():.2f} GeV   (mean {emu_f.mean():.2f})')
print(f'  th0   : {th0_deg.min():.1f} – {th0_deg.max():.1f} deg   (mean {th0_deg.mean():.1f})')
print(f'  r_lat : {r_lat.min():.1f} – {r_lat.max():.1f} m     (mean {r_lat.mean():.1f})')

# ---- Figure ----------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle(
    f'Cosmic muons surviving 200 m of water  (N = {n:,})',
    fontsize=13, fontweight='bold'
)

kw = dict(edgecolor='white', linewidth=0.4)

# 1. Initial momentum
ax = axes[0, 0]
bins = np.linspace(pmu0.min() * 0.95, pmu0.max() * 1.02, 40)
ax.hist(pmu0, bins=bins, color='steelblue', **kw)
ax.set_xlabel('Initial momentum $p_0$  (GeV/$c$)', fontsize=11)
ax.set_ylabel('Counts / bin', fontsize=10)
ax.set_title('Initial momentum')
ax.axvline(pmu0.mean(), color='navy', linestyle='--', lw=1.2,
           label=f'mean = {pmu0.mean():.1f} GeV/$c$')
ax.legend(fontsize=9)

# 2. Final energy
ax = axes[0, 1]
bins = np.linspace(0, emu_f.max() * 1.05, 50)
ax.hist(emu_f, bins=bins, color='tomato', **kw)
ax.set_xlabel('Final energy $E_f$  (GeV)', fontsize=11)
ax.set_ylabel('Counts / bin', fontsize=10)
ax.set_title('Final energy at 200 m depth')
ax.axvline(emu_f.mean(), color='darkred', linestyle='--', lw=1.2,
           label=f'mean = {emu_f.mean():.2f} GeV')
ax.legend(fontsize=9)

# 3. Initial zenith angle
ax = axes[1, 0]
bins = np.arange(0, th0_deg.max() + 3, 2)
ax.hist(th0_deg, bins=bins, color='mediumseagreen', **kw)
ax.set_xlabel(r'Initial zenith angle $\theta_0$  (deg)', fontsize=11)
ax.set_ylabel('Counts / bin', fontsize=10)
ax.set_title('Initial zenith angle')
ax.axvline(th0_deg.mean(), color='darkgreen', linestyle='--', lw=1.2,
           label=f'mean = {th0_deg.mean():.1f}°')
ax.legend(fontsize=9)

# 4. Lateral displacement
ax = axes[1, 1]
p99 = np.percentile(r_lat, 99)
bins = np.linspace(0, p99 * 1.05, 50)
ax.hist(r_lat, bins=bins, color='mediumpurple', **kw)
ax.set_xlabel('Lateral displacement $r = \\sqrt{x^2+y^2}$  (m)', fontsize=11)
ax.set_ylabel('Counts / bin', fontsize=10)
ax.set_title('Lateral displacement at 200 m depth')
ax.axvline(r_lat.mean(), color='indigo', linestyle='--', lw=1.2,
           label=f'mean = {r_lat.mean():.1f} m')
ax.legend(fontsize=9)

for ax in axes.flat:
    ax.grid(True, alpha=0.25, linestyle='--')
    ax.tick_params(labelsize=9)

plt.tight_layout()
outfile = 'muon_distributions.svg'
plt.savefig(outfile, bbox_inches='tight')
print(f'Saved {outfile}')
plt.show()
