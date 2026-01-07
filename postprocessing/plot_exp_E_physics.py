import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import os

# Configuration
input_file = 'experiments/EXP_E.nc'
output_file = 'fig_ExpE_snapshots.png'
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, input_file)

print(f"Loading {file_path}...")
ds = nc.Dataset(file_path)
ssh = ds.variables['sossheig'][:]  # Time, Y, X
nav_lon = ds.variables['nav_lon'][:]
nav_lat = ds.variables['nav_lat'][:]
time_steps = ssh.shape[0]

print(f"Data shape: {ssh.shape}, Timesteps: {time_steps}")

# Same timestamps as Exp A (2.0h, 5.5h, 9.0h, 12.5h, 16.0h)
# 60s timestep => 120, 330, 540, 750, 960 indices
indices = np.linspace(120, 960, 5, dtype=int)
times_hours = [i * 60 / 3600 for i in indices]

fig, axes = plt.subplots(1, 5, figsize=(15, 10), constrained_layout=True)
fig.suptitle('Experiment E: Physics Validation (f=0)', fontsize=20)

# Match Exp A scale (approx +/- 0.008m based on user input)
vmin = -0.008
vmax = 0.008

for ax, idx, th in zip(axes, indices, times_hours):
    data = ssh[idx, :, :]
    # Mask zero values if needed, or just plot
    # Assuming land is 0 or masked. In NEMO output, usually 0.
    
    im = ax.pcolormesh(nav_lon, nav_lat, data, cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
    
    ax.set_title(f"T = {th} h")
    ax.set_xlabel("Longitude (°E)")
    if ax == axes[0]:
        ax.set_ylabel("Latitude (°N)")
    else:
        ax.set_yticklabels([])
    
    ax.set_aspect('equal')

# Colorbar
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.05)
cbar.set_label('SSH Anomaly (m)')

output_path = os.path.join(script_dir, output_file)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Saved {output_path}")


# ds is already closed in original script logic, but let's keep it open or reopen
# Actually, the original script does not close it until the end.
# Let's add Hovmoller logic before closing.

# --- Hovmöller Diagram ---
# Path same as Exp A (East Coast, approx index 100)
# Data shape (1440, 102, 22) => (Time, Lat, Lon)
# East coast is at high X index. Let's pick index 20 (it is narrow).
# Wait, let's check nav_lon from previous run.
# Since domain is narrow (Exp A used nx-2).
nx = ssh.shape[2]
ny = ssh.shape[1]
east_idx = nx - 2

hovmoller = np.zeros((time_steps, ny))
for t in range(time_steps):
    for j in range(ny):
        hovmoller[t, j] = ssh[t, j, east_idx]

# Plot
fig2, ax2 = plt.subplots(figsize=(8, 6), constrained_layout=True)
# Time axis in hours
y_time = np.arange(time_steps) * 60 / 3600
# Space axis in indices (approx km)
x_space = np.arange(ny)

im2 = ax2.imshow(hovmoller, aspect='auto', origin='lower', cmap='RdBu_r',
                 extent=[0, ny, 0, y_time[-1]], vmin=vmin, vmax=vmax)
ax2.set_xlabel('Latitude Index (South -> North)')
ax2.set_ylabel('Time (hours)')
ax2.set_title('Hovmöller Diagram (East Coast) - Exp E (f=0)')
fig2.colorbar(im2, ax=ax2, label='SSH Anomaly (m)')

out_hov = os.path.join(script_dir, 'fig_ExpE_hovmoller.png')
plt.savefig(out_hov, dpi=300, bbox_inches='tight')
print(f"Saved {out_hov}")

ds.close()
