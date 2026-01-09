import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Configuration
input_file = 'experiments/EXP_C_slope_100_sigma.nc'
output_file = 'ssh_animation_expC_sigma.gif'
fps = 15  # Frames per second

# Get absolute paths
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, input_file)

if not os.path.exists(file_path):
    print(f"Error: Data file not found at {file_path}")
    exit(1)

print(f"Loading {file_path}...")
ds = nc.Dataset(file_path)
ssh = ds.variables['sossheig'][:] * 1000 # Time, Y, X (Convert to mm)
nav_lon = ds.variables['nav_lon'][:]
nav_lat = ds.variables['nav_lat'][:]

time_steps = ssh.shape[0]
print(f"Data shape: {ssh.shape}, Timesteps: {time_steps}")

# For Exp C Sigma, dt = 20s. Same as Exp C.
# 30 * 20s = 600s = 10 mins.

fig, ax = plt.subplots(figsize=(5, 8), constrained_layout=True)
fig.suptitle('SSH Evolution\nEXP C Sigma (Slope 1000m->100m)', fontsize=16)

# Scale adjustment for comparability
vmax = 0.4
vmin = -vmax

# Ideally mask land, but for quick animation pcolormesh handles defaults well
if np.ma.is_masked(ssh):
  ssh_filled = ssh.filled(0)
else:
  ssh_filled = ssh

im = ax.pcolormesh(nav_lon, nav_lat, ssh[0, :, :], cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_aspect('equal')
fig.colorbar(im, ax=ax, label='SSH (mm)')

# Add 'Sigma' annotation roughly
ax.text(0.05, 0.05, 'Sigma Coords', transform=ax.transAxes, color='black', fontsize=10, 
        bbox=dict(facecolor='white', alpha=0.7))

time_text = ax.text(0.5, 1.02, '', transform=ax.transAxes, ha='center', fontsize=12)

def update(frame):
    im.set_array(ssh[frame, :, :].ravel())
    
    # 20s timestep
    hours = frame * 20 / 3600
    time_text.set_text(f"Time: {hours:.2f} h")
    
    if frame % 300 == 0:
        print(f"Processing frame {frame}/{time_steps}")
    return im, time_text

print("Generating animation...")
step = 30 # 10 minutes per frame
# Limit to similar duration (~17h)
end_step = 3060
frames = list(range(0, min(end_step + 1, time_steps), step))

ani = animation.FuncAnimation(fig, update, frames=frames, blit=False, interval=50)

output_path = os.path.join(script_dir, output_file)
print(f"Saving to {output_path}...")
ani.save(output_path, writer='pillow', fps=fps)

print("Done!")
ds.close()
