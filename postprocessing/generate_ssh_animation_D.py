import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Configuration
input_file = 'experiments/EXP_D.nc'
output_file = 'ssh_animation_expD.gif'
fps = 15  # Frames per second

# Get absolute paths
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, input_file)

print(f"Loading {file_path}...")
ds = nc.Dataset(file_path)
ssh = ds.variables['sossheig'][:] * 1000 # Time, Y, X (Convert to mm)
nav_lon = ds.variables['nav_lon'][:]
nav_lat = ds.variables['nav_lat'][:]

time_steps = ssh.shape[0]
print(f"Data shape: {ssh.shape}, Timesteps: {time_steps}")

# Exp D setup: dt = 20s (4320 steps for 24h)
# Exp A had dt=60s and we skipped 10 frames (10 mins).
# To have 10 mins per frame in Exp D (20s), we need step = 30.
# 30 * 20s = 600s = 10 mins.

fig, ax = plt.subplots(figsize=(6, 8), constrained_layout=True)
fig.suptitle('SSH Evolution\nEXP D (High Res)', fontsize=16)

# Scale: Exp D is flat bottom like Exp A.
# Expected amplitudes should be similar to Exp A.
# Exp A Vmax was 0.8 mm. Let's use 0.8 mm here too for comparison.
vmax = 0.4
vmin = -vmax

# Use a mask if available, but for simplicity/robustness just plot all.

im = ax.pcolormesh(nav_lon, nav_lat, ssh[0, :, :], cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_aspect('equal')
fig.colorbar(im, ax=ax, label='SSH (mm)')

time_text = ax.text(0.5, 1.02, '', transform=ax.transAxes, ha='center', fontsize=12)

def update(frame):
    im.set_array(ssh[frame, :, :].ravel())
    
    # 20s timestep for Exp D
    hours = frame * 20 / 3600
    time_text.set_text(f"Time: {hours:.2f} h")
    
    if frame % 300 == 0:
        print(f"Processing frame {frame}/{time_steps}")
    return im, time_text

print("Generating animation...")
step = 30 # 10 minutes per frame (30 * 20s)
# Limit to similar duration (17h like Exp A)
# 17 * 3600 = 61200 s. 
# Total steps for 17h: 61200 / 20 = 3060 steps.
end_step = 3060
frames = list(range(0, min(end_step + 1, time_steps), step))

ani = animation.FuncAnimation(fig, update, frames=frames, blit=False, interval=50)

output_path = os.path.join(script_dir, output_file)
print(f"Saving to {output_path}...")
ani.save(output_path, writer='pillow', fps=fps)

print("Done!")
ds.close()
