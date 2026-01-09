import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Configuration
input_file = 'experiments/EXP_AMP_0.1m.nc'
output_file = 'ssh_animation_expA.gif'
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

fig, ax = plt.subplots(figsize=(5, 8), constrained_layout=True)
fig.suptitle('SSH Evolution\nEXP A (Baseline)', fontsize=16)

vmax = 0.4
vmin = -vmax

# Use a mask for land/water if possible, or just plot everything
# Baseline mask (zeros at boundaries) might be noisy, but let's just plot directly.

im = ax.pcolormesh(nav_lon, nav_lat, ssh[0, :, :], cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_aspect('equal')
fig.colorbar(im, ax=ax, label='SSH (mm)')

time_text = ax.text(0.5, 1.02, '', transform=ax.transAxes, ha='center', fontsize=12)

def update(frame):
    im.set_array(ssh[frame, :, :].ravel())
    
    # 60s timestep
    hours = frame * 60 / 3600
    time_text.set_text(f"Time: {hours:.2f} h")
    
    if frame % 100 == 0:
        print(f"Processing frame {frame}/{time_steps}")
    return im, time_text

print("Generating animation...")
step = 10
end_step = 1020 # 17 hours
frames = list(range(0, min(end_step + 1, time_steps), step))

ani = animation.FuncAnimation(fig, update, frames=frames, blit=False, interval=50)

output_path = os.path.join(script_dir, output_file)
print(f"Saving to {output_path}...")
ani.save(output_path, writer='pillow', fps=fps)

print("Done!")
ds.close()
