import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Configuration
input_file_1 = 'experiments/EXP_AMP_0.1m.nc'
input_file_2 = 'experiments/EXP_E.nc'
output_file = 'ssh_comparison.mp4'  # or .gif if mp4writers are not available
fps = 15  # Frames per second for the animation

# Get absolute paths
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path_1 = os.path.join(script_dir, input_file_1)
file_path_2 = os.path.join(script_dir, input_file_2)

print(f"Loading {file_path_1}...")
ds1 = nc.Dataset(file_path_1)
ssh1 = ds1.variables['sossheig'][:] * 100 # Time, Y, X (Convert to cm)
nav_lon1 = ds1.variables['nav_lon'][:]
nav_lat1 = ds1.variables['nav_lat'][:]

print(f"Loading {file_path_2}...")
ds2 = nc.Dataset(file_path_2)
ssh2 = ds2.variables['sossheig'][:] * 100 # Time, Y, X (Convert to cm)
nav_lon2 = ds2.variables['nav_lon'][:]
nav_lat2 = ds2.variables['nav_lat'][:]

# Ensure dimensions match (assuming same grid)
assert ssh1.shape == ssh2.shape, f"Shape mismatch: {ssh1.shape} vs {ssh2.shape}"
time_steps = ssh1.shape[0]

print(f"Data shape: {ssh1.shape}, Timesteps: {time_steps}")

# Setup Figure
fig, axes = plt.subplots(1, 2, figsize=(10, 10), constrained_layout=True)
fig.suptitle('SSH Evolution Comparison: EXP_AMP_0.1m vs EXP_E', fontsize=16)

# Plot settings
# Using same vmin/vmax for direct comparison
vmin = -0.5 # (1.0 = 1 cm)
vmax = 0.5

# Initial plots
# EXP_AMP_0.1m
im1 = axes[0].pcolormesh(nav_lon1, nav_lat1, ssh1[0, :, :], cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
axes[0].set_title('EXP_AMP_0.1m (Baseline)')
axes[0].set_xlabel("Longitude (°E)")
axes[0].set_ylabel("Latitude (°N)")
axes[0].set_aspect('equal')
fig.colorbar(im1, ax=axes[0], label='SSH (cm)')

# EXP_E
im2 = axes[1].pcolormesh(nav_lon2, nav_lat2, ssh2[0, :, :], cmap='RdBu_r', 
                       vmin=vmin, vmax=vmax, shading='auto')
axes[1].set_title('EXP_E (f=0)')
axes[1].set_xlabel("Longitude (°E)")
axes[1].set_yticklabels([]) # Hide Y labels for 2nd plot
axes[1].set_aspect('equal')
fig.colorbar(im2, ax=axes[1], label='SSH (cm)')

# Time text
time_text = fig.text(0.5, 0.02, '', ha='center', fontsize=12)

def update(frame):
    # Update data
    # pcolormesh returns a QuadMesh, we need to update the array
    # set_array expects a flattened array for QuadMesh
    im1.set_array(ssh1[frame, :, :].ravel())
    im2.set_array(ssh2[frame, :, :].ravel())
    
    # Update time text
    # Assuming 60s timestep
    hours = frame * 60 / 3600
    time_text.set_text(f"Time: {hours:.2f} h")
    
    if frame % 100 == 0:
        print(f"Processing frame {frame}/{time_steps}")
    
    return im1, im2, time_text

print("Generating animation...")
# Skip some frames if it's too heavy? Let's do step=2 or step=5 for speed if needed.
# For now, let's do step=10 (every 10 minutes) to make it reasonably fast to generate and watch
step = 10
# Limit to 17 hours
# 17 hours * 60 mins/hour = 1020 steps (since dt=60s)
end_step = 1020 
frames = list(range(0, end_step + 1, step))

ani = animation.FuncAnimation(fig, update, frames=frames, blit=False, interval=50, save_count=len(frames))

output_path = os.path.join(script_dir, output_file)
print(f"Saving to {output_path}...")

try:
    # Try saving as mp4
    ani.save(output_path, writer='ffmpeg', fps=fps)
except Exception as e:
    print(f"Error saving mp4: {e}")
    print("Falling back to GIF...")
    output_path = output_path.replace('.mp4', '.gif')
    ani.save(output_path, writer='pillow', fps=fps)

print("Done!")
ds1.close()
ds2.close()
