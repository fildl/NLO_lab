import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Definition of configurations to process
configurations = [
    {'input': 'experiments/EXP_C_slope_100.nc', 'output': 'ssh_animation_expC.gif', 'title_suffix': ''},
    {'input': 'experiments/EXP_C_slope_100_1core.nc', 'output': 'ssh_animation_expC_1core.gif', 'title_suffix': ' (1-Core)'}
]

def generate_animation(config):
    input_file = config['input']
    output_file = config['output']
    title_suffix = config['title_suffix']
    fps = 15  # Frames per second

    # Get absolute paths
    file_path = os.path.join(script_dir, input_file)
    output_path = os.path.join(script_dir, output_file)

    if not os.path.exists(file_path):
        print(f"Skipping {input_file}: File not found at {file_path}")
        return

    print(f"Loading {file_path}...")
    ds = nc.Dataset(file_path)
    ssh = ds.variables['sossheig'][:] * 1000 # Time, Y, X (Convert to mm)
    nav_lon = ds.variables['nav_lon'][:]
    nav_lat = ds.variables['nav_lat'][:]

    time_steps = ssh.shape[0]
    print(f"Data shape: {ssh.shape}, Timesteps: {time_steps}")

    # For Exp C, dt = 20s. We need to check frame skipping.
    # Exp A had dt=60s and we skipped 10 frames (10 mins).
    # To have 10 mins per frame in Exp C (20s), we need step = 30.
    # 30 * 20s = 600s = 10 mins.

    fig, ax = plt.subplots(figsize=(5, 8), constrained_layout=True)
    fig.suptitle(f'SSH Evolution\nEXP C (Slope 1000m->100m)', fontsize=16)

    # Scale adjustment: Previous plot_exp_C_bathymetry showed Vmax around 0.54 mm.
    # But shoaling might produce higher transient peaks.
    # User used 0.8 mm for Exp A.
    # Let's use 0.8 mm for Exp C too for direct comparability, or slightly higher if needed.
    # Given Vmax ~ 0.54 mm from snapshots, 0.8 mm is safe and comparable.
    vmax = 0.4
    vmin = -vmax

    # Use a mask if available, but for simplicity/robustness just plot all.
    # Ideally we should mask the land like in snapshots, but pcolormesh handles it ok usually.

    im = ax.pcolormesh(nav_lon, nav_lat, ssh[0, :, :], cmap='RdBu_r', 
                        vmin=vmin, vmax=vmax, shading='auto')
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, label='SSH (mm)')

    time_text = ax.text(0.5, 1.02, '', transform=ax.transAxes, ha='center', fontsize=12)

    def update(frame):
        im.set_array(ssh[frame, :, :].ravel())
        
        # 20s timestep for Exp C
        hours = frame * 20 / 3600
        time_text.set_text(f"Time: {hours:.2f} h")
        
        if frame % 300 == 0:
            print(f"Processing frame {frame}/{time_steps}")
        return im, time_text

    print(f"Generating animation for {output_file}...")
    step = 30 # 10 minutes per frame (30 * 20s)
    # Limit to similar duration (17h like Exp A)
    # 17 * 3600 = 61200 s. 
    # Total steps for 17h: 61200 / 20 = 3060 steps.
    end_step = 3060
    frames = list(range(0, min(end_step + 1, time_steps), step))

    ani = animation.FuncAnimation(fig, update, frames=frames, blit=False, interval=50)

    print(f"Saving to {output_path}...")
    ani.save(output_path, writer='pillow', fps=fps)

    print("Done!")
    ds.close()
    plt.close(fig) # Close figure to free memory

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for config in configurations:
        generate_animation(config)
