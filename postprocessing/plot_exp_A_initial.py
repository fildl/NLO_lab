import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_AMP_0.1m.nc'
filepath = os.path.join(work_dir, filename)

# Try to find mesh mask
meshpath = ""
meshpath_candidates = [
    os.path.join(script_dir, 'mesh', 'mesh_mask_A.nc'),
    os.path.join(work_dir, 'mesh_mask_A.nc'),
]
for p in meshpath_candidates:
    if os.path.exists(p):
        meshpath = p
        break

def plot_initial_state():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    ds = netCDF4.Dataset(filepath)
    ssh = ds.variables['sossheig'][:]

    # Get coordinates
    if 'nav_lon' in ds.variables:
        lon = ds.variables['nav_lon'][:]
        lat = ds.variables['nav_lat'][:]
    elif os.path.exists(meshpath):
        with netCDF4.Dataset(meshpath) as mds:
            lon = mds.variables['glamt'][0,:,:]
            lat = mds.variables['gphit'][0,:,:]
    else:
        print("Warning: No coordinates found. Using indices.")
        lat, lon = None, None

    # T=0 Snapshot
    t_step = 0
    ssh_t0 = ssh[t_step, :, :] * 1000 # Convert to mm

    # Plot settings
    plt.rcParams.update({'font.size': 12})
    fig, ax = plt.subplots(figsize=(6, 8), constrained_layout=True)

    vmax = 0.4
    
    if lon is not None:
        im = ax.pcolormesh(lon, lat, ssh_t0, cmap='RdBu_r', shading='auto', vmin=-vmax, vmax=vmax)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_aspect('equal')
    else:
        im = ax.imshow(ssh_t0, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        ax.set_xlabel('I')
        ax.set_ylabel('J')

    ax.set_title(f'Initial State (T=0)\nSSH Anomaly (mm)', fontsize=14)
    plt.colorbar(im, ax=ax, label='SSH (mm)', shrink=0.8)
    
    output_path = os.path.join(script_dir, 'fig_expA_initial_state.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path}")

    ds.close()

if __name__ == "__main__":
    plot_initial_state()
