import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')

# Exp A (Baseline)
file_A = 'EXP_AMP_0.1m.nc'  # Make sure this matches your local file
path_A = os.path.join(work_dir, file_A)

# Exp D (High Res)
file_D = 'EXP_D.nc' # Expected filename
path_D = os.path.join(work_dir, file_D)

def get_ssh_venice(filepath, label):
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None, None, None

    ds = netCDF4.Dataset(filepath)
    ssh = ds.variables['sossheig'][:]
    
    # Venice location (Northern end)
    # Generic approach: Find max latitude index
    nx = ssh.shape[2]
    ny = ssh.shape[1]
    
    # Grid coordinates if available, else indices
    # Venice is roughly at (i, j) = (nx-2, ny-2) for our idealized basin
    # For Exp A (GYRE=1), ny~22. For Exp D (GYRE=10), ny~220.
    # We select a point "at the north end"
    
    i_idx = nx - 2
    j_idx = ny - 5 # A bit off the wall to avoid boundary conditions
    
    print(f"[{label}] Extracting at index ({i_idx}, {j_idx}) from shape {ssh.shape}")
    
    timeseries = ssh[:, j_idx, i_idx]
    
    return timeseries, ds, (nx, ny)

def analyze_resolution():
    print("--- Experiment D: Resolution Analysis ---")
    
    # 1. Load Data
    ssh_A, ds_A, shape_A = get_ssh_venice(path_A, "Exp A")
    ssh_D, ds_D, shape_D = get_ssh_venice(path_D, "Exp D")
    
    if ssh_D is None:
        print("Exp D data missing. Please rename output file to 'EXP_D_highres.nc' and place in experiments/ folder.")
        return

    # Time Axes
    # Exp A: dt = 60s
    t_A = np.arange(len(ssh_A)) * 60.0 / 3600.0
    
    # Exp D: dt = 5s
    t_D = np.arange(len(ssh_D)) * 5.0 / 3600.0
    
    # 2. Compare Time Series at Venice
    plt.figure(figsize=(12, 6))
    
    if ssh_A is not None:
        plt.plot(t_A, ssh_A, 'k--', label=f'Exp A (dx~10km, dt=60s)')
        
    plt.plot(t_D, ssh_D, 'r-', linewidth=2, alpha=0.8, label=f'Exp D (dx~1km, dt=5s)')
    
    plt.title('Effect of Resolution on Kelvin Wave (Venice SSH)')
    plt.xlabel('Time (hours)')
    plt.ylabel('SSH (m)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(0, 24)
    
    out_compare = os.path.join(script_dir, 'fig_ExpD_resolution_comparison.png')
    plt.savefig(out_compare, dpi=300)
    print(f"Saved {out_compare}")
    
    # 3. Hovmöller for Exp D
    print("Generating Hovmöller for Exp D...")
    ny_D = shape_D[1]
    nx_D = shape_D[2]
    
    # Path along East coast
    idx_path = nx_D - 5
    
    # Full dataset extraction for Hov
    # SSH shape: (time, y, x)
    # Caution: High Res might be large in memory. Extract iteratively if needed.
    # 17280 steps x 1000 y is 17M points. Fine.
    
    ssh_full_D = ds_D.variables['sossheig']
    
    # Extract path slice: (time, y)
    hov_D = ssh_full_D[:, :, idx_path]
    
    plt.figure(figsize=(10, 8))
    vmax = np.percentile(np.abs(hov_D), 99)
    if vmax==0: vmax=0.05
    
    # dt=5s
    extent = [0, ny_D, 0, hov_D.shape[0]*5.0/3600.0] 
    
    plt.imshow(hov_D, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-vmax, vmax=vmax, extent=extent)
    
    plt.xlabel('Distance Index (South -> North)')
    plt.ylabel('Time (hours)')
    plt.title('Exp D: Hovmöller (High Res 1km)')
    plt.colorbar(label='SSH (m)')
    
    out_hov = os.path.join(script_dir, 'fig_ExpD_hovmoller.png')
    plt.savefig(out_hov, dpi=300)
    print(f"Saved {out_hov}")
    
    if ds_A: ds_A.close()
    if ds_D: ds_D.close()

if __name__ == "__main__":
    analyze_resolution()
