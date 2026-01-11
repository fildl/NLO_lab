import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
# Data is in 'experiments'
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_C_slope_100_1core.nc'
filepath = os.path.join(work_dir, filename)

# Mesh mask candidates for Exp C 1core
mesh_candidates = [
    os.path.join(script_dir, 'mesh', 'mesh_mask_C_1core.nc'),
    os.path.join(work_dir, 'mesh_mask_C_1core.nc'),
    os.path.join(script_dir, 'experiments', 'mesh_mask.nc')
]

def analyze_initial_state():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    print(f"Analyzing {filename} for Initial State (T=0)...")
    ds = netCDF4.Dataset(filepath)
    
    # Load SSH variable
    ssh_var = ds.variables['sossheig']
    # Initial state at t=0
    ssh_t0 = ssh_var[0, :, :] 
    
    # Try to load mask
    tmask = None
    meshpath = ""
    for p in mesh_candidates:
        if os.path.exists(p):
            meshpath = p
            break
            
    if meshpath:
        try:
            with netCDF4.Dataset(meshpath) as mds:
                if 'tmask' in mds.variables:
                    temp_mask = mds.variables['tmask'][0,0,:,:]
                    if temp_mask.shape == ssh_t0.shape:
                        tmask = temp_mask
                        print(f"Loaded tmask from {os.path.basename(meshpath)}")
                    else:
                        print(f"Warning: Mask dimensions mismatch. Data: {ssh_t0.shape}, Mask: {temp_mask.shape}")
        except Exception as e:
            print(f"Warning: Failed to load mesh mask: {e}")
            
    if tmask is None:
        print("Creating synthetic mask.")
        tmask = np.ones_like(ssh_t0)
        tmask[0,:] = 0; tmask[-1,:] = 0
        tmask[:,0] = 0; tmask[:,-1] = 0
    
    # Load Georeferenced Coordinates
    if 'nav_lon' in ds.variables:
        lon = ds.variables['nav_lon'][:]
        lat = ds.variables['nav_lat'][:]
    elif meshpath:
        try:
            with netCDF4.Dataset(meshpath) as mds:
                if 'glamt' in mds.variables:
                    lon = mds.variables['glamt'][0,:,:]
                    lat = mds.variables['gphit'][0,:,:]
                else:
                    lon, lat = None, None
        except:
            lon, lat = None, None
    else:
        lon, lat = None, None
        
    ds.close()

    # Convert to mm
    ssh_t0_mm = ssh_t0 * 1000.0
    vmax_init = 0.4 # Consistent with Exp A plots
    
    # Plotting
    fig_init, ax_init = plt.subplots(figsize=(6, 8), constrained_layout=True)
    
    if lon is not None and lat is not None:
        # Apply mask for pcolormesh if possible, or just plot regular
        # If using pcolormesh with mask, we usually mask the data
        ssh_plot = np.ma.masked_where(tmask == 0, ssh_t0_mm)
        
        im_init = ax_init.pcolormesh(lon, lat, ssh_plot, cmap='RdBu_r', shading='auto', vmin=-vmax_init, vmax=vmax_init)
        ax_init.set_xlabel('Longitude (°E)')
        ax_init.set_ylabel('Latitude (°N)')
        ax_init.set_aspect('equal')
    else:
        ssh_plot = np.ma.masked_where(tmask == 0, ssh_t0_mm)
        im_init = ax_init.imshow(ssh_plot, origin='lower', cmap='RdBu_r', vmin=-vmax_init, vmax=vmax_init)
        ax_init.set_xlabel('I')
        ax_init.set_ylabel('J')

    ax_init.set_title(f'Initial State (T=0)\nSSH Anomaly (mm)', fontsize=14)
    fig_init.colorbar(im_init, ax=ax_init, label='SSH (mm)', shrink=0.8)
    
    output_filename = 'fig_expC_1core_initial_state.png'
    output_path_init = os.path.join(script_dir, output_filename)
    plt.savefig(output_path_init, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path_init}")

if __name__ == "__main__":
    analyze_initial_state()
