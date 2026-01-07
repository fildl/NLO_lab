import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
# Data is now in 'experiments' directly
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_AMP_0.1m.nc' 
filepath = os.path.join(work_dir, filename)
# Mesh mask is expected in postprocessing/mesh/ or experiments/
# User said they have mesh_mask_C.nc. Let's try to look for it.
meshpath_candidates = [
    os.path.join(script_dir, 'mesh', 'mesh_mask_A.nc'),
    os.path.join(work_dir, 'mesh_mask_A.nc'),
]
meshpath = ""
for p in meshpath_candidates:
    if os.path.exists(p):
        meshpath = p
        break

def analyze_wave():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    ds = netCDF4.Dataset(filepath)
    
    ssh = ds.variables['sossheig'][:]
    
    # Try to load mask
    tmask = None
    if os.path.exists(meshpath):
        try:
            with netCDF4.Dataset(meshpath) as mds:
                temp_mask = mds.variables['tmask'][0,0,:,:]
            # Check dimensions
            if temp_mask.shape == ssh[0,:,:].shape:
                tmask = temp_mask
                print("Loaded tmask from mesh file.")
            else:
                print(f"Warning: Mesh mask dimensions {temp_mask.shape} do not match data {ssh[0,:,:].shape}. Ignoring mesh mask.")
        except Exception as e:
            print(f"Warning: Failed to load mesh mask: {e}")
            
    if tmask is None:
        print("Creating synthetic mask (assuming box domain).")
        tmask = np.ones_like(ssh[0,:,:]) 
        # Set boundaries to land (0)
        tmask[0,:] = 0; tmask[-1,:] = 0
        tmask[:,0] = 0; tmask[:,-1] = 0
    
    # Load Georeferenced Coordinates
    if 'nav_lon' in ds.variables:
        lon = ds.variables['nav_lon'][:]
        lat = ds.variables['nav_lat'][:]
    elif os.path.exists(meshpath):
        with netCDF4.Dataset(meshpath) as mds:
            lon = mds.variables['glamt'][0,:,:]
            lat = mds.variables['gphit'][0,:,:]
    else:
        # Fallback to simple indices if absolutely nothing available
        print("Warning: No coordinates found. Using indices.")
        lat, lon = np.meshgrid(np.arange(ssh.shape[1]), np.arange(ssh.shape[2]), indexing='ij')

    # --- 1. Hovmöller Diagram (Index-based) ---
    print("Generating Hovmöller Diagram...")
    ny, nx = ssh.shape[1], ssh.shape[2]
    
    # Robust East Coast Path: Fixed index due to idealized domain
    # Use fixed index near eastern boundary (nx-2) to avoid artifacts
    east_coast_idx = [nx - 2] * ny

    hovmoller = np.zeros((ssh.shape[0], ny))
    for t in range(ssh.shape[0]):
        for j in range(ny):
            idx = east_coast_idx[j]
            hovmoller[t, j] = ssh[t, j, int(idx)]
            
    # Plot Domain & Path (Requested "Similar to Exp C")
    if 'lon' in locals() and lon is not None:
        fig_path, ax_path = plt.subplots(figsize=(6, 8))
        # Simple domain mask
        if tmask is not None:
            ax_path.pcolormesh(lon, lat, tmask, cmap='Greys', alpha=0.3)
        
        # Path
        path_lons = [lon[j, nx-2] for j in range(ny)]
        path_lats = [lat[j, nx-2] for j in range(ny)]
        
        ax_path.plot(path_lons, path_lats, 'r-', linewidth=2, label='Hovmöller Path')
        ax_path.legend(loc='upper right')
        ax_path.set_title('Exp A: Extraction Path')
        ax_path.set_xlabel('Longitude (°E)')
        ax_path.set_ylabel('Latitude (°N)')
        ax_path.set_aspect('equal')
        plt.savefig(os.path.join(script_dir, 'fig_expA_path.png'), dpi=300, bbox_inches='tight')
        print("Saved fig_expA_path.png")

    plt.figure(figsize=(10, 8))
    # Use robust vmin/vmax for visibility
    h_vmax = np.percentile(np.abs(hovmoller), 99)
    if h_vmax == 0: h_vmax = 0.05
    
    plt.imshow(hovmoller, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, ny*10, 0, ssh.shape[0]*60/3600])
    plt.colorbar(label='SSH (m)')
    plt.xlabel('Distance along Coast (km approx)')
    plt.ylabel('Time (hours)')
    plt.title('Hovmöller Diagram (East Coast Path)')
    plt.savefig(os.path.join(script_dir, 'fig_expA_hovmoller_east.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expA_hovmoller_east.png')}")

    # --- 2. Snapshots (Index-based & Georeferenced) ---
    print("Generating Snapshots...")
    # Shifted times to avoid initial perturbation: 4h, 10h, 16h, 22h
    times_idx = [240, 600, 960, 1320]
    
    # Dynamic limit (exclude first 2 hours for scaling)
    vmax = np.max(np.abs(ssh[120:,:,:])) * 0.8
    if vmax < 1e-4: vmax = 0.05

    # Index-based Plot
    # Tall and narrow domain: Use square figure with constrained layout to center plots
    fig, axes = plt.subplots(1, 4, figsize=(10, 10), constrained_layout=True)
    
    for i, tidx in enumerate(times_idx):
        if tidx < ssh.shape[0]:
            ax = axes[i]
            data = ssh[tidx, :, :]
            data = np.ma.masked_where(tmask == 0, data)
            im = ax.imshow(data, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
            ax.set_title(f'T = {tidx*60/3600:.1f} h')
            
            # Horizontal colorbar for better use of space
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, orientation='horizontal')
            
            ax.set_xlabel('I')
            if i==0: ax.set_ylabel('J')
            else: ax.set_yticklabels([])
            
    plt.savefig(os.path.join(script_dir, 'fig_expA_ssh_snapshots_index.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expA_ssh_snapshots_index.png')}")

    # Georeferenced Plot
    print("Generating Georeferenced Snapshots...")
    # Plot 5 snapshots from 2h to 16h as requested
    # 2h = 120 steps, 16h = 960 steps
    steps = np.linspace(120, 960, 5, dtype=int)
    
    # Calculate robust Vmax from the selected frames ONLY
    # This ensures colorbar is not saturated by unshown transients
    selected_ssh = ssh[steps, :, :]
    vmax_local = np.percentile(np.abs(selected_ssh), 99.9) # Use 99.9 to avoid single pixel spikes
    print(f"Index-based Vmax: {vmax_local}")
    
    # Larger fonts
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14})
    
    fig, axes = plt.subplots(1, 5, figsize=(15, 10), constrained_layout=True)
    
    # Calculate bounds (re-check if lon is defined)
    if 'lon' not in locals() or lon is None: extent = None
    else:
        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()
        extent = [lon_min, lon_max, lat_min, lat_max]

    for i, t_step in enumerate(steps):
        if t_step >= ssh.shape[0]: continue
        ax = axes[i]
        data_step = ssh[t_step, :, :]
        
        if 'lon' in locals() and lon is not None:
            # Symmetic Vmin/Vmax
            im = ax.pcolormesh(lon, lat, data_step, cmap='RdBu_r', shading='auto', vmin=-vmax_local, vmax=vmax_local)
            ax.set_xlabel('Lon (°E)')
            if i == 0: ax.set_ylabel('Lat (°N)')
            if extent:
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
            # Adjust aspect ratio
            ax.set_aspect('equal')
            if i > 0: ax.set_yticklabels([])
        else:
            im = ax.imshow(data_step, origin='lower', cmap='RdBu_r', vmin=-vmax_local, vmax=vmax_local)
            ax.set_xlabel('I')
            if i > 0: ax.set_yticklabels([])
            
        time_hours = t_step * 60 / 3600
        ax.set_title(f"T = {time_hours:.1f} h")

    fig.suptitle('Kelvin Wave Propagation (SSH)', fontsize=20)
    # Add shared colorbar
    fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.02, label='SSH (m)')
    
    plt.savefig(os.path.join(script_dir, 'fig_expA_ssh_snapshots.png'), dpi=300, bbox_inches='tight')
    print(f"Saved {os.path.join(script_dir, 'fig_expA_ssh_snapshots.png')}")

    # --- 3. Variance Map (Geo) ---
    print("Generating Variance Map...")
    
    t_start = 120 
    if ssh.shape[0] > t_start:
        ssh_slice = ssh[t_start:, :, :]
    else:
        ssh_slice = ssh[:]
    
    ssh_var = np.var(ssh_slice, axis=0)
    vmax = np.percentile(ssh_var, 99)
    print(f"Variance Map Vmax: {vmax}")

    # Use subplots with constrained_layout
    # Adjusted figsize to be narrower to minimize whitespace for narrow basin
    fig, ax = plt.subplots(figsize=(5, 9), constrained_layout=True)
    
    if 'lon' in locals() and lon is not None:
        im = ax.pcolormesh(lon, lat, ssh_var, cmap='viridis', shading='auto', vmax=vmax)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        if extent:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        ax.set_aspect('equal')
    else:
        im = ax.imshow(ssh_var, origin='lower', cmap='viridis', vmax=vmax)
        ax.set_xlabel('I Index')
        ax.set_ylabel('J Index')
        
    fig.colorbar(im, ax=ax, label='SSH Variance ($m^2$)', shrink=0.8)
    ax.set_title('SSH Variance (Node Identification)')
    
    # bbox_inches='tight' removes extra whitespace around the plot
    plt.savefig(os.path.join(script_dir, 'fig_expA_ssh_variance.png'), dpi=300, bbox_inches='tight')
    print(f"Saved {os.path.join(script_dir, 'fig_expA_ssh_variance.png')}")

    ds.close()

if __name__ == "__main__":
    analyze_wave()
