import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments', 'zsigma_20km')
filename = 'EXP_AMP_0.5m.nc'
filepath = os.path.join(work_dir, filename)
meshpath = os.path.join(work_dir, '..', 'mesh_mask.nc') # Try parent directory

def analyze_wave():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return

    ds = netCDF4.Dataset(filepath)
    
    ssh = ds.variables['sossheig'][:]
    
    # Try to load mask
    if os.path.exists(meshpath):
        mds = netCDF4.Dataset(meshpath)
        tmask = mds.variables['tmask'][0,0,:,:]
        mds.close()
        print("Loaded tmask from mesh file.")
    else:
        print("Warning: mesh_mask.nc not found. Proceeding without land mask.")
        tmask = np.ones_like(ssh[0,:,:]) # All 1s (all water assumption)
    
    # Load Georeferenced Coordinates
    if 'nav_lon' in ds.variables:
        lon = ds.variables['nav_lon'][:]
        lat = ds.variables['nav_lat'][:]
    elif os.path.exists(meshpath):
        mds = netCDF4.Dataset(meshpath)
        lon = mds.variables['glamt'][0,:,:]
        lat = mds.variables['gphit'][0,:,:]
        mds.close()
    else:
        # Fallback to simple indices if absolutely nothing available
        print("Warning: No coordinates found. Using indices.")
        lat, lon = np.meshgrid(np.arange(ssh.shape[1]), np.arange(ssh.shape[2]), indexing='ij')

    # --- 1. Hovmöller Diagram (Index-based) ---
    print("Generating Hovmöller Diagram...")
    ny, nx = ssh.shape[1], ssh.shape[2]
    
    # Robust East Coast Path: Track max variance to find the wave
    # This avoids hitting land or zero-boundary cells
    ssh_var_2d = np.var(ssh, axis=0)
    east_coast_idx = []
    
    for j in range(ny):
        # Search in the eastern half of the domain for the maximum signal
        # This assumes the Kelvin wave is on the East coast
        start_search = nx // 2
        if start_search < nx:
            # Find index of max variance in this row
            local_max_idx = np.argmax(ssh_var_2d[j, start_search:])
            actual_i = start_search + local_max_idx
            east_coast_idx.append(actual_i)
        else:
            east_coast_idx.append(nx-2) # Fallback

    hovmoller = np.zeros((ssh.shape[0], ny))
    for t in range(ssh.shape[0]):
        for j in range(ny):
            idx = east_coast_idx[j]
            hovmoller[t, j] = ssh[t, j, int(idx)]

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
    plt.savefig(os.path.join(script_dir, 'fig_hovmoller_east.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_hovmoller_east.png')}")
    
    # --- 2. Snapshots (Index-based & Georeferenced) ---
    print("Generating Snapshots...")
    times_idx = [0, 240, 720, 1200]
    
    # Dynamic limit
    vmax = np.max(np.abs(ssh)) * 0.8
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
            
    plt.savefig(os.path.join(script_dir, 'fig_ssh_snapshots_index.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_ssh_snapshots_index.png')}")

    # Georeferenced Plot
    print("Generating Georeferenced Snapshots...")
    # Plot 4 snapshots
    steps = [0, 240, 720, 1200]
    fig, axes = plt.subplots(1, 4, figsize=(12, 10), constrained_layout=True)
    
    plt.rcParams.update({'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 12})
    
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
            im = ax.pcolormesh(lon, lat, data_step, cmap='RdBu_r', shading='auto', vmin=-0.1, vmax=0.1)
            ax.set_xlabel('Lon (°E)')
            if i == 0: ax.set_ylabel('Lat (°N)')
            if extent:
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
            # Adjust aspect ratio
            ax.set_aspect('equal')
            if i > 0: ax.set_yticklabels([])
        else:
            im = ax.imshow(data_step, origin='lower', cmap='RdBu_r', vmin=-0.1, vmax=0.1)
            ax.set_xlabel('I')
            
        time_hours = t_step * 60 / 3600
        ax.set_title(f"T = {time_hours:.1f} h")

    fig.suptitle('Kelvin Wave Propagation (SSH)', fontsize=16)
    # Add shared colorbar
    fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.02, label='SSH (m)')
    
    plt.savefig(os.path.join(script_dir, 'fig_ssh_snapshots.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_ssh_snapshots.png')}")

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
    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    
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
    
    plt.savefig(os.path.join(script_dir, 'fig_ssh_variance.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_ssh_variance.png')}")

    ds.close()
    if 'mds' in locals() and mds: mds.close()

if __name__ == "__main__":
    analyze_wave()
