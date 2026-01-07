import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_C_slope.nc'
filepath = os.path.join(work_dir, filename)

# Mesh mask path candidates
meshpath_candidates = [
    os.path.join(script_dir, 'mesh', 'mesh_mask_C.nc'),
    os.path.join(work_dir, 'mesh_mask_C.nc'),
    os.path.join(script_dir, 'experiments', 'mesh_mask.nc') # Fallback to generic
]
meshpath = ""
for p in meshpath_candidates:
    if os.path.exists(p):
        meshpath = p
        break

def analyze_exp_c():
    if not os.path.exists(filepath):
        print(f"Error: Data file not found at {filepath}")
        return

    print(f"Analyzing Experiment C: {filename}")
    ds = netCDF4.Dataset(filepath)
    ssh = ds.variables['sossheig'][:]
    
    # --- Robust Mask Loading ---
    tmask = None
    if os.path.exists(meshpath):
        try:
            mds = netCDF4.Dataset(meshpath)
            # Check if tmask exists
            if 'tmask' in mds.variables:
                temp_mask = mds.variables['tmask'][0,0,:,:]
                # Verify dimensions align with data
                if temp_mask.shape == ssh[0,:,:].shape:
                    tmask = temp_mask
                    print(f"Loaded valid tmask from {os.path.basename(meshpath)}")
                else:
                    print(f"Warning: Mask dimensions {temp_mask.shape} mismatch data {ssh[0,:,:].shape}. Ignoring file mask.")
            mds.close()
        except Exception as e:
            print(f"Warning: Failed to load mesh mask: {e}")

    if tmask is None:
        print("Using synthetic fallback mask (all water except 1-pixel boundary).")
        tmask = np.ones_like(ssh[0,:,:])
        tmask[0,:] = 0; tmask[-1,:] = 0
        tmask[:,0] = 0; tmask[:,-1] = 0

    # --- Load Coordinates ---
    lon, lat, extent = None, None, None
    if 'nav_lon' in ds.variables:
        lon = ds.variables['nav_lon'][:]
        lat = ds.variables['nav_lat'][:]
        extent = [lon.min(), lon.max(), lat.min(), lat.max()]
    

    # --- 1. Hovmöller Diagram (East Coast) ---
    print("Generating Hovmöller Diagram...")
    ny, nx = ssh.shape[1], ssh.shape[2]
    
    # Time settings
    dt_c = 20.0  # seconds (Exp C)
    dt_a = 60.0  # seconds (Exp A / Baseline)
    
    # Fixed East Coast Path
    # Use fixed index near eastern boundary (nx-2) to avoid artifacts
    east_coast_idx = [nx - 2] * ny

    hovmoller = np.zeros((ssh.shape[0], ny))
    for t in range(ssh.shape[0]):
        for j in range(ny):
            col = east_coast_idx[j]
            hovmoller[t, j] = ssh[t, j, int(col)]

    plt.figure(figsize=(10, 8))
    h_vmax = np.percentile(np.abs(hovmoller), 99)
    if h_vmax == 0: h_vmax = 0.05
    
    plt.imshow(hovmoller, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, ny*10, 0, ssh.shape[0]*dt_c/3600]) # Approx 10km grid
    plt.colorbar(label='SSH (m)')
    plt.xlabel('Distance South-North (km approx)')
    plt.ylabel('Time (hours)')
    plt.title('Exp C: Hovmöller Diagram (Shoaling Check)')
    plt.savefig(os.path.join(script_dir, 'fig_ExpC_hovmoller.png'), dpi=300)
    print("Saved fig_ExpC_hovmoller.png")

    # --- 2. Snapshots ---
    print("Generating Snapshots...")
    # Plot 5 snapshots from 2h to 16h
    # Need to calculate steps based on dt=20
    # 2h = 7200s = 360 steps
    # 16h = 57600s = 2880 steps
    steps = np.linspace(360, 2880, 5, dtype=int)
    
    # Robust Vmax from selected frames
    # Ensure steps don't exceed file length
    steps = [s for s in steps if s < ssh.shape[0]]
    if not steps: steps = [0]
    
    selected_ssh = ssh[steps, :, :]
    vmax_local = np.percentile(np.abs(selected_ssh), 99.9)
    print(f"Exp C Vmax: {vmax_local}")
    
    # Larger fonts
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14})
    
    fig, axes = plt.subplots(1, 5, figsize=(15, 10), constrained_layout=True)
    
    # Handle single plot case if only 1 step
    if len(steps) == 1: axes = [axes]

    for i, t_step in enumerate(steps):
        ax = axes[i]
        data_step = ssh[t_step, :, :]
        # Apply mask
        data_step = np.ma.masked_where(tmask == 0, data_step)
        
        if lon is not None:
            # Symmetric Vmin/Vmax
            im = ax.pcolormesh(lon, lat, data_step, cmap='RdBu_r', shading='auto', vmin=-vmax_local, vmax=vmax_local)
            if i==0: ax.set_ylabel('Lat (°N)')
            ax.set_xlabel('Lon (°E)')
            ax.set_aspect('equal')
        else:
            im = ax.imshow(data_step, origin='lower', cmap='RdBu_r', vmin=-vmax_local, vmax=vmax_local)
            ax.set_xlabel('I')
            
        ax.set_title(f"T = {t_step*dt_c/3600:.1f} h")
        if i > 0: ax.set_yticklabels([])

    fig.suptitle('Exp C: Kelvin Wave Propagation', fontsize=20)
    if 'im' in locals():
        fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.02, label='SSH (m)')
    plt.savefig(os.path.join(script_dir, 'fig_ExpC_snapshots.png'), dpi=300)
    print("Saved fig_ExpC_snapshots.png")

    ds.close()

    # --- 3. Shoaling Analysis (Comparison with Baseline) ---
    print("Generating Shoaling Analysis...")
    baseline_file = os.path.join(work_dir, 'EXP_AMP_0.1m.nc') # Use 0.1m as Baseline
    
    if os.path.exists(baseline_file):
        ds_base = netCDF4.Dataset(baseline_file)
        ssh_base = ds_base.variables['sossheig'][:]
        ds_base.close()
        
        # Compare Time Series at North
        j_north, i_north = -5, ssh.shape[2]//2
        # Auto-detect better indices if possible
        if east_coast_idx: i_north = int(east_coast_idx[-5])
        
        ts_c = ssh[:, j_north, i_north]
        ts_a = ssh_base[:, j_north, i_north]
        
        plt.figure(figsize=(10, 6))
        time_c = np.arange(len(ts_c)) * dt_c / 3600
        time_a = np.arange(len(ts_a)) * dt_a / 3600
        
        plt.plot(time_c, ts_c, 'r-', label='Exp C (Slope: 30m)', linewidth=2)
        plt.plot(time_a, ts_a, 'k--', label='Baseline (Flat: 100m)', linewidth=1.5)
        plt.title('Shoaling Effect: SSH at Northern Coast')
        plt.xlabel('Time (hours)')
        plt.ylabel('SSH (m)')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(script_dir, 'fig_ExpC_shoaling_comparison.png'), dpi=300)
        print("Saved fig_ExpC_shoaling_comparison.png")
    else:
        print(f"Baseline file {baseline_file} not found. Skipping shoaling comparison.")

    # --- 4. Bathymetry Map (2D) ---
    print("Generating Bathymetry Map...")
    if lon is not None:
        fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
        
        # Reconstruct Depth (Analytical: 1000m South -> 30m North)
        ny, nx = lat.shape
        depth_profile = np.linspace(1000, 30, ny)
        depth_grid = np.tile(depth_profile[:, np.newaxis], (1, nx))
        
        # Plot 2D Map
        # User requested "another colormap" (was terrain). 'viridis' is good for depth.
        cmap = 'viridis_r' # Reversed so deep (1000) is Purple/Dark, shallow (30) is Yellow/Light
        im = ax.pcolormesh(lon, lat, depth_grid, cmap=cmap, shading='auto')
        
        ax.set_title('3D Bathymetry\n(Experiment C)', fontsize=14)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.set_aspect('equal')
        
        # Colorbar with Min/Max
        dmin, dmax = np.min(depth_grid), np.max(depth_grid)
        cbar = fig.colorbar(im, ax=ax, label='Depth (m)')
        cbar.set_ticks([dmin, 250, 500, 750, dmax])
        cbar.set_ticklabels([f'{int(dmin)}m (North)', '250', '500', '750', f'{int(dmax)}m (South)'])
        
        # Annotation removed as per user request
        
        # --- Overlay Extraction Path ---
        # Convert grid indices to Lon/Lat coordinates
        path_lons = []
        path_lats = []
        for j in range(ny):
            i_idx = east_coast_idx[j]
            # Ensure index within bounds
            if i_idx >= nx: i_idx = nx - 1
            path_lons.append(lon[j, int(i_idx)])
            path_lats.append(lat[j, int(i_idx)])
            
        ax.plot(path_lons, path_lats, 'r-', linewidth=2, label='Hovmöller Path')
        ax.legend(loc='upper right')
        
        plt.savefig(os.path.join(script_dir, 'fig_ExpC_bathymetry_3D.png'), dpi=300) # Keep filename as requested by user flow implies replacing
        print("Saved fig_ExpC_bathymetry_3D.png (2D Map version with path)")
    else:
        print("Cannot plot Bathymetry without Longitude/Latitude data.")

    # --- 5. Hovmöller Comparison (A vs C) ---
    print("Generating Hovmöller Comparison (Exp A vs C)...")
    baseline_file = os.path.join(work_dir, 'EXP_AMP_0.1m.nc') # Use 0.1m as Baseline
    
    if os.path.exists(baseline_file):
        ds_base = netCDF4.Dataset(baseline_file)
        ssh_base = ds_base.variables['sossheig'][:]
        ds_base.close()
        
        # Helper to extract hovmoller
        def get_hov(data_ssh):
            ny, nx = data_ssh.shape[1], data_ssh.shape[2]
            # Fixed index path
            ec_idx = [nx - 2] * ny
            
            hov = np.zeros((data_ssh.shape[0], ny))
            for t in range(data_ssh.shape[0]):
                for j in range(ny):
                    hov[t, j] = data_ssh[t, j, int(ec_idx[j])]
            return hov
            
        hov_A = get_hov(ssh_base)
        hov_C = get_hov(ssh) # Recalculate or could reuse if refactored, quick recalc is fine
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 8), constrained_layout=True)
        
        # Common scaling
        vmax = np.percentile(np.abs(hov_A), 99)
        if vmax == 0: vmax = 0.05
    
        # 1. Exp A
        ax = axes[0]
        im = ax.imshow(hov_A, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       extent=[0, hov_A.shape[1]*10, 0, hov_A.shape[0]*dt_a/3600])
        ax.set_title('Exp A: Flat Bottom (100m)\nConstant Speed -> Straight Line')
        ax.set_xlabel('Distance South-North (km approx)')
        ax.set_ylabel('Time (hours)')
        
        # Visual guide roughly tracking the wave
        ax.plot([0, hov_A.shape[1]*10], [0, 18], 'k--', alpha=0.3, label='Linear Ref')
    
        # 2. Exp C
        ax = axes[1]
        vmax_c = np.percentile(np.abs(hov_C), 99)
        # Using symmetric vmin/vmax
        im2 = ax.imshow(hov_C, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-vmax_c, vmax=vmax_c,
                       extent=[0, hov_C.shape[1]*10, 0, hov_C.shape[0]*dt_c/3600])
        ax.set_title('Exp C: Sloping Bottom (1000m -> 30m)\nDecelerating -> Curved Line')
        ax.set_xlabel('Distance South-North (km approx)')
        ax.set_yticklabels([])
        
        fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.02, label='SSH (m)')
        
        plt.savefig(os.path.join(script_dir, 'fig_ExpC_hovmoller_comparison.png'), dpi=300)
        print("Saved fig_ExpC_hovmoller_comparison.png")
    else:
        print("Baseline file needed for Hovmöller Comparison not found.")

if __name__ == "__main__":
    analyze_exp_c()
