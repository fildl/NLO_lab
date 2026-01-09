import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')


def analyze_dataset(config):
    filename = config['filename']
    filepath = os.path.join(work_dir, filename)
    suffix = config['suffix']
    
    if not os.path.exists(filepath):
        print(f"Skipping {filename}: File not found at {filepath}")
        return

    print(f"Analyzing Experiment C: {filename} (Suffix: '{suffix}')")
    ds = netCDF4.Dataset(filepath)
    ssh = ds.variables['sossheig'][:]
    
    # --- Robust Mask Loading ---
    tmask = None
    meshpath = ""
    for p in config['mesh_candidates']:
        if os.path.exists(p):
            meshpath = p
            break
            
    if meshpath and os.path.exists(meshpath):
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
    scale_mm = 1000.0
    scale_mm = 1000.0
    h_vmax = 0.4 # Fixed scale

    plt.imshow(hovmoller * scale_mm, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, ny*10, 0, ssh.shape[0]*dt_c/3600]) # Approx 10km grid
    plt.colorbar(label='SSH (mm)')
    plt.xlabel('Distance South-North (km approx)')
    plt.ylabel('Time (hours)')
    plt.title(f'Exp C ({filename}): Hovmöller Diagram')
    plt.savefig(os.path.join(script_dir, f'fig_ExpC_hovmoller{suffix}.png'), dpi=300)
    print(f"Saved fig_ExpC_hovmoller{suffix}.png")

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
    scale_mm = 1000.0
    vmax_local = 0.4 # Fixed scale
    
    # Larger fonts
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14})
    
    # Reduced height to minimize whitespace (similar to Exp A refinement)
    # Matching Exp A exactly: (15, 8) with constrained_layout
    fig, axes = plt.subplots(1, 5, figsize=(15, 8), constrained_layout=True)
    
    # Handle single plot case if only 1 step
    if len(steps) == 1: axes = [axes]

    for i, t_step in enumerate(steps):
        ax = axes[i]
        data_step = ssh[t_step, :, :] * scale_mm
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

    # Match Exp A title position
    fig.suptitle(f'Exp C ({filename}): Kelvin Wave Propagation', fontsize=20, y=1.06)
    if 'im' in locals():
        # Match Exp A colorbar padding
        fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.03, label='SSH (mm)')
    plt.savefig(os.path.join(script_dir, f'fig_ExpC_snapshots{suffix}.png'), dpi=300, bbox_inches='tight')
    print(f"Saved fig_ExpC_snapshots{suffix}.png")

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
        
        scale_mm = 1000.0
        plt.plot(time_c, ts_c * scale_mm, 'r-', label=f'Exp C {suffix} (Slope: 100m)', linewidth=2)
        plt.plot(time_a, ts_a * scale_mm, 'k--', label='Baseline (Flat: 100m)', linewidth=1.5)
        plt.title(f'Shoaling Effect: SSH at Northern Coast ({filename})')
        plt.xlabel('Time (hours)')
        plt.ylabel('SSH (mm)')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(script_dir, f'fig_ExpC_shoaling_comparison{suffix}.png'), dpi=300)
        print(f"Saved fig_ExpC_shoaling_comparison{suffix}.png")
    else:
        print(f"Baseline file {baseline_file} not found. Skipping shoaling comparison.")

    # --- 4a. Bathymetry Map (2D) ---
    print("Generating 2D Bathymetry Map...")
    if lon is not None:
        fig_2d, ax_2d = plt.subplots(figsize=(8, 10), constrained_layout=True)
        
        # Use simple reconstruction for 2D map (or mesh mask)
        ny, nx = lat.shape
        depth_profile = np.linspace(1000, 100, ny)
        depth_grid = np.tile(depth_profile[:, np.newaxis], (1, nx))
        
        cmap = 'viridis_r'
        im = ax_2d.pcolormesh(lon, lat, depth_grid, cmap=cmap, shading='auto')
        
        ax_2d.set_title(f'Experiment C Bathymetry (2D) {suffix}', fontsize=14)
        ax_2d.set_xlabel('Longitude (°E)')
        ax_2d.set_ylabel('Latitude (°N)')
        ax_2d.set_aspect('equal')
        
        cbar = fig_2d.colorbar(im, ax=ax_2d, label='Depth (m)')
        dmin, dmax = np.min(depth_grid), np.max(depth_grid)
        cbar.set_ticks([dmin, 250, 500, 750, dmax])
        
        # Overlay Path
        path_lons = []
        path_lats = []
        for j in range(ny):
            i_idx = east_coast_idx[j]
            if i_idx >= nx: i_idx = nx - 1
            path_lons.append(lon[j, int(i_idx)])
            path_lats.append(lat[j, int(i_idx)])
            
        ax_2d.plot(path_lons, path_lats, 'r-', linewidth=2, label='Hovmöller Path')
        ax_2d.legend(loc='upper right')
        
        plt.savefig(os.path.join(script_dir, f'fig_ExpC_bathymetry_2D{suffix}.png'), dpi=300, bbox_inches='tight')
        print(f"Saved fig_ExpC_bathymetry_2D{suffix}.png")

    # --- 4b. Bathymetry Map (3D Surface) ---
    print("Generating 3D Bathymetry Map from Mesh Mask...")
    
    if meshpath and os.path.exists(meshpath):
        try:
            ds_mesh = netCDF4.Dataset(meshpath)
            gdepw = ds_mesh.variables['gdepw_0'][0, :, :, :] 
            mbathy = ds_mesh.variables['mbathy'][0, :, :]
            ny, nx = mbathy.shape
            bathy_2d = np.zeros((ny, nx))
            
            for j in range(ny):
                for i in range(nx):
                    level = mbathy[j, i]
                    if level > 0:
                        bathy_2d[j, i] = gdepw[int(level), j, i]
                    else:
                        bathy_2d[j, i] = 0.0
            ds_mesh.close()
            
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib import cm
            
            fig = plt.figure(figsize=(10, 8), constrained_layout=True)
            ax = fig.add_subplot(111, projection='3d')
            
            if lon is not None:
                X, Y = lon, lat
            else:
                X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
            
            Z = -bathy_2d 
            
            # Plot Surface
            surf = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=False, alpha=0.9)
            
            # Add Red Path Line in 3D
            # Need (x, y, z) for the line. 
            # We have path_lons (x) and path_lats (y). Need to lookup depth (z) at those indices.
            path_z = []
            for j in range(ny):
                i_idx = int(east_coast_idx[j])
                if i_idx >= nx: i_idx = nx - 1
                # Find depth at this point
                d_val = bathy_2d[j, i_idx]
                path_z.append(-d_val + 10) # Lift slightly so it is visible above surface
            
            ax.plot(path_lons, path_lats, path_z, 'r-', linewidth=3, label='Hovmöller Path', zorder=10)
            
            ax.set_title(f'Experiment C Bathymetry', fontsize=16)
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
            ax.set_zlabel('Depth (m)')
            ax.legend()
            
            m = cm.ScalarMappable(cmap='viridis')
            m.set_array(bathy_2d)
            fig.colorbar(m, ax=ax, shrink=0.6, label='Depth (m)')
            
            ax.view_init(elev=30, azim=225)
            
            plt.savefig(os.path.join(script_dir, f'fig_ExpC_bathymetry_3D{suffix}.png'), dpi=300, bbox_inches='tight')
            print(f"Saved fig_ExpC_bathymetry_3D{suffix}.png")
            
        except Exception as e:
            print(f"Error plotting 3D bathymetry: {e}")
    else:
        print(f"Mesh mask not found at {meshpath}. Skipping 3D Bathymetry.")

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
        vmax = 0.4 # Fixed scale
        print(f"Comparison Vmax (mm): {vmax}")
    
        # 1. Exp A
        ax = axes[0]
        im = ax.imshow(hov_A * 1000, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       extent=[0, hov_A.shape[1]*10, 0, hov_A.shape[0]*dt_a/3600])
        ax.set_title('Exp A: Flat Bottom (100m)\nConstant Speed -> Straight Line')
        ax.set_xlabel('Distance South-North (km approx)')
        ax.set_ylabel('Time (hours)')
        
        # Visual guide roughly tracking the wave
        # c = sqrt(9.8 * 100) = 31.3 m/s = 112 km/h. Dist ~ 850km -> T ~ 7.5h
        ax.plot([0, hov_A.shape[1]*10], [0, hov_A.shape[1]*10000 / 31.3 / 3600], 'k--', alpha=0.5, label='Theory c=31m/s')
        ax.legend()
    
        # 2. Exp C
        ax = axes[1]
        vmax_c = 0.4 # Fixed scale
        # Using symmetric vmin/vmax
        im2 = ax.imshow(hov_C * 1000, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-vmax_c, vmax=vmax_c,
                       extent=[0, hov_C.shape[1]*10, 0, hov_C.shape[0]*dt_c/3600])
        ax.set_title(f'Exp C ({suffix}): Sloping Bottom\nDecelerating -> Curved Line')
        ax.set_xlabel('Distance South-North (km approx)')
        ax.set_yticklabels([])
        
        # Theoretical Curve for Exp C (Slope)
        # H(y) = 1000 - (1000-100) * y / L
        # c(y) = sqrt(g * H(y))
        # t(y) = Integral dy / c(y)
        L = hov_C.shape[1] * 10000.0 # Length in meters
        y_pts = np.linspace(0, L, 100)
        h_pts = 1000.0 - (1000.0 - 100.0) * y_pts / L
        c_pts = np.sqrt(9.81 * h_pts)
        t_pts = np.cumsum(1.0 / c_pts) * (y_pts[1]-y_pts[0]) # Simple integration
        t_hrs = t_pts / 3600.0
        
        # Plot Theory Curve
        ax.plot(y_pts / 1000.0, t_hrs, 'k--', linewidth=2, label='Theory (variable depth)')
        ax.legend()
        
        fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.02, label='SSH (mm)')
        
        plt.savefig(os.path.join(script_dir, f'fig_ExpC_hovmoller_comparison{suffix}.png'), dpi=300)
        print(f"Saved fig_ExpC_hovmoller_comparison{suffix}.png")
    else:
        print("Baseline file needed for Hovmöller Comparison not found.")

def main():
    analyze_exp_c()

if __name__ == "__main__":
    # Definition of configurations to process
    configurations = [
        {
            'filename': 'EXP_C_slope_100.nc',
            'mesh_candidates': [
                os.path.join(script_dir, 'mesh', 'mesh_mask_C.nc'),
                os.path.join(work_dir, 'mesh_mask_C.nc'),
                os.path.join(script_dir, 'experiments', 'mesh_mask.nc')
            ],
            'suffix': ''
        },
        {
            'filename': 'EXP_C_slope_100_1core.nc',
            'mesh_candidates': [
                os.path.join(script_dir, 'mesh', 'mesh_mask_C_1core.nc'),
                os.path.join(work_dir, 'mesh_mask_C_1core.nc')
            ],
            'suffix': '_1core'
        }
    ]

    for config in configurations:
        analyze_dataset(config)
