import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_C_slope_100_sigma.nc'
filepath = os.path.join(work_dir, filename)
mesh_filename = 'mesh_mask_C_sigma.nc' 
# Assuming mesh file is in the same experiments dir or mesh dir. 
# User said "ho caricato... mesh_mask_C_sigma.nc". 
# I'll check usual locations.

# Mesh mask path candidates
meshpath_candidates = [
    os.path.join(script_dir, 'mesh', mesh_filename),
    os.path.join(work_dir, mesh_filename),
]
meshpath = ""
for p in meshpath_candidates:
    if os.path.exists(p):
        meshpath = p
        break

def analyze_exp_c_sigma():
    if not os.path.exists(filepath):
        print(f"Error: Data file not found at {filepath}")
        return

    print(f"Analyzing Experiment C (Sigma): {filename}")
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
    print("Generating Hovmöller Diagram (Sigma)...")
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
    h_vmax = 0.4 # Fixed scale

    plt.imshow(hovmoller * scale_mm, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, ny*10, 0, ssh.shape[0]*dt_c/3600]) # Approx 10km grid
    plt.colorbar(label='SSH (mm)')
    plt.xlabel('Distance South-North (km approx)')
    plt.ylabel('Time (hours)')
    plt.title('Exp C (Sigma): Hovmöller Diagram (Smooth Bathymetry)')
    plt.savefig(os.path.join(script_dir, 'fig_ExpC_hovmoller_sigma.png'), dpi=300)
    print("Saved fig_ExpC_hovmoller_sigma.png")

    # --- 2. Snapshots ---
    print("Generating Snapshots (Sigma)...")
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
    print(f"Exp C Sigma Vmax (mm): {vmax_local}")
    
    # Larger fonts
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14})
    
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
    fig.suptitle('Exp C (Sigma): Kelvin Wave Propagation', fontsize=20, y=1.06)
    if 'im' in locals():
        # Match Exp A colorbar padding
        fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.03, label='SSH (mm)')
    plt.savefig(os.path.join(script_dir, 'fig_ExpC_snapshots_sigma.png'), dpi=300, bbox_inches='tight')
    print("Saved fig_ExpC_snapshots_sigma.png")

    ds.close()

    # --- 3. Shoaling Analysis (Comparison with Baseline) ---
    # Optional: we could compare Z-coord vs Sigma-coord here, but User said "non cancellare i plot vecchi".
    # I will generate specific comparison: Exp A vs Exp C (Sigma)
    print("Generating Shoaling Analysis (Exp A vs C-Sigma)...")
    baseline_file = os.path.join(work_dir, 'EXP_AMP_0.1m.nc') # Use 0.1m as Baseline
    
    if os.path.exists(baseline_file):
        ds_base = netCDF4.Dataset(baseline_file)
        ssh_base = ds_base.variables['sossheig'][:]
        ds_base.close()
        
        # Compare Time Series at North
        j_north, i_north = -5, ssh.shape[2]//2
        if east_coast_idx: i_north = int(east_coast_idx[-5])
        
        ts_c = ssh[:, j_north, i_north]
        ts_a = ssh_base[:, j_north, i_north]
        
        plt.figure(figsize=(10, 6))
        time_c = np.arange(len(ts_c)) * dt_c / 3600
        time_a = np.arange(len(ts_a)) * dt_a / 3600
        
        scale_mm = 1000.0
        plt.plot(time_c, ts_c * scale_mm, 'r-', label='Exp C Sigma (Slope)', linewidth=2)
        plt.plot(time_a, ts_a * scale_mm, 'k--', label='Baseline (Flat)', linewidth=1.5)
        plt.title('Shoaling Effect (Sigma Coord): SSH at Northern Coast')
        plt.xlabel('Time (hours)')
        plt.ylabel('SSH (mm)')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(script_dir, 'fig_ExpC_shoaling_comparison_sigma.png'), dpi=300)
        print("Saved fig_ExpC_shoaling_comparison_sigma.png")
    else:
        print(f"Baseline file {baseline_file} not found. Skipping shoaling comparison.")


    # --- 4. Bathymetry Map (3D Surface) - Sigma ---
    print("Generating 3D Bathymetry Map (Sigma) from Mesh Mask...")
    # NOTE: mesh_mask_C_sigma.nc should show smooth levels if Sigma is active and mesh is correct.
    
    # We used global meshpath found at start function
    mesh_mask_file = meshpath
    
    if os.path.exists(mesh_mask_file):
        try:
            ds_mesh = netCDF4.Dataset(mesh_mask_file)
            # gdepw_0 usually contains the depths
            # In Sigma, gdepw_0[level, y, x] varies with y.
            gdepw = ds_mesh.variables['gdepw_0'][0, :, :, :] 
            mbathy = ds_mesh.variables['mbathy'][0, :, :]
            # For Sigma, mbathy is likely N-1 everywhere (full column)
            
            # Reconstruction logic: We want the Bottom Depth.
            # Usually bottom T-level depth or W-level.
            # Let's take the last W-level which represents the seafloor?
            # Or use mbathy indices.
            
            ny, nx = mbathy.shape
            bathy_2d = np.zeros((ny, nx))
            
            for j in range(ny):
                for i in range(nx):
                    level = mbathy[j, i] # Last WET level index (Fortran index usually? NetCDF is 0-based?)
                    # NEMO masks: mbathy gives the number of wet levels.
                    # So bottom depth is at w-level index = mbathy.
                    # e.g. if mbathy=30, levels 0..29 are wet. w-depth at 30 is bottom.
                    if level > 0:
                        # Safety check index
                        idx = int(level)
                        if idx >= gdepw.shape[0]: idx = gdepw.shape[0] - 1
                        bathy_2d[j, i] = gdepw[idx, j, i]
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
            path_lons = []
            path_lats = []
            path_z = []
            for j in range(ny):
                i_idx = int(east_coast_idx[j])
                if i_idx >= nx: i_idx = nx - 1
                if lon is not None:
                    path_lons.append(lon[j, int(i_idx)])
                    path_lats.append(lat[j, int(i_idx)])
                else:
                    path_lons.append(i_idx)
                    path_lats.append(j)
                
                # Find depth at this point
                d_val = bathy_2d[j, i_idx]
                path_z.append(-d_val + 10) # Lift slightly
            
            ax.plot(path_lons, path_lats, path_z, 'r-', linewidth=3, label='Hovmöller Path', zorder=10)
            
            ax.set_title('Experiment C (Sigma): Bathymetry (Smooth)', fontsize=16)
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
            ax.set_zlabel('Depth (m)')
            ax.legend()
            
            m = cm.ScalarMappable(cmap='viridis')
            m.set_array(bathy_2d)
            fig.colorbar(m, ax=ax, shrink=0.6, label='Depth (m)')
            
            ax.view_init(elev=30, azim=225)
            
            plt.savefig(os.path.join(script_dir, 'fig_ExpC_bathymetry_3D_sigma.png'), dpi=300, bbox_inches='tight')
            print("Saved fig_ExpC_bathymetry_3D_sigma.png")
            
        except Exception as e:
            print(f"Error plotting 3D bathymetry: {e}")
    else:
        print(f"Mesh mask not found at {mesh_mask_file}. Skipping 3D Bathymetry.")

if __name__ == "__main__":
    analyze_exp_c_sigma()
