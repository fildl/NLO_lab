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

    # --- 0. Initial State (T=0) ---
    print("Generating Initial State Plot...")
    t_step_0 = 0
    ssh_t0 = ssh[t_step_0, :, :] * 1000 # Convert to mm
    vmax_init = 0.4
    
    fig_init, ax_init = plt.subplots(figsize=(6, 8), constrained_layout=True)
    if 'lon' in locals() and lon is not None:
        im_init = ax_init.pcolormesh(lon, lat, ssh_t0, cmap='RdBu_r', shading='auto', vmin=-vmax_init, vmax=vmax_init)
        ax_init.set_xlabel('Longitude (°E)')
        ax_init.set_ylabel('Latitude (°N)')
        ax_init.set_aspect('equal')
    else:
        im_init = ax_init.imshow(ssh_t0, origin='lower', cmap='RdBu_r', vmin=-vmax_init, vmax=vmax_init)
        ax_init.set_xlabel('I')
        ax_init.set_ylabel('J')

    ax_init.set_title(f'Initial State (T=0)\nSSH Anomaly (mm)', fontsize=14)
    fig_init.colorbar(im_init, ax=ax_init, label='SSH (mm)', shrink=0.8)
    
    output_path_init = os.path.join(script_dir, 'fig_expA_initial_state.png')
    plt.savefig(output_path_init, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path_init}")

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
    
    # Convert to mm
    hovmoller_mm = hovmoller * 1000
    
    # Use fixed scale for consistency
    h_vmax = 0.4
    
    # 1a. Standard Hovmöller
    plt.imshow(hovmoller_mm, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, ny*10, 0, ssh.shape[0]*60/3600])
    plt.colorbar(label='SSH (mm)')
    plt.xlabel('Distance along Coast (km)')
    plt.ylabel('Time (hours)')
    plt.title('Hovmöller Diagram (East Coast Path)')
    plt.savefig(os.path.join(script_dir, 'fig_expA_hovmoller_east.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expA_hovmoller_east.png')}")

    # 1b. Hovmöller with Theory Line
    print("Generating Hovmöller Diagram with Theory Line...")
    plt.figure(figsize=(10, 8))
    
    # Time/Dist info
    total_time_hours_theory = ssh.shape[0] * 60.0 / 3600
    total_dist_km_theory = ny * 10 
    
    plt.imshow(hovmoller_mm, aspect='auto', origin='lower', cmap='RdBu_r', 
               vmin=-h_vmax, vmax=h_vmax,
               extent=[0, total_dist_km_theory, 0, total_time_hours_theory])
    
    # Theoretical Line: c = sqrt(g * H)
    g_grav = 9.81
    H_depth = 100.0
    c_theory = np.sqrt(g_grav * H_depth) # ~ 31.32 m/s
    
    x_km_theory = np.linspace(0, total_dist_km_theory, 100)
    x_m_theory = x_km_theory * 1000.0
    t_sec_theory = x_m_theory / c_theory
    t_hours_theory = t_sec_theory / 3600.0
    
    plt.plot(x_km_theory, t_hours_theory, 'k--', linewidth=2, label=f'Theory c={c_theory:.1f} m/s')
    
    plt.colorbar(label='SSH (mm)')
    plt.xlabel('Distance along Coast (km)')
    plt.ylabel('Time (hours)')
    plt.title('Hovmöller Diagram (East Coast Path)')
    plt.legend(loc='upper right')
    plt.savefig(os.path.join(script_dir, 'fig_expA_hovmoller_east_theory.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expA_hovmoller_east_theory.png')}")


    # Georeferenced Snapshots
    print("Generating Georeferenced Snapshots...")
    # Plot 5 snapshots from 2h to 16h as requested
    # 2h = 120 steps, 16h = 960 steps
    steps = np.linspace(120, 960, 5, dtype=int)
    
    # Calculate robust Vmax from the selected frames ONLY
    # This ensures colorbar is not saturated by unshown transients
    selected_ssh = ssh[steps, :, :]
    
    # Convert to mm for plotting
    selected_ssh_mm = selected_ssh * 1000
    
    vmax_local = 0.4 # Fixed scale
    print(f"Using Fixed Vmax (mm): {vmax_local}")
    
    # Larger fonts
    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14})
    
    fig, axes = plt.subplots(1, 5, figsize=(15, 8), constrained_layout=True)
    
    # Calculate bounds (re-check if lon is defined)
    if 'lon' not in locals() or lon is None: extent = None
    else:
        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()
        extent = [lon_min, lon_max, lat_min, lat_max]

    for i, t_step in enumerate(steps):
        if t_step >= ssh.shape[0]: continue
        ax = axes[i]
        
        # Get data and convert to mm
        data_step = ssh[t_step, :, :] * 1000
        
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

    # Reduced top spacing by bringing title closer, but ensuring no overlap
    fig.suptitle('Kelvin Wave Propagation (SSH)', fontsize=20, y=1.06)
    
    # Add shared colorbar
    # Pad reduced to bring it closer
    fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.05, pad=0.03, label='SSH (mm)')
    
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
    vmax_var = np.percentile(ssh_var, 99)
    print(f"Variance Map Vmax: {vmax_var}")

    # 3a. Standard Variance Plot
    # Use standardized figsize (5.5, 9)
    fig, ax = plt.subplots(figsize=(5.5, 9), constrained_layout=True)
    
    if 'lon' in locals() and lon is not None:
        im = ax.pcolormesh(lon, lat, ssh_var, cmap='viridis', shading='auto', vmax=vmax_var)
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        if extent:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        ax.set_aspect('equal')
    else:
        im = ax.imshow(ssh_var, origin='lower', cmap='viridis', vmax=vmax_var)
        ax.set_xlabel('I Index')
        ax.set_ylabel('J Index')
        
    fig.colorbar(im, ax=ax, label='SSH Variance ($m^2$)', shrink=0.8, format='%.1e')
    ax.set_title('SSH Variance', fontsize=14)
    
    plt.savefig(os.path.join(script_dir, 'fig_expA_ssh_variance.png'), dpi=300, bbox_inches='tight')
    print(f"Saved {os.path.join(script_dir, 'fig_expA_ssh_variance.png')}")

    # 3b. Variance Contour Plot
    print("Generating Variance Contour Plot...")
    # Use standardized figsize (5.5, 9)
    fig_contour, ax_contour = plt.subplots(figsize=(5.5, 9), constrained_layout=True)
    
    # Ensure 0 is included
    levels = np.linspace(0, vmax_var, 15)

    if 'lon' in locals() and lon is not None:
        # Filled contours
        cf = ax_contour.contourf(lon, lat, ssh_var, levels=levels, cmap='viridis', extend='max')
        # Line contours
        c = ax_contour.contour(lon, lat, ssh_var, levels=levels, colors='white', linewidths=0.5, alpha=0.5)
        ax_contour.set_xlabel('Longitude (°E)')
        ax_contour.set_ylabel('Latitude (°N)')
        ax_contour.set_aspect('equal')
    else:
        cf = ax_contour.contourf(ssh_var, levels=levels, cmap='viridis', extend='max')
        c = ax_contour.contour(ssh_var, levels=levels, colors='white', linewidths=0.5, alpha=0.5)
        ax_contour.set_xlabel('I Index')
        ax_contour.set_ylabel('J Index')

    ax_contour.set_title('SSH Variance', fontsize=14)
    # Use scientific notation for small variance values
    fig_contour.colorbar(cf, ax=ax_contour, label='SSH Variance ($m^2$)', orientation='vertical', shrink=0.8, format='%.1e')
    
    output_path_contour = os.path.join(script_dir, 'fig_expA_ssh_variance_contour.png')
    plt.savefig(output_path_contour, dpi=300, bbox_inches='tight')
    print(f"Saved {output_path_contour}")

    ds.close()

if __name__ == "__main__":
    analyze_wave()
