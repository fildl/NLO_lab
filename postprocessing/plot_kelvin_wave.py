import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
work_dir = '/Users/filippodiludovico/Library/Mobile Documents/com~apple~CloudDocs/Uni/NLO/03_tides/postprocessing/experiments'
filename = 'EXP_REF_Baseline.nc'
filepath = os.path.join(work_dir, filename)
# meshpath = os.path.join(work_dir, 'mesh_mask_REF.nc') # Temporarily commented out if missing

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
    east_coast_idx = []
    
    for j in range(ny):
        found = False
        for i in range(nx-1, -1, -1):
            if tmask[j, i] == 1:
                east_coast_idx.append(i)
                found = True
                break
        if not found:
            east_coast_idx.append(np.nan)

    hovmoller = np.zeros((ssh.shape[0], ny))
    for t in range(ssh.shape[0]):
        for j in range(ny):
            idx = east_coast_idx[j]
            if not np.isnan(idx):
                hovmoller[t, j] = ssh[t, j, int(idx)]
            else:
                hovmoller[t, j] = np.nan

    plt.figure(figsize=(10, 8))
    plt.imshow(hovmoller, aspect='auto', origin='lower', cmap='RdBu_r', 
               extent=[0, ny*10, 0, ssh.shape[0]*60/3600])
    plt.colorbar(label='SSH (m)')
    plt.xlabel('Distance along Coast (km approx)')
    plt.ylabel('Time (hours)')
    plt.title('Hovmöller Diagram (East Coast)')
    plt.savefig('hovmoller_east.png')
    print("Saved hovmoller_east.png")
    
    # --- 2. Snapshots (Index-based & Georeferenced) ---
    print("Generating Snapshots...")
    times_idx = [0, 200, 500, 900]
    
    # Dynamic limit
    vmax = np.max(np.abs(ssh)) * 0.8
    if vmax < 1e-4: vmax = 0.05

    # Index-based Plot
    fig, axes = plt.subplots(1, 4, figsize=(16, 5))
    for i, tidx in enumerate(times_idx):
        if tidx < ssh.shape[0]:
            ax = axes[i]
            data = ssh[tidx, :, :]
            data = np.ma.masked_where(tmask == 0, data)
            im = ax.imshow(data, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
            ax.set_title(f'T = {tidx*60/3600:.1f} h')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig('ssh_snapshots_index.png')
    print("Saved ssh_snapshots_index.png")

    # Georeferenced Plot
    fig, axes = plt.subplots(1, 4, figsize=(18, 6), sharey=True)
    for i, tidx in enumerate(times_idx):
        if tidx < ssh.shape[0]:
            ax = axes[i]
            data = ssh[tidx, :, :]
            data = np.ma.masked_where(tmask == 0, data)
            # Use pcolormesh for georeferencing
            pcm = ax.pcolormesh(lon, lat, data, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
            
            # Coastline removed as requested
            # ax.contour(lon, lat, tmask, levels=[0.5], colors='black', linewidths=1.5)
            
            ax.set_title(f'T = {tidx*60/3600:.1f} h')
            ax.set_xlabel('Longitude (°E)')
            if i == 0: ax.set_ylabel('Latitude (°N)')
            ax.grid(True, linestyle='--', alpha=0.5)
    
    # Common colorbar for geo plot
    # Create a dummy axis for colorbar if needed, or stick to raveled axes
    # The previous method using axes.ravel().tolist() is robust
    cb = plt.colorbar(pcm, ax=axes.ravel().tolist(), fraction=0.046, pad=0.04)
    cb.set_label('SSH (m)')
    plt.savefig('ssh_snapshots_geo.png')
    print("Saved ssh_snapshots_geo.png")

    # --- 3. Variance Map (Index-based & Georeferenced) ---
    print("Generating Variance Maps...")
    ssh_var = np.var(ssh, axis=0)
    ssh_var = np.ma.masked_where(tmask == 0, ssh_var)
    
    # Index-based
    plt.figure(figsize=(6, 8))
    plt.imshow(ssh_var, origin='lower', cmap='viridis')
    plt.colorbar(label='SSH Variance (m^2)')
    plt.title('SSH Variance (Index coordinates)')
    plt.savefig('ssh_variance_index.png')
    print("Saved ssh_variance_index.png")
    
    # Georeferenced
    plt.figure(figsize=(8, 10))
    pcm = plt.pcolormesh(lon, lat, ssh_var, cmap='viridis', shading='auto')
    # plt.contour(lon, lat, tmask, levels=[0.5], colors='black', linewidths=1.5) # Coastline removed
    plt.colorbar(pcm, label='SSH Variance (m^2)')
    plt.title('SSH Variance (Minima = Potential Nodes)')
    plt.xlabel('Longitude (°E)')
    plt.ylabel('Latitude (°N)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig('ssh_variance_geo.png')
    print("Saved ssh_variance_geo.png")

    ds.close()
    mds.close()

if __name__ == "__main__":
    analyze_wave()
