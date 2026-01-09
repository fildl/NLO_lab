import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
exp_dir = os.path.join(script_dir, 'experiments')

files = {
    'AMP0.1': ('EXP_AMP_0.1m.nc', 0.1),
    'AMP0.5': ('EXP_AMP_0.5m.nc', 0.5),
    'AMP1.0': ('EXP_AMP_1.0m.nc', 1.0)
}

def compare_runs_20km():
    data = {}
    
    # Load Data
    for key, (fname, amp0) in files.items():
        fpath = os.path.join(exp_dir, fname)
        if not os.path.exists(fpath):
            print(f"Missing {fname} at {fpath}")
            continue
            
        ds = netCDF4.Dataset(fpath)
        ssh = ds.variables['sossheig'][:]
        data[key] = {
            'ssh': ssh,
            'norm_ssh': ssh / amp0, 
            'amp0': amp0
        }
        ds.close()

    # Load Coordinates from the first available file
    lon, lat = None, None
    for key, (fname, _) in files.items():
        fpath = os.path.join(exp_dir, fname)
        if os.path.exists(fpath):
            ds = netCDF4.Dataset(fpath)
            if 'nav_lon' in ds.variables:
                lon = ds.variables['nav_lon'][:]
                lat = ds.variables['nav_lat'][:]
            ds.close()
            break
            
    if lon is None:
        print("Warning: Could not load coordinates. Maps will use indices.")
        # Fallback to indices
        ssh_shape = data[list(data.keys())[0]]['ssh'].shape
        lat, lon = np.meshgrid(np.arange(ssh_shape[1]), np.arange(ssh_shape[2]), indexing='ij')
        extent = None
    else:
        # Calculate bounds for tight plotting
        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()
        extent = [lon_min, lon_max, lat_min, lat_max]

    # Global Style Settings
    plt.rcParams.update({'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 12})

    # --- 1. Linearity Check (Time Series at Northern End) ---
    j_venice = data[list(data.keys())[0]]['ssh'].shape[1] - 5
    i_venice = data[list(data.keys())[0]]['ssh'].shape[2] // 2 
    
    print(f"Monitoring Point (North): J={j_venice}, I={i_venice}")

    # --- 0. Location Map ---
    # Requested similar to fig_expA_path.png
    plt.figure(figsize=(4, 8))
    
    # Create mask based on activity (where waves propagate)
    # This avoids masking water as land if it starts at 0
    ssh_all = data[list(data.keys())[0]]['ssh']
    water_mask = np.any(np.abs(ssh_all) > 1e-9, axis=0)
    
    # Mask the land (False values)
    mask_plot = np.ma.masked_where(~water_mask, np.ones_like(water_mask))
    
    if lon is not None:
        # Plot water with a light color
        plt.pcolormesh(lon, lat, mask_plot, cmap='Blues', vmin=0, vmax=1.5, shading='auto')
        # Plot point
        pt_lon = lon[j_venice, i_venice]
        pt_lat = lat[j_venice, i_venice]
        plt.plot(pt_lon, pt_lat, 'r*', markersize=15, markeredgecolor='black', label='Extraction Point')
        
        plt.xlabel('Longitude (°E)')
        plt.ylabel('Latitude (°N)')
        if extent:
            plt.xlim(extent[0], extent[1])
            plt.ylim(extent[2], extent[3])
    else:
        plt.imshow(mask_plot, origin='lower', cmap='Blues', vmin=0, vmax=1.5)
        plt.plot(i_venice, j_venice, 'r*', markersize=15, markeredgecolor='black', label='Extraction Point')
        plt.xlabel('I Index')
        plt.ylabel('J Index')
        
    plt.title('Exp B: Monitoring Location')
    plt.legend(loc='lower left')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, 'fig_expB_location.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expB_location.png')}")

    plt.figure(figsize=(10, 6))
    for key, d in data.items():
        ts = d['norm_ssh'][:, j_venice, i_venice]
        time_axis = np.arange(len(ts)) * 60 / 3600 # hours
        plt.plot(time_axis, ts, label=f"{key} (Normalized)", linewidth=2)
    
    plt.title('Normalized SSH at Northern Basin End')
    plt.xlabel('Time (hours)')
    plt.ylabel('SSH / A0')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, 'fig_expB_compare_linearity.png'), dpi=300)
    print(f"Saved {os.path.join(script_dir, 'fig_expB_compare_linearity.png')}")

    # --- 2. Difference Map (Non-Linearity) ---
    if 'AMP0.1' in data and 'AMP1.0' in data:
        # Difference between Normalized 1.0m and Normalized 0.1m
        diff = data['AMP1.0']['norm_ssh'] - data['AMP0.1']['norm_ssh']
        diff_var = np.var(diff, axis=0)
        
        plt.figure(figsize=(8, 10))
        plt.pcolormesh(lon, lat, diff_var, cmap='inferno', shading='auto')
        plt.colorbar(label='Variance of Norm. Diff ($m^2$)')
        plt.title('Non-Linearity Map (1.0m vs 0.1m)')
        plt.xlabel('Longitude (°E)')
        plt.ylabel('Latitude (°N)')
        if extent:
            plt.xlim(extent[0], extent[1])
            plt.ylim(extent[2], extent[3])
        plt.axis('equal') # Aspect ratio
        plt.tight_layout()
        plt.savefig(os.path.join(script_dir, 'fig_expB_compare_nonlinear_map.png'), dpi=300, bbox_inches='tight')
        print(f"Saved {os.path.join(script_dir, 'fig_expB_compare_nonlinear_map.png')}")

    # --- 3. Variance Map (Amphidromic Point Search) ---
    ref_run = 'AMP0.5' if 'AMP0.5' in data else 'AMP0.1'
    
    # Calculate Variance EXCLUDING the first 2 hours
    t_start = 120 
    if data[ref_run]['ssh'].shape[0] > t_start:
        ssh_slice = data[ref_run]['ssh'][t_start:, :, :]
        print(f"Calculating variance from step {t_start} onwards")
    else:
        ssh_slice = data[ref_run]['ssh'][:]
        
    ssh_var = np.var(ssh_slice, axis=0)
    
    # Robust scaling: Use 99th percentile
    vmax = np.percentile(ssh_var, 99)
    print(f"Variance Map Vmax: {vmax}")
    
    plt.figure(figsize=(8, 10))
    plt.pcolormesh(lon, lat, ssh_var, cmap='viridis', vmax=vmax, shading='auto')
    plt.colorbar(label='SSH Variance ($m^2$)')
    plt.title(f'Amphidromic Point Check ({ref_run})')
    plt.xlabel('Longitude (°E)')
    plt.ylabel('Latitude (°N)')
    if extent:
        plt.xlim(extent[0], extent[1])
        plt.ylim(extent[2], extent[3])
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, 'fig_expB_compare_amphidromic.png'), dpi=300, bbox_inches='tight') 
    print(f"Saved {os.path.join(script_dir, 'fig_expB_compare_amphidromic.png')}")

if __name__ == "__main__":
    compare_runs_20km()
