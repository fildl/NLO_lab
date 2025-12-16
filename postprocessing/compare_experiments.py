import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
exp_dir = '/Users/filippodiludovico/Library/Mobile Documents/com~apple~CloudDocs/Uni/NLO/03_tides/postprocessing/experiments'
files = {
    'REF': ('EXP_REF_Baseline.nc', 0.05),
    'AMP05': ('EXP_AMP_0.5m.nc', 0.5),
    'AMP10': ('EXP_AMP_1.0m.nc', 1.0)
}
meshpath = os.path.join(exp_dir, '../mesh_mask.nc') # Assuming mesh mask is in postprocessing root

def compare_runs():
    data = {}
    
    # Load Data
    for key, (fname, amp0) in files.items():
        fpath = os.path.join(exp_dir, fname)
        if not os.path.exists(fpath):
            print(f"Missing {fname}")
            continue
            
        ds = netCDF4.Dataset(fpath)
        ssh = ds.variables['sossheig'][:]
        data[key] = {
            'ssh': ssh,
            'norm_ssh': ssh / amp0, # Normalized by initial amplitude
            'amp0': amp0
        }
        ds.close()

    if len(data) < 2:
        print("Not enough data to compare.")
        return

    # Load Mask if available
    try:
        mds = netCDF4.Dataset(meshpath)
        tmask = mds.variables['tmask'][0,0,:,:]
        mds.close()
    except:
        tmask = np.ones_like(data['REF']['ssh'][0])

    # --- 1. Linearity Check (Time Series at "Venice" - North End) ---
    # Find a point at the North end (max J, mid I)
    ny, nx = tmask.shape
    j_venice = ny - 5
    i_venice = nx // 2 
    # Refine point: ensure it is wet
    while tmask[j_venice, i_venice] == 0 and j_venice > 0:
        j_venice -= 1
    
    print(f"Monitoring Point (Venice): J={j_venice}, I={i_venice}")

    plt.figure(figsize=(12, 6))
    for key, d in data.items():
        ts = d['norm_ssh'][:, j_venice, i_venice]
        time_axis = np.arange(len(ts)) * 60 / 3600 # hours
        plt.plot(time_axis, ts, label=f"{key} (Normalized by A0={d['amp0']}m)")
    
    plt.title('Normalized SSH at Northern End (Venice)')
    plt.xlabel('Time (hours)')
    plt.ylabel('SSH / A0')
    plt.grid(True)
    plt.legend()
    plt.savefig('comparison_linearity_timeseries.png')
    print("Saved comparison_linearity_timeseries.png")


    # --- 2. Hovmöller Comparison (East Coast) ---
    print("Generating Hovmoller...")
    # Find East Coast Indices
    east_coast_idx = []
    for j in range(ny):
        found = False
        for i in range(nx-1, -1, -1):
            if tmask[j, i] == 1:
                east_coast_idx.append(i)
                found = True
                break
        if not found:
            # Fallback if no wet point found (shouldn't happen with ones_like)
            east_coast_idx.append(nx-2) # Use nx-2 instead of nx-1 (avoid likely halo/boundary)

    # Check extracted indices
    print(f"East Coast indices sample (mid-basin): {east_coast_idx[len(east_coast_idx)//2]}")

    plt.figure(figsize=(10, 8))
    colors = {'REF': 'black', 'AMP05': 'blue', 'AMP10': 'red'}
    linestyles = {'REF': '-', 'AMP05': '--', 'AMP10': ':'}

    # Plot contours
    for key, d in data.items():
        ssh = d['norm_ssh']
        hov = np.zeros((ssh.shape[0], ny))
        for t in range(ssh.shape[0]):
            for j in range(ny):
                idx = int(east_coast_idx[j])
                hov[t, j] = ssh[t, j, idx]
        
        # Debug data range
        print(f" Run {key} Hovmoller Max: {np.nanmax(hov)}, Min: {np.nanmin(hov)}")
        
        time_axis = np.arange(ssh.shape[0]) * 60 / 3600
        dist_axis = np.arange(ny) * 10 # km
        
        # Use simple plot of crest position or a low threshold contour
        # If signal is clean, contour at 0.05 (normalized) should work
        try:
            plt.contour(dist_axis, time_axis, hov, levels=[0.05, 0.5], 
                       colors=colors[key], linestyles=linestyles[key], linewidths=2)
        except:
            print(f"Could not plot contour for {key}")

    # Add dummy legend
    lines = [plt.Line2D([0], [0], color=c, linestyle=l) for c, l in zip(colors.values(), linestyles.values())]
    plt.legend(lines, list(data.keys()), title="Wave Front (SSH=0.05/0.5 * A0)")
    
    plt.xlabel('Distance along Coast (km)')
    plt.ylabel('Time (hours)')
    plt.title('Wave Front Propagation Comparison')
    plt.savefig('comparison_speed_hovmoller.png')
    print("Saved comparison_speed_hovmoller.png")

    # --- 3. Difference Map (Non-Linearity spatial structure) ---
    # SSH_1.0_norm - SSH_REF_norm
    if 'REF' in data and 'AMP10' in data:
        diff = data['AMP10']['norm_ssh'] - data['REF']['norm_ssh']
        diff_var = np.var(diff, axis=0)
        
        plt.figure(figsize=(6, 8))
        plt.imshow(diff_var, origin='lower', cmap='inferno')
        plt.colorbar(label='Variance of Normalized Difference')
        plt.title('Non-Linearity Map (Where A=1.0 differs from Linear)')
        plt.savefig('comparison_nonlinearity_map.png')
        print("Saved comparison_nonlinearity_map.png")

if __name__ == "__main__":
    compare_runs()
