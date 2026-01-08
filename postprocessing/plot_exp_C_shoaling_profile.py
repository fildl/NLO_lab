
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os

# Settings
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(script_dir, 'experiments')
filename = 'EXP_C_slope_100.nc'
filepath = os.path.join(work_dir, filename)

def plot_shoaling_profile():
    if not os.path.exists(filepath):
        print(f"Error: Data file not found at {filepath}")
        return

    print(f"Analyzing Shoaling Profile for Experiment C: {filename}")
    ds = netCDF4.Dataset(filepath)
    ssh = ds.variables['sossheig'][:]
    
    # Load Dimensions
    nt, ny, nx = ssh.shape
    
    # Load Coords if available to compute distance roughly
    if 'nav_lat' in ds.variables:
        lat = ds.variables['nav_lat'][:]
        # Calculate approximate distance vector (South to North)
        # Assuming ~111km per degree, verify with grid spacing
        # Or just use model grid steps (dy=10km) which is known from config
        dist_km = np.arange(ny) * 10.0
    else:
        dist_km = np.arange(ny) * 10.0 # Default fallback

    # --- 1. Extract Max Amplitude along Propagation Path ---
    # Path: East Coast (same as Hovmöller)
    # i_idx = nx - 2 (Fixed index near right boundary)
    i_idx = nx - 2
    
    amp_observed = np.zeros(ny)
    
    print("Extracting maximum amplitudes along path...")
    for j in range(ny):
        # Time series at this point
        ts = ssh[:, j, i_idx]
        # Max amplitude (peak of the wave)
        # Use simple max(abs) or just max if wave is positive. 
        # Kelvin wave hump is positive. But might have slight negative rebound.
        # We generally care about the max elevation.
        amp_observed[j] = np.max(ts)
        
        amp_observed[j] = np.max(ts)
    
    # Define Depth Profile H(y) for the whole domain first
    h_start = 1000.0
    h_end = 100.0
    h_profile = np.linspace(h_start, h_end, ny)

    # --- 2. Filter / Crop Initial Transient (User Request: ignore first 400km) ---
    # User Request (Iter 2): Cut also at 1000km to remove boundary effects
    start_dist_km = 400.0
    end_dist_km = 1000.0
    mask = (dist_km >= start_dist_km) & (dist_km <= end_dist_km)
    
    dist_km_sub = dist_km[mask]
    amp_observed_sub = amp_observed[mask]
    h_profile_sub = h_profile[mask]
    
    # --- 3. Calculate Theoretical Amplitude (Green's Law) ---
    # Green's Law: A ~ H^(-1/4)  => A(y) = A_ref * (H_ref / H(y))^(1/4)
    # Re-anchor the theory to the first valid observed point after the cut
    
    A_ref = amp_observed_sub[0]
    H_ref = h_profile_sub[0]
    
    # Green's Law Prediction anchored at 100km
    amp_theory_sub = A_ref * (H_ref / h_profile_sub)**0.25
    
    # --- 4. Plotting ---
    plt.figure(figsize=(10, 6))
    
    # Convert to cm for plotting
    plt.plot(dist_km_sub, amp_observed_sub * 100, 'b-', linewidth=2, label='Observed Amplitude')
    plt.plot(dist_km_sub, amp_theory_sub * 100, 'r--', linewidth=2, label="Green's Law (Theory, anchored at 400km)")
    
    # Add Bathymetry scale on twin axis for context
    ax1 = plt.gca()
    ax1.set_xlabel('Distance from South (km)')
    ax1.set_ylabel('SSH Amplitude (cm)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True, alpha=0.3)
    
    ax2 = ax1.twinx()
    ax2.plot(dist_km_sub, h_profile_sub, 'k:', alpha=0.5, label='Depth Profile')
    ax2.set_ylabel('Depth (m)', color='k')
    ax2.set_ylim(0, 700)
    ax2.invert_yaxis() # Depth downwards
    
    # Highlight Start and End
    plt.title('Experiment C: Shoaling Effect Analysis\nAmplitude Evolution vs Depth')
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center')
    
    save_path = os.path.join(script_dir, 'fig_ExpC_shoaling_profile.png')
    plt.savefig(save_path, dpi=300)
    plt.show()
    print(f"Saved Shoaling Profile to {save_path}")
    
    ds.close()

if __name__ == "__main__":
    plot_shoaling_profile()
