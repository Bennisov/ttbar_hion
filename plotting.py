import numpy as np
import matplotlib.pyplot as plt

def plot_met_distributions(gen_met, reco_met, gen_phi, reco_phi, 
                           channel_name, prefix=''):
    """Plot both magnitude and angular distributions for MET."""
    
    # Filter out invalid events
    valid_mask = (gen_met >= 0) & (reco_met >= 0)
    gen_met = gen_met[valid_mask]
    reco_met = reco_met[valid_mask]
    gen_phi = gen_phi[valid_mask]
    reco_phi = reco_phi[valid_mask]
    
    if len(gen_met) == 0:
        print(f"No valid events for {channel_name}")
        return
    
    # Create figure with 2 subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # --- Left plot: MET Magnitude ---
    bins_met = np.linspace(0, 100, 50)
    
    ax1.hist(reco_met, bins=bins_met, histtype='step', 
             label='Reco MET', color='blue', linewidth=2)
    ax1.hist(gen_met, bins=bins_met, histtype='step', 
             label='Gen MET (Neutrinos)', color='red', linestyle='--', linewidth=2)
    
    ax1.set_xlabel("Missing Transverse Energy [GeV]", fontsize=12)
    ax1.set_ylabel("Events", fontsize=12)
    ax1.set_title(f"{channel_name}: MET Magnitude", fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)
    
    # Add statistics
    mean_gen = np.mean(gen_met)
    mean_reco = np.mean(reco_met)
    ax1.text(0.65, 0.85, f'Gen Mean: {mean_gen:.1f} GeV\nReco Mean: {mean_reco:.1f} GeV',
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # --- Right plot: MET Phi (Angular Distribution) ---
    bins_phi = np.linspace(-np.pi, np.pi, 50)
    
    ax2.hist(reco_phi, bins=bins_phi, histtype='step', 
             label='Reco MET φ', color='blue', linewidth=2)
    ax2.hist(gen_phi, bins=bins_phi, histtype='step', 
             label='Gen MET φ (Neutrinos)', color='red', linestyle='--', linewidth=2)
    
    ax2.set_xlabel("MET Azimuthal Angle φ [rad]", fontsize=12)
    ax2.set_ylabel("Events", fontsize=12)
    ax2.set_title(f"{channel_name}: MET Angular Distribution", fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(-np.pi, np.pi)
    
    # Add π labels on x-axis
    ax2.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    ax2.set_xticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
    
    plt.tight_layout()
    
    filename = f'{prefix}met_distributions_{channel_name.lower().replace(" ", "_").replace("+", "")}.png'
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close()


def plot_met_phi_correlation(gen_phi, reco_phi, channel_name, prefix=''):
    """Plot 2D correlation between Gen and Reco MET phi angles."""
    
    valid_mask = (gen_phi != 0) | (reco_phi != 0)  # Keep events with any phi info
    gen_phi = gen_phi[valid_mask]
    reco_phi = reco_phi[valid_mask]
    
    if len(gen_phi) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # 2D histogram
    h = ax.hist2d(gen_phi, reco_phi, bins=40, range=[[-np.pi, np.pi], [-np.pi, np.pi]], 
                  cmap='viridis', cmin=1)
    
    # Add diagonal line (perfect agreement)
    ax.plot([-np.pi, np.pi], [-np.pi, np.pi], 'r--', linewidth=2, label='Perfect Agreement')
    
    ax.set_xlabel("Gen MET φ [rad]", fontsize=12)
    ax.set_ylabel("Reco MET φ [rad]", fontsize=12)
    ax.set_title(f"{channel_name}: Gen vs Reco MET φ Correlation", fontsize=14, fontweight='bold')
    
    # Add π labels
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    ax.set_xticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
    ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    ax.set_yticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
    
    ax.legend(fontsize=10)
    plt.colorbar(h[3], ax=ax, label='Events')
    plt.tight_layout()
    
    filename = f'{prefix}met_phi_correlation_{channel_name.lower().replace(" ", "_").replace("+", "")}.png'
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close()


def main():
    try:
        data = np.load('analysis_output.npz', allow_pickle=True)
        print("Loaded analysis_output.npz")
        print(f"Available keys: {list(data.keys())}")
    except FileNotFoundError:
        print("Error: analysis_output.npz not found. Run analysis script first.")
        return
    
    # Channel mapping
    channels = {
        'diel': 'Dielectron',
        'dimu': 'Dimuon',
        'elmu': 'Electron-Muon',
        'mu1j': 'Muon+≥2Jets',
        'el1j': 'Electron+≥2Jets'
    }
    
    print("\nGenerating plots for each channel...")
    
    for ch_key, ch_name in channels.items():
        met_key = f'{ch_key}_reco_met'
        gen_key = f'{ch_key}_gen_met'
        met_phi_key = f'{ch_key}_reco_met_phi'
        gen_phi_key = f'{ch_key}_gen_met_phi'
        
        if all(key in data for key in [met_key, gen_key, met_phi_key, gen_phi_key]):
            print(f"\nProcessing {ch_name}...")
            
            # Plot magnitude and angular distributions
            plot_met_distributions(
                data[gen_key], 
                data[met_key],
                data[gen_phi_key],
                data[met_phi_key],
                ch_name
            )
            
            # Plot phi correlation
            plot_met_phi_correlation(
                data[gen_phi_key],
                data[met_phi_key],
                ch_name
            )
        else:
            print(f"Missing data for {ch_name}")
    
    # Create summary plot with all channels (magnitude only)
    print("\nCreating summary plot...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    bins_met = np.linspace(0, 100, 40)
    
    for idx, (ch_key, ch_name) in enumerate(channels.items()):
        met_key = f'{ch_key}_reco_met'
        gen_key = f'{ch_key}_gen_met'
        
        if met_key in data and gen_key in data:
            gen_met = data[gen_key]
            reco_met = data[met_key]
            
            valid_mask = (gen_met >= 0) & (reco_met >= 0)
            gen_met = gen_met[valid_mask]
            reco_met = reco_met[valid_mask]
            
            if len(gen_met) > 0:
                axes[idx].hist(reco_met, bins=bins_met, histtype='step',
                              label='Reco MET', color='blue', linewidth=1.5)
                axes[idx].hist(gen_met, bins=bins_met, histtype='step',
                              label='Gen MET', color='red', linestyle='--', linewidth=1.5)
                
                axes[idx].set_xlabel("MET [GeV]", fontsize=10)
                axes[idx].set_ylabel("Events", fontsize=10)
                axes[idx].set_title(f"{ch_name} (N={len(gen_met)})", fontsize=11, fontweight='bold')
                axes[idx].legend(fontsize=9)
                axes[idx].grid(alpha=0.3)
    
    # Hide the last subplot if we have fewer than 6 channels
    if len(channels) < 6:
        axes[-1].axis('off')
    
    plt.suptitle("MET Distributions: All Channels", fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('met_summary_all_channels.png', dpi=150)
    print("Saved met_summary_all_channels.png")
    plt.close()
    
    print("\n✓ All plots generated successfully!")


if __name__ == "__main__":
    main()