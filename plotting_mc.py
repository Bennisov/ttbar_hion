import numpy as np
import matplotlib.pyplot as plt

# Load data
print("Loading data...")
data = np.load('analysis_output.npz', allow_pickle=True)

# Define channels
channels = {
    'el': 'Dielectron',
    'elmu': 'Electron-Muon',
    'mu': 'Dimuon',
    'mu1jet_mu': 'Muon+Jets',
    'mu1jet_el': 'Electron+Jets'
}

# Plot for each channel
for prefix, name in channels.items():
    gen_met_key = f'{prefix}_gen_met'
    gen_met_phi_key = f'{prefix}_gen_met_phi'
    
    # Skip if no generator MET
    if gen_met_key not in data or len(data[gen_met_key]) == 0:
        print(f"Skipping {name} - no generator MET data")
        continue
    
    gen_met = np.asarray(data[gen_met_key]).ravel() / 1000.0  # Convert to GeV
    gen_met_phi = np.asarray(data[gen_met_phi_key]).ravel() if gen_met_phi_key in data else None
    
    # Create figure with 5 subplots for magnitude (one for each MET definition)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    for i in range(5):
        met_key = f'{prefix}_met{i+1}_pts'
        
        if met_key in data and len(data[met_key]) > 0:
            reco_met = np.asarray(data[met_key]).ravel()
            
            # Match array lengths
            min_len = min(len(reco_met), len(gen_met))
            reco_met = reco_met[:min_len]
            gen_met_plot = gen_met[:min_len]
            
            # Plot
            ax = axes[i]
            bins = np.linspace(0, max(np.max(reco_met), np.max(gen_met_plot)), 50)
            ax.hist(gen_met_plot, bins=bins, alpha=0.5, label='Gen MET', 
                   color='red', histtype='step', linewidth=2)
            ax.hist(reco_met, bins=bins, alpha=0.5, label=f'MET{i+1}', 
                   color='blue', histtype='step', linewidth=2)
            
            ax.set_xlabel('MET [GeV]')
            ax.set_ylabel('Events')
            ax.set_title(f'MET{i+1} vs Generator MET')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    # Remove empty subplot
    axes[5].axis('off')
    
    fig.suptitle(f'{name} Channel: MET Magnitude Comparisons', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'met_comparison_{prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved magnitude plot for {name}")
    
    # Create figure for angular (phi) comparison
    if gen_met_phi is not None:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.ravel()
        
        for i in range(5):
            met_phi_key = f'{prefix}_met{i+1}_etas'  # Note: stored as 'etas' but actually phi
            
            if met_phi_key in data and len(data[met_phi_key]) > 0:
                reco_met_phi = np.asarray(data[met_phi_key]).ravel()
                
                # Match array lengths
                min_len = min(len(reco_met_phi), len(gen_met_phi))
                reco_met_phi = reco_met_phi[:min_len]
                gen_met_phi_plot = gen_met_phi[:min_len]
                
                # Plot
                ax = axes[i]
                bins = np.linspace(-np.pi, np.pi, 50)
                ax.hist(gen_met_phi_plot, bins=bins, alpha=0.5, label='Gen MET φ', 
                       color='red', histtype='step', linewidth=2)
                ax.hist(reco_met_phi, bins=bins, alpha=0.5, label=f'MET{i+1} φ', 
                       color='blue', histtype='step', linewidth=2)
                
                ax.set_xlabel('MET φ [rad]')
                ax.set_ylabel('Events')
                ax.set_title(f'MET{i+1} φ vs Generator MET φ')
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        # Remove empty subplot
        axes[5].axis('off')
        
        fig.suptitle(f'{name} Channel: MET Angular (φ) Comparisons', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'met_phi_comparison_{prefix}.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved angular plot for {name}")
    
    # Create figure for MET magnitude differences (resolution)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    for i in range(5):
        met_key = f'{prefix}_met{i+1}_pts'
        
        if met_key in data and len(data[met_key]) > 0:
            reco_met = np.asarray(data[met_key]).ravel()
            
            # Match array lengths
            min_len = min(len(reco_met), len(gen_met))
            reco_met = reco_met[:min_len]
            gen_met_plot = gen_met[:min_len]
            
            # Calculate difference (resolution)
            diff = reco_met - gen_met_plot
            mean_diff = np.mean(diff)
            std_diff = np.std(diff)
            
            # Plot
            ax = axes[i]
            ax.hist(diff, bins=50, alpha=0.7, color='purple', edgecolor='black')
            ax.axvline(mean_diff, color='red', linestyle='--', linewidth=2, 
                      label=f'Mean: {mean_diff:.2f} GeV')
            ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
            
            ax.set_xlabel('MET - Gen MET [GeV]')
            ax.set_ylabel('Events')
            ax.set_title(f'MET{i+1} Resolution (μ={mean_diff:.2f}, σ={std_diff:.2f})')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    # Remove empty subplot
    axes[5].axis('off')
    
    fig.suptitle(f'{name} Channel: MET Resolution (Reco - Gen)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'met_resolution_{prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved resolution plot for {name}")
    
    # Create figure for MET phi differences
    if gen_met_phi is not None:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.ravel()
        
        for i in range(5):
            met_phi_key = f'{prefix}_met{i+1}_etas'
            
            if met_phi_key in data and len(data[met_phi_key]) > 0:
                reco_met_phi = np.asarray(data[met_phi_key]).ravel()
                
                # Match array lengths
                min_len = min(len(reco_met_phi), len(gen_met_phi))
                reco_met_phi = np.asarray(reco_met_phi[:min_len], dtype=float)
                gen_met_phi_plot = np.asarray(gen_met_phi[:min_len], dtype=float)
                
                # Calculate phi difference (handle wraparound)
                diff_phi = reco_met_phi - gen_met_phi_plot
                diff_phi = np.arctan2(np.sin(diff_phi), np.cos(diff_phi))  # Wrap to [-π, π]
                mean_diff_phi = np.mean(diff_phi)
                std_diff_phi = np.std(diff_phi)
                
                # Plot
                ax = axes[i]
                ax.hist(diff_phi, bins=50, alpha=0.7, color='orange', edgecolor='black')
                ax.axvline(mean_diff_phi, color='red', linestyle='--', linewidth=2, 
                          label=f'Mean: {mean_diff_phi:.3f} rad')
                ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
                
                ax.set_xlabel('Δφ (MET - Gen MET) [rad]')
                ax.set_ylabel('Events')
                ax.set_title(f'MET{i+1} φ Resolution (μ={mean_diff_phi:.3f}, σ={std_diff_phi:.3f})')
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        # Remove empty subplot
        axes[5].axis('off')
        
        fig.suptitle(f'{name} Channel: MET φ Resolution (Reco - Gen)', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'met_phi_resolution_{prefix}.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved φ resolution plot for {name}")

print("\nDone!")