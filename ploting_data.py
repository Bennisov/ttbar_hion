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

# Plot lepton pT for all channels
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

for idx, (prefix, name) in enumerate(channels.items()):
    pt_key = f'{prefix}_pts'
    if pt_key in data and len(data[pt_key]) > 0:
        pts = data[pt_key]
        axes[idx].hist(pts, bins=50, range=(0, 100), alpha=0.7, color='blue', edgecolor='black')
        axes[idx].set_xlabel('Lepton pT [GeV]')
        axes[idx].set_ylabel('Events')
        axes[idx].set_title(f'{name} - Lepton pT')
        axes[idx].grid(True, alpha=0.3)

axes[5].axis('off')
plt.tight_layout()
plt.savefig('lepton_pt_all_channels.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved lepton pT plot")

# Plot lepton eta for all channels
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

for idx, (prefix, name) in enumerate(channels.items()):
    eta_key = f'{prefix}_etas'
    if eta_key in data and len(data[eta_key]) > 0:
        etas = data[eta_key]
        axes[idx].hist(etas, bins=50, range=(-5, 5), alpha=0.7, color='green', edgecolor='black')
        axes[idx].set_xlabel('Lepton η')
        axes[idx].set_ylabel('Events')
        axes[idx].set_title(f'{name} - Lepton η')
        axes[idx].grid(True, alpha=0.3)

axes[5].axis('off')
plt.tight_layout()
plt.savefig('lepton_eta_all_channels.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved lepton eta plot")

# Plot jet multiplicity for all channels
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

for idx, (prefix, name) in enumerate(channels.items()):
    jet_key = f'{prefix}_jet_mults'
    if jet_key in data and len(data[jet_key]) > 0:
        jets = data[jet_key]
        axes[idx].hist(jets, bins=range(0, 15), alpha=0.7, color='red', edgecolor='black')
        axes[idx].set_xlabel('Jet Multiplicity')
        axes[idx].set_ylabel('Events')
        axes[idx].set_title(f'{name} - Jet Multiplicity')
        axes[idx].grid(True, alpha=0.3)

axes[5].axis('off')
plt.tight_layout()
plt.savefig('jet_multiplicity_all_channels.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved jet multiplicity plot")

# Plot all 5 MET definitions for each channel
for prefix, name in channels.items():
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    has_data = False
    for i in range(5):
        met_key = f'{prefix}_met{i+1}_pts'
        if met_key in data and len(data[met_key]) > 0:
            has_data = True
            met = data[met_key]
            axes[i].hist(met, bins=50, range=(0, 200), alpha=0.7, 
                        color=['blue', 'green', 'red', 'orange', 'purple'][i], 
                        edgecolor='black')
            axes[i].set_xlabel('MET [GeV]')
            axes[i].set_ylabel('Events')
            axes[i].set_title(f'MET{i+1}')
            axes[i].grid(True, alpha=0.3)
    
    axes[5].axis('off')
    
    if has_data:
        fig.suptitle(f'{name} - All MET Definitions', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'met_definitions_{prefix}.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved MET definitions plot for {name}")
    else:
        plt.close()

# Plot centrality distribution
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

for idx, (prefix, name) in enumerate(channels.items()):
    cent_key = f'{prefix}_centrality'
    if cent_key in data and len(data[cent_key]) > 0:
        centrality = data[cent_key]
        axes[idx].hist(centrality, bins=[-0.5, 0.5, 1.5], alpha=0.7, 
                      color='purple', edgecolor='black')
        axes[idx].set_xlabel('Centrality (0=Peripheral, 1=Central)')
        axes[idx].set_ylabel('Events')
        axes[idx].set_title(f'{name} - Centrality')
        axes[idx].set_xticks([0, 1])
        axes[idx].grid(True, alpha=0.3)

axes[5].axis('off')
plt.tight_layout()
plt.savefig('centrality_all_channels.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved centrality plot")

# Create summary plot comparing MET1 across all channels
fig, ax = plt.subplots(figsize=(10, 6))

colors = ['blue', 'green', 'red', 'orange', 'purple']
for idx, (prefix, name) in enumerate(channels.items()):
    met_key = f'{prefix}_met1_pts'
    if met_key in data and len(data[met_key]) > 0:
        met = data[met_key]
        ax.hist(met, bins=50, range=(0, 200), alpha=0.5, 
               label=name, color=colors[idx], histtype='step', linewidth=2)

ax.set_xlabel('MET1 [GeV]')
ax.set_ylabel('Events')
ax.set_title('MET1 Comparison Across All Channels', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('met1_comparison_all_channels.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved MET1 comparison plot")

print("\nDone! Generated plots:")
print("  - lepton_pt_all_channels.png")
print("  - lepton_eta_all_channels.png")
print("  - jet_multiplicity_all_channels.png")
print("  - met_definitions_{channel}.png (for each channel)")
print("  - centrality_all_channels.png")
print("  - met1_comparison_all_channels.png")