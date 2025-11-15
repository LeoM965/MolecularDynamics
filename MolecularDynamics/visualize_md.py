import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

df = pd.read_csv("trajectory.csv")

atom_colors = {'C': '#2A4D3A', 'H': '#E8F4FD', 'O': '#E74C3C', 'N': '#3498DB', 'Be': '#2ECC71', 'Cl': '#F39C12'}
atom_sizes = {'C': 450, 'H': 250, 'O': 350, 'N': 320, 'Be': 290, 'Cl': 380}

bonds = []
try:
    with open("topology.csv", 'r') as f:
        lines = f.readlines()
    in_bonds = False
    for line in lines:
        if "Bond details" in line:
            in_bonds = True
            continue
        if in_bonds and ',' in line and not 'molecule_id' in line:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                try:
                    bonds.append((int(parts[2]), int(parts[3])))
                except:
                    pass
        if "Angle details" in line:
            break
except:
    pass

steps = sorted(df['step'].unique())

molecules = {
    'CH4': [0, 1, 2, 3, 4],      # la (10, 10, 10)
    'H2O': [5, 6, 7],            # la (13, 10, 10)
    'CO2': [8, 9, 10],           # la (16, 10, 10)
    'NH3': [11, 12, 13, 14],     # la (10, 13, 10)
    'C2H2': [15, 16, 17, 18],    # la (13, 13, 10)
    'HCl': [19, 20]              # la (16, 13, 10)
}

mol_colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#E67E22']

plt.style.use('seaborn-v0_8-whitegrid')
fig = plt.figure(figsize=(24, 16))
fig.patch.set_facecolor('#F8FAFC')

ax_3d = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=3, projection='3d')
ax_bonds = plt.subplot2grid((3, 4), (0, 2), colspan=2)
ax_angles = plt.subplot2grid((3, 4), (1, 2), colspan=2)
ax_legend = plt.subplot2grid((3, 4), (2, 2), colspan=2)

try:
    bonds_df = pd.read_csv("bonds_angles.csv")
    has_bonds_data = True
except:
    has_bonds_data = False

bond_data = {mol: [] for mol in molecules.keys()}
angle_data = {mol: [] for mol in molecules.keys()}
step_data = []

plt.ion()

for i, step in enumerate(steps):
    ax_3d.clear()
    ax_3d.set_facecolor('#FFFFFF')

    data = df[df['step'] == step]

    for _, atom in data.iterrows():
        color = atom_colors.get(atom['type'], '#95A5A6')
        size = atom_sizes.get(atom['type'], 200)
        ax_3d.scatter(atom['x'], atom['y'], atom['z'], c=color, s=size, alpha=0.95,
                      edgecolors='#34495E', linewidth=2.5, depthshade=True)
        ax_3d.text(atom['x'] + 0.35, atom['y'] + 0.35, atom['z'] + 0.35,
                   f"{atom['type']}{int(atom['atom_id'])}", fontsize=11, weight='bold', color='#2C3E50')

    for bond in bonds:
        try:
            atom1 = data[data['atom_id'] == bond[0]].iloc[0]
            atom2 = data[data['atom_id'] == bond[1]].iloc[0]
            ax_3d.plot([atom1['x'], atom2['x']], [atom1['y'], atom2['y']], [atom1['z'], atom2['z']],
                       color='#7F8C8D', linewidth=7, alpha=0.9, solid_capstyle='round')
        except:
            pass

    ax_3d.set_xlim(8, 18)
    ax_3d.set_ylim(8, 15)
    ax_3d.set_zlim(8, 12)
    ax_3d.set_xlabel('X (Å)', fontsize=15, fontweight='bold', color='#2C3E50')
    ax_3d.set_ylabel('Y (Å)', fontsize=15, fontweight='bold', color='#2C3E50')
    ax_3d.set_zlabel('Z (Å)', fontsize=15, fontweight='bold', color='#2C3E50')
    ax_3d.set_title(f'MD Simulation - Step {step} ({(i + 1) / len(steps) * 100:.0f}%)',
                    fontsize=20, weight='bold', pad=25, color='#2C3E50')
    ax_3d.view_init(elev=35, azim=60)
    ax_3d.grid(True, alpha=0.4, color='#BDC3C7')

    ax_3d.xaxis.pane.fill = False
    ax_3d.yaxis.pane.fill = False
    ax_3d.zaxis.pane.fill = False
    ax_3d.xaxis.pane.set_edgecolor('#E5E5E5')
    ax_3d.yaxis.pane.set_edgecolor('#E5E5E5')
    ax_3d.zaxis.pane.set_edgecolor('#E5E5E5')

    if has_bonds_data:
        current_data = bonds_df[bonds_df['step'] == step]
        step_data.append(step)

        for mol_id, (mol_name, atom_ids) in enumerate(molecules.items()):
            mol_bonds = current_data[(current_data['molecule_id'] == mol_id) &
                                     (current_data['measurement_type'] == 'bond')]
            mol_angles = current_data[(current_data['molecule_id'] == mol_id) &
                                      (current_data['measurement_type'] == 'angle')]

            if len(mol_bonds) > 0:
                avg_bond = mol_bonds['value'].mean()
                bond_data[mol_name].append(avg_bond)
            else:
                bond_data[mol_name].append(None)

            if len(mol_angles) > 0:
                avg_angle = mol_angles['value'].mean()
                angle_data[mol_name].append(avg_angle)
            else:
                angle_data[mol_name].append(None)

        ax_bonds.clear()
        ax_angles.clear()

        for j, (mol_name, color) in enumerate(zip(molecules.keys(), mol_colors)):
            valid_bonds = [b for b in bond_data[mol_name] if b is not None]
            valid_angles = [a for a in angle_data[mol_name] if a is not None]
            valid_steps_bonds = step_data[-len(valid_bonds):] if valid_bonds else []
            valid_steps_angles = step_data[-len(valid_angles):] if valid_angles else []

            if valid_bonds:
                ax_bonds.plot(valid_steps_bonds, valid_bonds, color=color,
                              linewidth=4.5, label=mol_name, marker='o', markersize=8,
                              alpha=0.9, markeredgecolor='white', markeredgewidth=2)

            if valid_angles:
                ax_angles.plot(valid_steps_angles, valid_angles, color=color,
                               linewidth=4.5, label=mol_name, marker='s', markersize=8,
                               alpha=0.9, markeredgecolor='white', markeredgewidth=2)

        ax_bonds.set_title('Bond Lengths Evolution', fontsize=17, weight='bold', pad=18, color='#2C3E50')
        ax_bonds.set_ylabel('Length (Å)', fontsize=14, weight='bold', color='#2C3E50')
        ax_bonds.grid(True, alpha=0.5, linestyle='--', color='#BDC3C7')
        ax_bonds.legend(fontsize=12, ncol=2, frameon=True, fancybox=True, shadow=True)
        ax_bonds.set_facecolor('#FEFEFE')

        ax_angles.set_title('Bond Angles Evolution', fontsize=17, weight='bold', pad=18, color='#2C3E50')
        ax_angles.set_xlabel('Step', fontsize=14, weight='bold', color='#2C3E50')
        ax_angles.set_ylabel('Angle (°)', fontsize=14, weight='bold', color='#2C3E50')
        ax_angles.grid(True, alpha=0.5, linestyle='--', color='#BDC3C7')
        ax_angles.legend(fontsize=12, ncol=2, frameon=True, fancybox=True, shadow=True)
        ax_angles.set_facecolor('#FEFEFE')

    ax_legend.clear()
    ax_legend.axis('off')

    ax_legend.text(0.5, 0.95, 'Atom Types', ha='center', fontsize=17, weight='bold',
                   transform=ax_legend.transAxes, color='#2C3E50')

    atom_positions = [(0.08, 0.8), (0.32, 0.8), (0.56, 0.8), (0.8, 0.8), (0.08, 0.6), (0.32, 0.6)]
    for k, (atom_type, color) in enumerate(atom_colors.items()):
        if k < len(atom_positions):
            x, y = atom_positions[k]
            ax_legend.scatter(x, y, c=color, s=350, edgecolors='#34495E', linewidth=2.5,
                              transform=ax_legend.transAxes, alpha=0.95)
            ax_legend.text(x + 0.06, y, atom_type, fontsize=15, weight='bold',
                           transform=ax_legend.transAxes, va='center', color='#2C3E50')

    ax_legend.text(0.5, 0.35, 'Molecules', ha='center', fontsize=17, weight='bold',
                   transform=ax_legend.transAxes, color='#2C3E50')

    mol_y_positions = [0.25, 0.19, 0.13, 0.07, 0.01, -0.05]
    for j, (mol_name, color) in enumerate(zip(molecules.keys(), mol_colors)):
        if j < len(mol_y_positions):
            ax_legend.plot([0.08, 0.22], [mol_y_positions[j], mol_y_positions[j]],
                           color=color, linewidth=6, transform=ax_legend.transAxes, alpha=0.9)
            ax_legend.text(0.26, mol_y_positions[j], mol_name, fontsize=13, weight='bold',
                           transform=ax_legend.transAxes, va='center', color='#2C3E50')

    plt.tight_layout()
    plt.draw()
    plt.pause(0.3)

plt.ioff()
plt.show()

print(f"\nMolecular dynamics visualization completed!")
print(f"Processed {len(steps)} simulation steps")
print(f"Tracked {len(molecules)} different molecules")