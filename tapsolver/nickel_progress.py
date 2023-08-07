import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

data = np.random.rand(14, 8) * 40
data = np.zeros((14, 8))
print(data)
data[0][:] = 15
data[1][:] = 15
data[2][7] = 15
data[3][7] = 15
data[4][7] = 15
# create discrete colormap
cmap = colors.ListedColormap(['red', 'yellow', 'green'])
bounds = [0,10,20,30]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(4/1.2,7/1.2))
ax.imshow(data, cmap=cmap, norm=norm)

surface = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,'']
species = ['CH3CH3*','CH3CH2*','CH3CH*','CH2CH2*','CH2CH*','H*','H2*','*','']

# draw gridlines
ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
ax.set_xticks(np.arange(-0.5, 8, 1));
ax.set_yticks(np.arange(-0.5, 14, 1));

ax.set_ylabel('Surface Configuration')
ax.set_yticklabels(surface)
#ax.set_xlabel('Surface Species')
ax.set_xticklabels(species)
plt.xticks(rotation=90)
fig.tight_layout()
plt.show()

data = np.random.rand(14, 8) * 40
data = np.zeros((14, 8))
print(data)
# create discrete colormap
cmap = colors.ListedColormap(['red', 'yellow', 'green'])
bounds = [0,10,20,30]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(4/1.2,7/1.2))
ax.imshow(data, cmap=cmap, norm=norm)

surface = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,'']
species = ['CH3CH3*','CH3CH2*','CH3CH*','CH2CH2*','CH2CH*','H*','H2*','*','']

# draw gridlines
ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
ax.set_xticks(np.arange(-0.5, 8, 1));
ax.set_yticks(np.arange(-0.5, 14, 1));

ax.set_ylabel('Surface Configuration')
ax.set_yticklabels(surface)
#ax.set_xlabel('Surface Species')
ax.set_xticklabels(species)
plt.xticks(rotation=90)
fig.tight_layout()
plt.show()

data = np.random.rand(14, 5) * 40
data = np.zeros((14, 5))
print(data)
# create discrete colormap
cmap = colors.ListedColormap(['red', 'yellow', 'green'])
bounds = [0,10,20,30]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(4/1.4,7/1.2))
ax.imshow(data, cmap=cmap, norm=norm)

surface = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,'']
species = ['2','3','5','7','8','']

# draw gridlines
ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
ax.set_xticks(np.arange(-0.5, 5, 1));
ax.set_yticks(np.arange(-0.5, 14, 1));

ax.set_ylabel('Surface Considered')
ax.set_yticklabels(surface)
ax.set_xlabel('Transition State')
ax.set_xticklabels(species)
fig.tight_layout()

plt.show()

data = np.random.rand(14, 5) * 40
data = np.zeros((14, 5))
print(data)
# create discrete colormap
cmap = colors.ListedColormap(['red', 'yellow', 'green'])
bounds = [0,10,20,30]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(4/1.4,7/1.2))
ax.imshow(data, cmap=cmap, norm=norm)

surface = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,'']
species = ['2','3','5','7','8','']

# draw gridlines
ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
ax.set_xticks(np.arange(-0.5, 5, 1));
ax.set_yticks(np.arange(-0.5, 14, 1));

ax.set_ylabel('Surface Considered')
ax.set_yticklabels(surface)
ax.set_xlabel('Transition State')
ax.set_xticklabels(species)
fig.tight_layout()

plt.show()