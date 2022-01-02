def sphere_volume(radius):
	return (4/3)*3.14159*(radius**3)

def sphere_area(radius):
	return (4)*3.14159*(radius**2)

cluster_diameter = 90
lattice_edge_length = 8.173
cr_atoms_per_lattice = 16
o_atoms_per_lattice = 24
cr_radius = 2.6/2
o_radius = 0.74/2
latices_in_cluster = sphere_volume(90)/sphere_volume(lattice_edge_length)
cluster_volume = sphere_volume(cluster_diameter)
cluster_area = sphere_area(cluster_diameter)

cr_al_surf_ratio = 0.075
cr6_cr3_surf_ratio = 0.25
cr6_cr3_bulk_ratio = 0.075
#cr6_concentration = 
chromium_atoms_per = 3
unit_cell_area = lattice_edge_length**2
atoms_per_cell = chromium_atoms_per/unit_cell_area
#total_area_m_g = 

#print(cr_atoms_per_lattice*sphere_volume(90)/sphere_volume(lattice_edge_length))

