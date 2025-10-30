import pymeshlab
import sys
ms = pymeshlab.MeshSet()

# input file name of .ply point cloud
ms.load_new_mesh(sys.argv[1])
ms.compute_normal_for_point_clouds()
ms.generate_surface_reconstruction_screened_poisson()

m = ms.current_mesh()
#print(len(m.face_matrix()))

# close small holes for a watertight mesh
if len(m.face_matrix()) > 0:
    ms.meshing_close_holes(maxholesize=100)

out_dict = ms.get_geometric_measures()
#print(out_dict['mesh_volume'])

# record volume, 0 if mesh not watertight
if 'mesh_volume' in out_dict:
    v = out_dict['mesh_volume']
else:
    v = 0

#print(v)

