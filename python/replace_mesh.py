import os

#if os.path.isdir("logical_meshes"):
#   os.system("git mv logical_meshes meshes")
#if os.path.isdir("logical_meshes_multipatch"):
#   os.system("git mv logical_meshes_multipatch meshes_multipatch")

def recursive_replace( root, pattern, replace) :
    for dir, subdirs, names in os.walk( root ):
        for name in names:
            if name[-3:] == "F90":
               path = os.path.join( dir, name )
               text = open( path ).read()
               if pattern in text:
                   print 'occurence in :' + name
                   open( path, 'w' ).write( text.replace( pattern, replace ) )

#recursive_replace('.','sll_logical','sll_cartesian')
#recursive_replace('.','new_logical_mesh','new_cart#esian_mesh')
#recursive_replace('.','get_logical_mesh','get_cart#esian_mesh')
#recursive_replace('.','SLL_LOGICAL_MESH','SLL_CART#ESIAN_MESH')
#recursive_replace('.','logical_mesh2d','cartesian_#mesh2d')
#recursive_replace('.','logical_mesh4d','cartesian_#mesh4d')
#recursive_replace('.','logical_mesh4D','cartesian_#mesh4d')
#recursive_replace('.','_logical_mesh_1d','_cartesi#an_mesh_1d')
#recursive_replace('.','_logical_mesh_2d','_cartesi#an_mesh_2d')
#recursive_replace('.','_logical_mesh_3d','_cartesi#an_mesh_3d')
#recursive_replace('.','_logical_mesh_4d','_cartesi#an_mesh_4d')
#recursive_replace('.','new_logical_mesh_2d','new_c#artesian_mesh_2d')
#recursive_replace('.','_logical_mesh1d','_cartesia#n_mesh1d')
#recursive_replace('.','unit_test_logical','unit_test')
recursive_replace('.','logical_meshes_multipatch','meshes_multipatch')
