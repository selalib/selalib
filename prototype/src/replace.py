import os

def recursive_replace( root, pattern, replace ) :
    for dir, subdirs, names in os.walk( root ):
        for name in names:
            if name[-3:] == 'F90' or name[-3:] == 'txt':
               path = os.path.join( dir, name )
               text = open( path ).read()
               if pattern in text:
                   print 'occurence in :' + name
                   open( path, 'w' ).write( text.replace( pattern, replace ) )

recursive_replace('.','sll_logical','sll_cartesian')
recursive_replace('.','new_logical_mesh','new_cartesian_mesh')
recursive_replace('.','get_logical_mesh','get_cartesian_mesh')
recursive_replace('.','SLL_LOGICAL_MESH','SLL_CARTESIAN_MESH')

def recursive_replace_in_cmake( root, pattern, replace ) :
    for dir, subdirs, names in os.walk( root ):
        for name in names:
            if name[-4:] == '.txt':
               path = os.path.join( dir, name )
               text = open( path ).read()
               if pattern in text:
                   print 'occurence in :' + name
                   open( path, 'w' ).write( text.replace( pattern, replace ) )

recursive_replace('.','sll_logical_meshes','sll_meshes')
