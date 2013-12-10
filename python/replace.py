import os

def recursive_replace( root, pattern, replace ) :
    for dir, subdirs, names in os.walk( root ):
        for name in names:
            path = os.path.join( dir, name )
            text = open( path ).read()
            if pattern in text:
                open( path, 'w' ).write( text.replace( pattern, replace ) )
