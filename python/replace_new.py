"""
  Just a simple function to replace a string inside a directory
   root : directory
   pattern : searched string
   replace "pattern" by "replace"
"""
import os
import re

def recursive_replace_new( root ) :
 for dir, subdirs, names in os.walk( root ):
  for name in names:
   path = os.path.join( dir, name )
   text = open( path ).read()
   if name[-3:] == 'F90':
    for pattern in text.split():
     if pattern.startswith("sll_f_new"):
         old_function = (pattern.split('(')[0]).split(',')[0]
         new_function = "sll_f"+old_function[9:]+"_new"
         print name + " : " + old_function + " -> " + new_function
         open( path,'w').write(text.replace(old_function, new_function))

recursive_replace_new("src")
recursive_replace_new("simulations")
recursive_replace_new("thirdparty")
