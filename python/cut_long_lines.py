import os

def recursive_replace(root, pattern, replace) :
 for directory, subdirs, names in os.walk(root):
  for name in names:
   if name[-8:] == "_pgi.F90":
    path = os.path.join(directory, name)
    text = open(path).read()
    if pattern in text:
     open(path,'w').write(text.replace(pattern,replace))

recursive_replace('.','#','!')
recursive_replace('.',';','\n')
