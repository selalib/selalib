#python $sll_py/map2dir.py "command" run1 run2 ...
# does "command" in dirs run1,run2 ... 
import os
import sys

if len(sys.argv)<3:
  print("Usage python map2dir.py \"command\" run1 run2")
#pwd = os.environ.get('PWD')
#if(pwd=='None'):
#  print("PWD not detected")
#  sys.exit()

cmd = sys.argv[1]
for i in sys.argv[2:]:
  os.system("cd "+i+";"+cmd)
