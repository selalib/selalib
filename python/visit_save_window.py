#visit -nowin -cli -s "$sll_py/visit_save_window.py test_visit $PWD  filnam1 ... filnamN"
#filnam1.xmf ... filnamN.xmf are used
#filnam1_0000.png ... filnamN_0000.png are generated
#if they do not already exist
#or other format depending on the SaveWindowAttributes
#of your configuration of visit
#otherwise filnam1_0001.png is generated and so on
#this depends on the preferences of visit
#if family is not set on then filnam1_.png is created
#instead of filnam1_0001.png
#in the user specific file test_visit.py should be
#x1_min,x1_max,x2_min,x2_max,fmin,fmax	
#seems to work as expected only for one filename

import sys
import math
import importlib
print(sys.argv)
cnt=0
sys.path.insert(0, sys.argv[2]) #os.getcwd())
themodule = importlib.import_module(sys.argv[1])

att=visit.View2DAttributes()
att.fullFrameActivationMode=att.On
x1_min = themodule.x1_min
x1_max = themodule.x1_max
x2_min = themodule.x2_min
x2_max = themodule.x2_max
att.windowCoords = (x1_min, x1_max, x2_min, x2_max)
att.viewportCoords = (0.2, 0.95, 0.15, 0.95)
visit.SetView2D(att)
fmin = themodule.fmin
fmax = themodule.fmax

for filname in sys.argv[3:]: 
  cnt = cnt+1
  if(cnt==1):
    visit.OpenDatabase(filname+".xmf")
    visit.AddPlot("Pseudocolor", "values")
    att=visit.PseudocolorAttributes()
    att.min = fmin
    att.max = fmax
    visit.SetPlotOptions(att)
    visit.DrawPlots()
  else:
    visit.OpenDatabase(filname+".xmf")
    visit.ReplaceDatabase(filname+".xmf")  
  s = visit.SaveWindowAttributes()
  s.fileName = filname+"_"
  visit.SetSaveWindowAttributes(s)
  visit.SaveWindow()
