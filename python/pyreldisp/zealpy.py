import sys
from pyreldisp_params import params

#-------------------------------------------------------------
# Class 'zealpy'  
#-------------------------------------------------------------
class zealpy(object):
    def __init__(self):
        self.zeros         = []
        self.zealpy_module = None        
        self.params        = params()
    #end def init_


    #-------------------------------------------------------------
    # Associate the fortran module to the corresponding
    #  python package 
    #-------------------------------------------------------------
    def setup_pymodule(self,module_name):
        self.zealpy_module         = module_name
        self.math_tools_module     = self.zealpy_module.math_tools_module
        self.zeal_input_module     = self.zealpy_module.zeal_input_module
        self.split_module          = self.zealpy_module.split_module
        self.zout                  = self.zealpy_module.zout
        self.function_input_module = self.zealpy_module.function_input_module


    #-------------------------------------------------------------
    # Fill all the variables used in Fortran modules
    #-------------------------------------------------------------
    def fill_fortran_var(self):
        #--> Check the validity of the parameter required
        #-->  for function_input_module
        if (not self.params.check_validity()):
            raise pyreldisp_exception("Fille_fortran_var: Parameters not valid")


        #--> Fill variables for function_input_module
        #-->  example for one variable:
        #-->   self.function_input_module.Ti = self.params.dic['Ti']  
        for param_key in self.params.dic.keys():
            exec "self.function_input_module.%s = " \
                "self.params.dic['%s']"%(param_key,param_key)        
    #end def fill_fortran_var


    #-------------------------------------------------------------
    #  Find the zeros
    #-------------------------------------------------------------
    def get_zeros(self,xmin,xmax,ymin,ymax):
 
        # check for overflow at xmin,ymin
        fout = self.math_tools_module.check_overflow(xmin+1j*ymin)
        if abs(fout)>1.0e300:
            print 'zealpy has failed due to overflow' 
            print 'Try using a larger value of ymin'
            sys.exit()
            global zeros
        self.zeal_input_module.lv=[xmin,ymin]
        h =[float(xmax-xmin),float(ymax-ymin)]
        self.zeal_input_module.h = h
        self.zout.zeal()
        #print 'zeros ',zout.zeros
        # ZEAL failed try again with refined box
        if self.zout.fail:
            if self.zout.information == 0 :
                print 'Improper input parameters'
                sys.exit()
            if self.zout.information == 2 :
                self.split_module.trap1()
                self.split_module.trap2()
                self.split_module.trap3()
                self.split_module.trap4()
                if self.split_module.fail1 | self.split_module.fail2 | self.split_module.fail3 | self.split_module.fail4:
                    if self.split_module.fail1:
                        # we move the lower edge
                        self.zeal_input_module.lv=[xmin,ymin+0.00002]
                        h =[float(xmax-xmin),float(ymax-ymin-0.00002)]
                        self.zeal_input_module.h = h
                        self.zout.zeal()
                        if self.zout.fail:
                            if self.zout.information == 2 :
                                self.split_module.trap1()
                                if self.split_module.fail1:
                                    print 'zealpy has failed due to the zeros are too close on the lower edge'
                                    print 'Try using larger value of ymin or smaller box'
                                    sys.exit()
                                else:
                                    print 'zealpy has failed because there is a zero on the lower edge'
                                    print 'Try using larger value of ymin or smaller box'
                                    sys.exit()
                            else :
                                print 'zealpy has failed because there is a zero on the lower edge'
                                print 'Try using larger value of ymin or smaller box'
                                sys.exit()
                        else:
                            if not(self.zout.fzeros == None):
                                if (self.zout.totalnumber > 0):
                                    #   print xmin,xmax,ymin,ymax, max(abs(zout.fzeros))
                                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                        # pb with computation of zeros. Refine box
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                    else:
                                        self.zeros=self.zeros+self.zout.zeros.tolist()
                    if self.split_module.fail2:
                        # we move the right edge
                        self.zeal_input_module.lv=[xmin,ymin]
                        h =[float(xmax-xmin-0.00002),float(ymax-ymin)]
                        self.zeal_input_module.h = h
                        self.zout.zeal()
                        if self.zout.fail:
                            if self.zout.information == 2 :
                                self.split_module.trap2()
                                if self.split_module.fail2:
                                    print 'zealpy has failed due to the zeros are too close on the right edge'
                                    print 'Try using smaller box'
                                    sys.exit()
                                else:
                                    print 'zealpy has failed because there is a zero on the right edge'
                                    print 'Try using smaller box'
                                    sys.exit()
                            else :
                                print 'zealpy has failed because there is a zero on the right edge'
                                print 'Try using larger value of ymin or smaller box'
                                sys.exit()
                        else:
                            if not(self.zout.fzeros == None):
                                if (self.zout.totalnumber > 0):
                                    #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                        # pb with computation of zeros. Refine box
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                    else:
                                        self.zeros=self.zeros+self.zout.zeros.tolist()
                    if self.split_module.fail3:
                        # we move the upper edge
                        self.zeal_input_module.lv=[xmin,ymin]
                        h =[float(xmax-xmin),float(ymax-ymin-0.00002)]
                        self.zeal_input_module.h = h
                        self.zout.zeal()
                        if self.zout.fail:
                            if self.zout.information == 2 :
                                self.split_module.trap3()
                                if self.split_module.fail3:
                                    print 'zealpy has failed due to the zeros are too close on the upper edge'
                                    print 'Try using smaller box'
                                    sys.exit()
                                else:
                                    print 'zealpy has failed because there is a zero on the upper edge'
                                    print 'Try using smaller box'
                                    sys.exit()
                            else :
                                print 'zealpy has failed because there is a zero on the upper edge'
                                print 'Try using smaller box'
                                sys.exit()
                        else:
                            if not(self.zout.fzeros == None):
                                if (self.zout.totalnumber > 0):
                                    #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                        # pb with computation of zeros. Refine box
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                    else:
                                        self.zeros=self.zeros+self.zout.zeros.tolist()
                    if self.split_module.fail4:
                        # we move the right edge
                        self.zeal_input_module.lv=[xmin+0.00002,ymin]
                        h =[float(xmax-xmin-0.00002),float(ymax-ymin)]
                        self.zeal_input_module.h = h
                        self.zout.zeal()
                        if self.zout.fail:
                            if self.zout.information == 2 :
                                self.split_module.trap4()
                                if self.split_module.fail4:
                                    print 'zealpy has failed due to the zeros are too close on the left edge'
                                    print 'Try using larger value of xmin or smaller box'
                                    sys.exit()
                                else:
                                    print 'zealpy has failed because there is a zero on the left edge'
                                    print 'Try using larger value of xmin or smaller box'
                                    sys.exit()
                            else :
                                print 'zealpy has failed because there is a zero on the left edge'
                                print 'Try using larger value of xmin or smaller box'
                                sys.exit()
                        else:
                            if not(self.zout.fzeros == None):
                                if (self.zout.totalnumber > 0):
                                    #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                        # pb with computation of zeros. Refine box
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                    else:
                                        self.zeros=self.zeros+self.zout.zeros.tolist()
                else:
                    print 'Zealpy has failed because the procedure for the calculation of the total number of zeros has failed'
                    sys.exit()
                
                
            if self.zout.information == 3:
                # we move the lower edge
                self.zeal_input_module.lv=[xmin,ymin+0.0002]
                h =[float(xmax-xmin),float(ymax-ymin-0.0002)]
                self.zeal_input_module.h = h
                self.zout.zeal()
                if self.zout.fail:
                    # we move the right edge
                    self.zeal_input_module.lv=[xmin,ymin]
                    h =[float(xmax-xmin-0.0002),float(ymax-ymin)]
                    self.zeal_input_module.h = h
                    self.zout.zeal()
                    if self.zout.fail:
                        # we move the upper edge
                        self.zeal_input_module.lv=[xmin,ymin]
                        h =[float(xmax-xmin),float(ymax-ymin-0.0002)]
                        self.zeal_input_module.h = h
                        self.zout.zeal()
                        if self.zout.fail:
                            # we move the left edge
                            self.zeal_input_module.lv=[xmin+0.0002,ymin]
                            h =[float(xmax-xmin-0.0002),float(ymax-ymin)]
                            self.zeal_input_module.h = h
                            self.zout.zeal()
                            if self.zout.fail:
                                print 'Zealpy has failed because the procedure for the isolation of the zeros has failed'
                                sys.exit()
                            else:
                                if not(self.zout.fzeros == None):
                                    if (self.zout.totalnumber > 0):
                                        #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                        if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                            xmin=xmin+0.0002
                                            # pb with computation of zeros. Refine box
                                            self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                            self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                            self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                            self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                        else:
                                            self.zeros=self.zeros+self.zout.zeros.tolist()
                        else:
                            if not(self.zout.fzeros == None):
                                if (self.zout.totalnumber > 0):
                                    #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                        # pb with computation of zeros. Refine box
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                    else:
                                        self.zeros=self.zeros+self.zout.zeros.tolist()
                
                    else:
                        if not(self.zout.fzeros == None):
                            if (self.zout.totalnumber > 0):
                                #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                                if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                    # pb with computation of zeros. Refine box
                                    self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                    self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                    self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                    self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                                else:
                                    self.zeros=self.zeros+self.zout.zeros.tolist()
                else:
                    if not(self.zout.fzeros == None):
                        if (self.zout.totalnumber > 0):
                            #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                            if (max(abs(self.zout.fzeros)) > 1.0e-10):
                                ymin=ymin+0.0002
                                # pb with computation of zeros. Refine box
                                self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                                self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                                self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                                self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                            else:
                                self.zeros=self.zeros+self.zout.zeros.tolist()
                
            if self.zout.information == 4:
                print 'Zealpy has failed because the procedure for the computation of the zeros has failed'
                sys.exit()
        else:
            if not(self.zout.fzeros == None):
                if (self.zout.totalnumber > 0):
                    #   print xmin,xmax,ymin,ymax, max(abs(self.zout.fzeros))
                    if (max(abs(self.zout.fzeros)) > 1.0e-10):
                        # pb with computation of zeros. Refine box
                        self.get_zeros(xmin,xmin+0.5*h[0],ymin,ymin+0.5*h[1])
                        self.get_zeros(xmin+0.5*h[0],xmax,ymin,ymin+0.5*h[1])
                        self.get_zeros(xmin,xmin+0.5*h[0],ymin+0.5*h[1],ymax)
                        self.get_zeros(xmin+0.5*h[0],xmax,ymin+0.5*h[1],ymax)
                    else:
                        self.zeros=self.zeros+self.zout.zeros.tolist()
                                    

