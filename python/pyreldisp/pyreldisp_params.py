import sys
import numpy as np
import copy

class pyreldisp_exception(Exception):
    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

#---------------------------------------------------------------
# Class params
#---------------------------------------------------------------
class params(object):
    def __init__(self,required_params=None,user_values=None):
        self.required_params   = copy.copy(required_params)
        self.non_assigned_int  = -10000000
        self.non_assigned_real = -10000000.
        self.dic = {}

        if (required_params != None):
            for param in required_params:
                if required_params[param] == int:
                    self.dic[param] = self.non_assigned_int
                elif required_params[param] == float:
                    self.dic[param] = self.non_assigned_real
                elif required_params[param] == np.ndarray :
                    self.dic[param] = np.array([])
            #end for
        #end if

        if (user_values != None):
            for param in user_values: 
                self.set_value(param,user_values[param])
        #end if
    #end def __init__


    #-----------------------------------------------------
    # Set the value 'val' to a key 'key'
    #-----------------------------------------------------
    def set_value(self,key,val):
        if not self.is_valid_entry(key, val):  
            raise pyreldisp_exception('Error set_value: key=%s, val=%s'%(key,str(val))) 
        self.dic[key] = val

        
    #-----------------------------------------------------
    # Get the value of the key 'key'
    #-----------------------------------------------------
    def get_value(self,key):
        if not self.is_valid(key):
            raise pyreldisp_exception('Error get_value: key=%s'%(key)) 
        return self.dic[key]


    #-----------------------------------------------------
    # Check if all the required values are filled
    #-----------------------------------------------------
    def check_validity(self):
        ret = True
        for param_key in self.dic:
            ret = ret and self.is_valid(param_key)

        return ret
        #end for
    #end def check_validity

    #-----------------------------------------------------
    # Check if all the required values are of the good type and are initialized
    #-----------------------------------------------------
    def is_valid(self,param):
        if param not in self.dic:
            print 'Param %s not in parameter set, parameter available: %s'%(param,self.dic.keys()) 
            return False
        return self.is_valid_entry(param,self.dic[param])

    def is_valid_entry(self,param,val):
        if param not in self.required_params:
            print 'Param %s not in parameter set, parameter available: %s'%(param,self.required_params.keys()) 
            return False
        required_type = self.required_params[param]
        if not isinstance(val,required_type):
            print 'Param %s is of type %s but %s required'%(param,type(val),required_type)
            return False
        if required_type == int and val == self.non_assigned_int:
            print 'Param %s not initialized or wrong type'%param
            return False
        if required_type == float and val == self.non_assigned_real:
            print 'Param %s not initialized or wrong type'%param
            return False
        if required_type == np.ndarray and np.size(val) == 0:
            print 'Param %s not initialized or wrong type'%param
            return False

        return True
#end class params


#**************************************************************
#  UNIT TESTS
#**************************************************************
import unittest
class test_params(unittest.TestCase):
  
  def are_dic_equal(self,a,b):
    for k in a:
      if k not in b:
        print 'not equal because of k ',(k)
        return False
      if a[k] != b[k]:
        print 'not equal because of val[%s]= %s and %s ',(k,a[k],b[k])
        return False

    for k in b:
      if k not in a:
        print 'not equal because of k ',(k)
        return False
      if a[k] != b[k]:
        print 'not equal because of val[%s]= %s and %s ',(k,a[k],b[k])
        return False

    return True  
      
  def test_constructor0(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    expected = {'nut': -10000000.0, 'bean': np.array([]), 'apple': -10000000}
    p = params(required_params)
    res = self.are_dic_equal(expected,p.dic)
    self.assertTrue(res, "Expected=%s but got %s"%(expected,p.dic))

  def test_constructor1(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    user_params = {'apple':1}
    expected = {'nut': -10000000.0, 'bean': np.array([]), 'apple': 1}
    p = params(required_params, user_params)
    res = self.are_dic_equal(expected,p.dic)
    self.assertTrue(res, "Expected=%s but got %s"%(expected,p.dic))

  def is_valid_entry(self, key, val, expected):
    required_params = {'apple':np.int32, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)

    res = p.is_valid_entry(key,val)
    self.assertTrue(res==expected, "Expected=%s but got %s"%(expected,res))

  def test_is_valid_entry0(self):
    a = np.int32(1)
    self.is_valid_entry('apple', a, True)

  def test_is_valid_entry1(self):
    self.is_valid_entry('bean', np.ones(10), True)

  def test_is_valid_entry2(self):
    self.is_valid_entry('nut', 2.0, True)

  #check if numpy.float64 are recognized as float
  def test_is_valid_entry7(self):
    a = np.ones(10)
    self.is_valid_entry('nut', a[0], True)

  #check if numpy.int32 are recognized as int
  def test_is_valid_entry8(self):
    a = np.ones(10, dtype=np.int32)
    self.is_valid_entry('apple', a[0], True)

  #check wrong type 
  def test_is_valid_entry3(self):
    self.is_valid_entry('apple', 'wrong_type', False)

  def test_is_valid_entry4(self):
    self.is_valid_entry('bean', 'wrong_type', False)

  def test_is_valid_entry5(self):
    self.is_valid_entry('nut', 'wrong_type', False)
     
  #check if numpy.float32 are not recognized as float
  def test_is_valid_entry9(self):
    a = np.ones(10, dtype=np.float32)
    self.is_valid_entry('nut', a[0], False)

  #check if numpy.int64 are recognized as int
  def test_is_valid_entry10(self):
    a = np.ones(10, dtype=np.int64)
    self.is_valid_entry('apple', a[0], False)

  #check wrong key
  def test_is_valid_entry6(self):
    self.is_valid_entry('wrong_key', 'wrong_type', False)

  def is_valid(self, key, set_user, expected):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    user_params = {'apple':1, 'nut':2.0, 'bean':np.ones(10)}
    if set_user:
      p = params(required_params, user_params)
    else:
      p = params(required_params)

    res = p.is_valid(key)
    self.assertTrue(res==expected, "Expected=%s but got %s"%(expected,res))

  def test_is_valid0(self):
    self.is_valid('apple', True, True)

  def test_is_valid1(self):
    self.is_valid('bean', True, True)

  def test_is_valid2(self):
    self.is_valid('nut', True, True)
  #Check initialization
  def test_is_valid3(self):
    self.is_valid('apple', False, False)

  def test_is_valid4(self):
    self.is_valid('bean', False, False)

  def test_is_valid5(self):
    self.is_valid('nut', False, False)

  def test_check_validity0(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    user_params = {'apple':1, 'nut':2.0, 'bean': np.ones(10)}
    p = params(required_params, user_params)
    res = p.check_validity()
    expected = True
    self.assertTrue(res==expected, "Expected=%s but got %s"%(expected,res))

  def test_check_validity1(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    user_params = {'apple':1, 'bean': np.ones(10)}
    p = params(required_params, user_params)
    res = p.check_validity()
    expected = False
    self.assertTrue(res==expected, "Expected=%s but got %s"%(expected,res))

  def test_set_value0(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    p.set_value('apple',1)
    expected = {'nut': -10000000.0, 'bean': np.array([]), 'apple': 1}
    res = self.are_dic_equal(expected,p.dic)
    self.assertTrue(res, "Expected=%s but got %s"%(expected,p.dic))

  def test_set_value1(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    self.assertRaises(pyreldisp_exception, p.set_value, 'apple',1.0)

  def test_set_value2(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    self.assertRaises(pyreldisp_exception, p.set_value, 'wrong_key',1)

  def test_set_value3(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    self.assertRaises(pyreldisp_exception, p.set_value, 'apple',-10000000)


  def test_get_value0(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    p.set_value('apple',1)
    res = p.get_value('apple')
    expected = 1
    self.assertTrue(res==expected, "Expected=%s but got %s"%(expected,p.dic))


  def test_get_value1(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    p = params(required_params)
    self.assertRaises(pyreldisp_exception, p.get_value, 'apple')

  def test_set_value2(self):
    required_params = {'apple':int, 'nut':float, 'bean':np.ndarray}
    user_params = {'apple':1, 'nut':2.0, 'bean': np.ones(10)}
    p = params(required_params, user_params)
    self.assertRaises(pyreldisp_exception, p.get_value, 'wrong_key')

  
def runTest(method_to_test):
  if len(method_to_test) == 0:
    suite = unittest.TestLoader().loadTestsFromTestCase(test_params)
  else:
    suite = unittest.TestSuite()
    for i in method_to_test:
      suite.addTest(test_params(i))

  unittest.TextTestRunner(verbosity=2).run(suite)
  

  
if __name__ == "__main__":
  runTest(sys.argv[1:])
  
