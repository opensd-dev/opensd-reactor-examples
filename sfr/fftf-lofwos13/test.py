import sys,os
import numpy as np
from scipy import interpolate
import unittest

a_path = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0,a_path)
from PINET import main
from PINET import initializer


class Test(unittest.TestCase):

    def test_case21(self): #FFTF core kinetics problem
        c_path = "model.py"
        initializer.initialize(trans_sim=False,inputpath=c_path)
        circuits,tim_sim,calcs,hslabs = main.main(trans_sim=False,verbosity=1,flag_write=False)
        data_sim = ( [getattr(calc.res[0],R) for calc in calcs if calc.res[0].__name__ == "rho_fb" 
                        for R in ['TRFL','TRCL','TRNA','TRDOP','RGEM','RBMF','RG','RC']] )
        data_sim = [x * 100000. for x in data_sim] # to pcm
        data_sim = np.asarray(data_sim).squeeze()
        data_val = [-203.6716738, 5.376634063, -15.19708581, -317.4654366, -0.302180031, 92.59463444, -119.7540521, -4.846233884] #pcm #FNX simulation trial40.proj
        np.testing.assert_allclose(data_sim, data_val, rtol=0.0, atol=3.) #pcm 

        initializer.initialize(trans_sim=True,inputpath=c_path,flag_write=False)
        circuits,tim_sim,calcs,hslabs = main.main(trans_sim=True,verbosity=1,flag_write=False)
        tim_sim=np.insert(tim_sim,0,0)
        data_sim=calcs[2].timseries[1:]
        data_sim = np.asarray(data_sim).squeeze()
        tim_val = [2.,3,5,10,15,16,18,20,25,30,35,40,50] #,75,100,125,150,200]
        tim_val = [x-2. for x in tim_val] #to overcome the time offset in FNX simulation
        data_val = [0,-5.07422E-06,-0.000113929,-0.001309529,-0.002713691,-0.002908516,-0.003123225,-0.003175427,-0.00303914,-0.002944642,-0.00292163,-0.002924204,-0.00294703] #,-0.002974931,-0.00299743,-0.003016899,-0.002935774,-0.002759539] #bar #FNX simulation trial40.proj
        data_val = [x*1.E5 for x in data_val]
        func = interpolate.interp1d(tim_sim, data_sim)
        data_sim_new = func(tim_val)
        np.testing.assert_allclose(data_sim_new, data_val, rtol=0.1, atol=10.) #pcm

if __name__ == '__main__':
    unittest.main()
    