from solidframe import solidframe

class solid(solidframe):
    def __init__(self,opt='all'): #properties are classified as thermophysical and mechanical properties. Depending on the component where it is assigned, the corresponding properties are alone loaded into the object.
        self.opt = opt
        if opt == 'ther' or opt == 'all':
            self._rhomass = 10600.
            self._cpmass = 430.*0.5
            self._conductivity = 3.
        if opt == 'mech' or opt == 'all':
            pass
            # self._youngs_modulus=1.E11
            # self._poissons_ratio=0.26
        