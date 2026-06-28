def fun2(flow_elem,WallTemp):
    Pe = flow_elem.ther_gues.rhomass()*flow_elem.velocity*flow_elem.diameter*flow_elem.ther_gues.cpmass()/flow_elem.ther_gues.conductivity()
    Nu = 5.0 + 0.025*Pe**0.8
    h = Nu * flow_elem.ther_gues.conductivity() / flow_elem.diameter
    return h

def fun1(flow_elem,WallTemp):
    Pe = flow_elem.ther_gues.rhomass()*flow_elem.velocity*flow_elem.diameter*flow_elem.ther_gues.cpmass()/flow_elem.ther_gues.conductivity()
    P = 0.00726    #Pin Pitch
    D = 0.005842   #Pin Diameter
    Nu = 4.0 + 0.16*(P/D)**5.0 + 0.33*(P/D)**3.8*(Pe/100)**0.86
            #13.066               #7.0961
    b = Nu * flow_elem.ther_gues.conductivity() / flow_elem.diameter
    return b