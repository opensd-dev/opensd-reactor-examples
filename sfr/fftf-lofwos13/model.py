import opensd

circuit1=opensd.Circuit(identifier="circuit1")
circuit1.assign_fluid("Na23",fltype="incompressible",fllib="User")

circuit1.add_node("node1",elevation=0)
circuit1.add_node("node7",elevation=2.4519)

bc1=circuit1.add_BC("bc1","node1",'T',598.899)
bc2=circuit1.add_BC("bc2","node1",'P',100.E5)
global bc3
bc3=circuit1.add_BC("bc3","node7",'msource',-2191.7992)
#bc3=circuit1.add_BC("bc3","node7",'P',376464.8406894)

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

FA=0.00433
D=0.003238
E=[0.0000,0.9144,0.9347,1.0795,2.1719]#elevation  of nodes in channel
PDC=[9.7689,10.23113,9.77031,9.94184,11.84798,25.21945]#pressure drop coefficient
UN3FA=0.0052128
UN3D=0.07648207121336499
DF = [D,D,D,D,D,UN3D]
NINCF = [1,10,1,1,1,1]
PF=[1,6,9,16,18,30]#no. of parallel flows in a channel
L=[0.0203,0.9144,0.0203,0.1448,1.0924,0.28]#length of pipe
FALAB=[0.0808,0.4851,0.7276,1.2936,1.4553,2.4254]
FALABI=[0.0703,0.4218,0.6327,1.1248,1.2654,2.1091]
GIFALAB=[0.066,0.399,0.598,1.063,1.196,1.993]
FAAC=[3.642,21.850,32.776,58.268,65.551,109.252]
FAACI=[3.167,19,28.5,50.668,65.551,95.002]
GIFAAC=[3.079,18.477,27.715,49.271,55.430,92.383]
FAUAB=[0.081,0.485,0.728,1.294,1.455,2.425]
FAUABI=[0.070,0.422,0.633,1.125,1.265,2.109]
GIFAUAB=[0.066,0.399,0.598,1.063,1.196,1.993]
FAUN1=[0.577,3.460,5.189,9.226,10.379,17.298]
FAUN2=[4.351,26.104,39.156,69.610,78.311,130.519]
FAUN3=[0.076336,0.458016,0.687024,1.221376,1.374048,2.29008]
global HIUAB
HIUAB=[14643,84355,121930,188410,176920,231170]
global HIAC
HIAC=[3316500,20275000,28279000,41711000,43479000,54461000]
global HILAB
HILAB=[9752.7,54413,70593,102010,88518,109640]
AFF=[[0.077776,0.093363,0.108730,0.118943,0.123089,0.120939,0.112795,0.099371,0.081658,0.063335],
     [0.078353,0.093682,0.108895,0.119016,0.123074,0.120839,0.112640,0.099237,0.081493,0.062772],
     [0.080316,0.094493,0.109211,0.119091,0.122948,0.120487,0.112116,0.098895,0.081216,0.061226],
     [0.081924,0.095950,0.110592,0.120322,0.123891,0.120782,0.111400,0.097186,0.078868,0.059086],
     [0.084594,0.098532,0.113178,0.122793,0.126394,0.122085,0.109379,0.092788,0.074502,0.055754],
     [0.088890,0.101690,0.115791,0.124753,0.126694,0.120845,0.106322,0.088871,0.071330,0.054815]]
global p
p=[]

for i in range(6):
    for j in range(5):
        a=str(i+1)+str(j+1)
        circuit1.add_node("node"+a,elevation=E[j])

    for k in range(6):
        b=str(i+1)+str(k+1)
        c=str(i+1)+str(k)
        if k==0:    circuit1.add_pipe("pipe"+b,DF[k],L[k],"node1", "node"+b,0.03,1,NINCF[k],cfarea=FA,   npar=PF[i],Kforward=PDC[i])
        elif k==1:  circuit1.add_pipe("pipe"+b,DF[k],L[k],"node"+c,"node"+b,0.03,1,NINCF[k],cfarea=FA,   npar=PF[i])
        elif 1<k<5: circuit1.add_pipe("pipe"+b,DF[k],L[k],"node"+c,"node"+b,0.03,1,NINCF[k],cfarea=FA,   npar=PF[i])
        else:       circuit1.add_pipe("pipe"+b,DF[k],L[k],"node"+c,"node7", 0.03,1,NINCF[k],cfarea=UN3FA,npar=PF[i])

    for k in range(6):
        d=str(i+1)+str(k+1)
        opensd.SNode("snode"+d)
        if k==0:
            g=opensd.HSlab("hslab"+d,"snode"+d,"hflux",0.0,"pipe"+d,"pipe",[fun1],GIFALAB[i],nlayers=3)
            p.append(g.add_layer(0.0024,0.0203,2,GIFALAB[i],'MOX23','User',heat_input=HILAB[i]))
            g.add_layer(1.4E-4,0.0203,2,FALABI[i],'gap24','User')
            g.add_layer(3.81E-4,0.0203,2,FALAB[i],'SS23','User')
        elif k==1:
            g=opensd.HSlab("hslab"+d,"snode"+d,"hflux",0.0,"pipe"+d,"pipe",[fun1],GIFAAC[i],nlayers=3)
            p.append(g.add_layer(0.00247,0.9144,2,GIFAAC[i],'MOX23','User',10,heat_input=HIAC[i],AFF=AFF[i]))
            g.add_layer(7.E-5,0.9144,2,FAACI[i],'gap23','User')
            g.add_layer(3.81E-4,0.9144,2,FAAC[i],'SS23','User')
        elif k==2:
            g=opensd.HSlab("hslab"+d,"snode"+d,"hflux",0.0,"pipe"+d,"pipe",[fun1],GIFAUAB[i],nlayers=3)
            p.append(g.add_layer(0.0024,0.0203,2,GIFAUAB[i],'MOX23','User',heat_input=HIUAB[i]))
            g.add_layer(1.4E-4,0.0203,2,FAUABI[i],'gap24','User')
            g.add_layer(3.81E-4,0.0203,2,FAUAB[i],'SS23','User')
        elif k==3:
            g=opensd.HSlab("hslab"+d,"pipe"+d,"pipe",[fun1],"snode"+d,"hflux",0.0,FAUN1[i],nlayers=1)
            g.add_layer(0.0018497,0.1448,2,FAUN1[i],'SS23','User')
        elif k==4:
            g=opensd.HSlab("hslab"+d,"pipe"+d,"pipe",[fun1],"snode"+d,"hflux",0.0,FAUN2[i],nlayers=1)
            g.add_layer(8.625E-4,1.0924,2,FAUN2[i],'SS23','User')
        else:
            g=opensd.HSlab("hslab"+d,"pipe"+d,"pipe",[fun2],"snode"+d,"hflux",0.0,FAUN3[i],nlayers=1)
            g.add_layer(0.027,0.28,2,FAUN3[i],'SS23','User')

global pipe7
pipe7=circuit1.add_pipe("pipe7",0.003238,2.4722,"node1","node7",0.03,1,1,cfarea=0.12982,Kforward=147.6114,heat_input=7113740.)

opensd.LumpedMass("MassDTRG",4.7572736,"node1")
opensd.LumpedMass("MassDTRV",27.314107,"node1")
opensd.LumpedMass("MassDTRC",22.0,"node7")

#_____________________________________________________________________________________________________________________________

def csv_reader(filename):
    import pandas as pd
    from scipy.interpolate import interp1d
    df = pd.read_csv(filename,header=None)
    t = df.iloc[:,0]
    y = df.iloc[:,1]
    y = interp1d(t,y)
    return y

global series1
series1=csv_reader("decay.csv")
def power_trans(time,delt):
    if time == 0:
        powpk.PO = 0.986
    for m in range(6):
        for n in range(3):
            if n==0:
                p[m*3+n].heat_input = powpk.PO*HILAB[m]
            elif n==1:
                p[m*3+n].heat_input = powpk.PO*HIAC[m]
            else:
                p[m*3+n].heat_input = powpk.PO*HIUAB[m]

    pipe7.heat_input = powpk.PO * 7113740.
# action_setup.Action(None,None,power_trans)

# times  = [0.0, 20.0, 50.0]
# values = [-753.6, 0.0, 0.0]
# tdist = opensd.Tabular(times, values)

# a1 = opensd.Action("ramp_bc2", "bc2", "bval", tdist)



global series2
series2=csv_reader("msource.csv")
def mdot(time,delt):
    y=series2(time)
    return y
# action_setup.Action("bc3","bval",mdot)


settings = opensd.Settings()

settings.conv_crit_ht = 1.E-7
settings.conv_crit_flow = 1.E-7
settings.conv_crit_temp_trans = 1.E-7

settings.no_main_iter = 500
settings.verbosity = 1
settings.temp_solve = False
settings.run_mode = "transient"
settings.tim_slot = [[1., 51.0]]
settings.flag_write = True
settings.export_to_xml()

global rho_fb
def rho_fb():

    # reading components
    from PINET.project import get_comp
    DHSNZLAB  = []
    DHSNZAC   = []
    DHSNZUAB  = []
    PipeNZLAB = []
    PipeNZAC  = []
    PipeNZUAB = []
    
    # radial zone1
    DHSNZLAB .append(get_comp("hslab11"))
    DHSNZAC  .append(get_comp("hslab12"))
    DHSNZUAB .append(get_comp("hslab13"))
    PipeNZLAB.append(get_comp("pipe11") )
    PipeNZAC .append(get_comp("pipe12") )
    PipeNZUAB.append(get_comp("pipe13") )

    DHSNZLAB .append(get_comp("hslab21"))
    DHSNZAC  .append(get_comp("hslab22"))
    DHSNZUAB .append(get_comp("hslab23"))
    PipeNZLAB.append(get_comp ("pipe21") )
    PipeNZAC .append(get_comp ("pipe22") )
    PipeNZUAB.append(get_comp ("pipe23") )
    DHSNZLAB .append(get_comp("hslab31"))
    DHSNZAC  .append(get_comp("hslab32"))
    DHSNZUAB .append(get_comp("hslab33"))
    PipeNZLAB.append(get_comp ("pipe31") )
    PipeNZAC .append(get_comp ("pipe32") )
    PipeNZUAB.append(get_comp ("pipe33") )
    DHSNZLAB .append(get_comp("hslab41"))
    DHSNZAC  .append(get_comp("hslab42"))
    DHSNZUAB .append(get_comp("hslab43"))
    PipeNZLAB.append(get_comp ("pipe41") )
    PipeNZAC .append(get_comp ("pipe42") )
    PipeNZUAB.append(get_comp ("pipe43") )
    DHSNZLAB .append(get_comp("hslab51"))
    DHSNZAC  .append(get_comp("hslab52"))
    DHSNZUAB .append(get_comp("hslab53"))
    PipeNZLAB.append(get_comp ("pipe51") )
    PipeNZAC .append(get_comp ("pipe52") )
    PipeNZUAB.append(get_comp ("pipe53") )
    DHSNZLAB .append(get_comp("hslab61"))
    DHSNZAC  .append(get_comp("hslab62"))
    DHSNZUAB .append(get_comp("hslab63"))
    PipeNZLAB.append(get_comp ("pipe61") )
    PipeNZAC .append(get_comp ("pipe62") )
    PipeNZUAB.append(get_comp ("pipe63") )

    from rpdat import TOREF
    
    from rpdat import TC1AC, TC1LAB, TC1UAB, TC2AC, TC2LAB, TC2UAB, TC3AC, TC3LAB, TC3UAB, TCDAC, TCDLAB, TCDUAB, TCD2AC, TCD2LAB, TCD2UAB, gemW, gemh

    # calculate rho_fb.TRFL, rho_fb.TRCL, rho_fb.TRNA, rho_fb.TRDOP
    rho_fb.TRFL  = 0.0
    rho_fb.TRCL  = 0.0
    rho_fb.TRNA  = 0.0
    rho_fb.TRDOP = 0.0

    NZNR = 6 #no. of radial neutronic zones
    for I in range(NZNR):
        RFL  = 0.0
        RCL  = 0.0
        RNA  = 0.0
        RDOP = 0.0

        for K in range(10): #pending automate
            # fuel region 
            TNAAC = PipeNZAC[I].faces[K].dnode.stemp_gues
            RNA = RNA + TC3AC[I][K-1]*(TNAAC-TOREF)
            
            TDA3 = DHSNZAC[I].layers[0].nodes[10+K].temp_gues
            TDA4 = (DHSNZAC[I].layers[0].nodes[10+K].temp_gues+DHSNZAC[I].layers[1].nodes[K].temp_gues)/2.
            TDA6=0.5*(TDA3+TDA4)
            QDNA=TDA6/TOREF
            if (TDA6<=1000.0 and TDA6 > 473.0):#pending automate
                RDOP = RDOP + TCDAC[I][K-1]*math.log(QDNA)
            elif (TDA6 > 1000.0 and TDA6 <= 2000.0):
                RDOP = RDOP + TCD2AC[I][K-1]*math.log(QDNA)
            else:
                print("warning: Doppler out of range. Time = " + time)
                sys.exit()

            TAC1 = DHSNZAC[I].layers[2].nodes[K].temp_gues
            TAC2 = DHSNZAC[I].layers[2].nodes[K].eface.temp_gues
            TAC=0.5*(TAC1+TAC2)
            RCL = RCL + TC2AC[I][K-1]*(TAC-TOREF)

            IFR = 0
            if (IFR==1): # fuel stuck to clad
                RFL = RFL + TC1AC[I][K-1]*(TAC-TOREF)
            elif (IFR==0): # fuel free to expand
                RFL = RFL + TC1AC[I][K-1]*(TDA6-TOREF)

        # bottom axial blanket region
        TLAB1 = DHSNZLAB[I].layers[2].nodes[0].temp_gues
        TLAB2 = DHSNZLAB[I].layers[2].nodes[0].eface.temp_gues
        TLAB=0.5*(TLAB1+TLAB2)
        RCL = RCL + TC2LAB[I]*(TLAB-TOREF)

        TNALAB = PipeNZLAB[I].dnode.stemp_gues
        RNA = RNA + TC3LAB[I]*(TNALAB-TOREF)

        TDL3 = DHSNZLAB[I].layers[0].nodes[1].temp_gues
        TDL4 = (DHSNZLAB[I].layers[0].nodes[1].temp_gues+DHSNZLAB[I].layers[1].nodes[0].temp_gues)/2.
        TDL6=0.5*(TDL3+TDL4)
        QDNL=TDL6/TOREF
        if (TDL6<=1000.0 and TDL6 > 473.0):
            RDOP = RDOP + TCDLAB[I]*math.log(QDNL)
        elif (TDL6 > 1000.0 and TDL6 <= 2000.0):
            RDOP = RDOP + TCD2LAB[I]*math.log(QDNL)
        else:
          print("warning: Doppler out of range. time = " + Time)
          sys.exit()

        if (IFR==1): #fuel stuck to clad
            RFL = RFL + TC1LAB[I]*(TLAB-TOREF)
        elif (IFR==0): #fuel free to expand
            RFL = RFL + TC1LAB[I]*(TDL6-TOREF)

        # upper axial blanket
        TUAB1 = DHSNZUAB[I].layers[2].nodes[0].temp_gues
        TUAB2 = DHSNZUAB[I].layers[2].nodes[0].eface.temp_gues
        TUAB = 0.5*(TUAB1+TUAB2)
        RCL = RCL + TC2UAB[I]*(TUAB-TOREF)

        TNAUAB = PipeNZUAB[I].dnode.stemp_gues
        RNA = RNA + TC3UAB[I]*(TNAUAB-TOREF)

        TDU3 = DHSNZUAB[I].layers[0].nodes[1].temp_gues
        TDU4 = (DHSNZUAB[I].layers[0].nodes[1].temp_gues+DHSNZUAB[I].layers[1].nodes[0].temp_gues)/2.
        TDU6=0.5*(TDU3+TDU4)
        QDNU=TDU6/TOREF
        if (TDU6<=1000.0 and TDU6 > 473.0):
            RDOP = RDOP + TCDUAB[I]*math.log(QDNU)
        elif (TDU6 > 1000.0 and TDU6 <= 2000.0):
            RDOP = RDOP + TCD2UAB[I]*math.log(QDNU)
        else:
            print("warning: Doppler out of range. time = " + Time)
            sys.exit()

        if (IFR==1): #fuel stuck to clad
            RFL = RFL + TC1UAB[I]*(TUAB-TOREF)
        elif (IFR==0): #fuel free to expand
            RFL = RFL + TC1UAB[I]*(TDU6-TOREF)

        rho_fb.TRCL = rho_fb.TRCL + RCL
        rho_fb.TRFL = rho_fb.TRFL + RFL
        rho_fb.TRNA = rho_fb.TRNA + RNA
        rho_fb.TRDOP = rho_fb.TRDOP + RDOP
    
    # calculate rho_fb.RBMF
    rho_fb.RBMF = 0.0
    for I in range(NZNR):
        SUM1 = 0.0
        for K in range(10): # fuel region
            if (IFR==1):
                TAC1 = DHSNZAC[I].layers[2].nodes[K].temp_gues
                TAC2 = DHSNZAC[I].layers[2].nodes[K].eface.temp_gues
                TLAB1 = DHSNZLAB[I].layers[2].nodes[0].temp_gues
                TLAB2 = DHSNZLAB[I].layers[2].nodes[0].eface.temp_gues
                TUAB1 = DHSNZUAB[I].layers[2].nodes[0].temp_gues
                TUAB2 = DHSNZUAB[I].layers[2].nodes[0].eface.temp_gues
                TAC=0.5*(TAC1+TAC2)
                SUM1=SUM1+TAC-TOREF
                TLAB=0.5*(TLAB1+TLAB2)
                TUAB=0.5*(TUAB1+TUAB2)
            else:
                TDA3 = DHSNZAC[I].layers[0].nodes[10+K].temp_gues
                TDA4 = (DHSNZAC[I].layers[0].nodes[10+K].temp_gues+DHSNZAC[I].layers[1].nodes[K].temp_gues)/2.
                TDA6 = 0.5*(TDA3+TDA4)
                TDL3 = DHSNZLAB[I].layers[0].nodes[1].temp_gues
                TDL4 = (DHSNZLAB[I].layers[0].nodes[1].temp_gues+DHSNZLAB[I].layers[1].nodes[0].temp_gues)/2.
                TDL6 = 0.5*(TDL3+TDL4)
                TDU3 = DHSNZUAB[I].layers[0].nodes[1].temp_gues
                TDU4 = (DHSNZUAB[I].layers[0].nodes[1].temp_gues+DHSNZUAB[I].layers[1].nodes[0].temp_gues)/2.
                TDU6 = 0.5*(TDU3+TDU4)
                SUM1 = SUM1+TDA6-TOREF

        if (IFR==1): #upper and lower axial blanket region
            SUM1=SUM1+TLAB-TOREF
            SUM1=SUM1+TUAB-TOREF
        else:
            SUM1=SUM1+TDL6-TOREF
            SUM1=SUM1+TDU6-TOREF

        rho_fb.RBMF = rho_fb.RBMF + TC1AC[I][9]*SUM1*(-1.0);

    # calculate RGEM
    QR = sum([pipe.mflow for pipe in PipeNZLAB]) + pipe7.mflow
    geml = 265.0-539504/(2440.13+QR*QR/(734.08*3*734.08*3)*100*100)
    geml = geml-18.59
    if (geml <= 25.3 or geml >= 207.35):
        print("warning: geml out of range time = " + Time)

    for K in range(1,22+1):
        if (geml <= gemh[K-1]):
            rho_fb.RGEM = gemW[K-1]*450.0/497.59601-(gemW[K-1]*450.0/497.59601-gemW[K-1-1]*450.0/497.59601)/(gemh[K-1]-gemh[K-1-1])*(gemh[K-1]-geml)
            break

    # calculate TCGR, TCCR
    # QSDTRC = get_comp("QSDTRC")
    mass = get_comp("MassDTRC")
    DTC = mass.DTST
    mass = get_comp("MassDTRG")
    DTG = mass.DTST
    mass = get_comp("MassDTRV")
    DTR = mass.DTST

    TCGR = -0.952410817E-5
    WDCR = -7.82E-5*1000.0
    rho_fb.RG = TCGR*(DTG-273.15)
    ALCR = 1.60e-5
    ALRV = 1.60e-5
    EFLCR = 4.89
    EFLRV = 7.341
    DLCR = ALCR*EFLCR*(DTC-273.15)
    DLRV = ALRV*EFLRV*(DTR-273.15)
    rho_fb.RC = WDCR*(DLCR-DLRV)

    rho_fb.TFR = rho_fb.TRFL+rho_fb.TRCL+rho_fb.TRNA+rho_fb.TRDOP+rho_fb.RGEM+rho_fb.RBMF+rho_fb.RG+rho_fb.RC
    return rho_fb.TFR*1.E5

global powpk
def powpk():
    from PINET.scheduler import time,delt
    
    # inputs from simulation
    if (time==0):
        powpk.REXTSS = -rho_fb.TFR
        powpk.RTOT = 0
    else:
        powpk.REXT = powpk.REXTSS #+REXTCR.Value+REXTSC.Value;
        powpk.RTOT = powpk.REXT + rho_fb.TFR
    DKK = powpk.RTOT

    # Q2 calculation
    DT = delt

    # local variable initialization
    QALPHA = [0,0,0,0,0,0]
    QGAMMA = [0,0,0,0,0,0]
    AUXCI1 = [0,0,0,0,0,0]
    AUXCI2 = [0,0,0,0,0,0]

    # user inputs
    QBETA = (8.98184E-05,7.51112E-04,6.56706E-04,1.22693E-03,5.30909E-04,1.69559E-04)
    QLAMBD = (1.29701E-02,3.13017E-02,1.34437E-01,3.41507E-01,1.35304E+00,3.70998E+00)
    QL = 4.780164E-07


    #transient preparatory calculation
    QB = 0.0
    for M in range(6):
        QB = QB + QBETA[M]
    
    for M in range(6):
        AUXCI1[M] =  2.0+DT*QLAMBD[M]
        AUXCI2[M] = (2.0-DT*QLAMBD[M])/AUXCI1[M]
    AUXCI1 = tuple(i for i in AUXCI1)
    AUXCI2 = tuple(i for i in AUXCI2)

    if (time==0): # SS calculation
        powpk.Q2SS = 0.936418

        powpk.Q2d = powpk.Q2SS
        powpk.Q2  = powpk.Q2SS

        powpk.CM1d = (QBETA[0]*powpk.Q2d)/(QL*QLAMBD[0])
        powpk.CM2d = (QBETA[1]*powpk.Q2d)/(QL*QLAMBD[1])
        powpk.CM3d = (QBETA[2]*powpk.Q2d)/(QL*QLAMBD[2])
        powpk.CM4d = (QBETA[3]*powpk.Q2d)/(QL*QLAMBD[3])
        powpk.CM5d = (QBETA[4]*powpk.Q2d)/(QL*QLAMBD[4])
        powpk.CM6d = (QBETA[5]*powpk.Q2d)/(QL*QLAMBD[5])
        powpk.CM1  = powpk.CM1d
        powpk.CM2  = powpk.CM2d
        powpk.CM3  = powpk.CM3d
        powpk.CM4  = powpk.CM4d
        powpk.CM5  = powpk.CM5d
        powpk.CM6  = powpk.CM6d


    else: #transient calculation
        A1 = DT*(1.0+DKK)/QL
        A11 = 0.0
        for M in range(6):
            QALPHA[M] = QBETA[M]*A1/AUXCI1[M]
            A11       = A11 + QALPHA[M]*QLAMBD[M]

        QGAMMA[0] = AUXCI2[0]*powpk.CM1d + powpk.Q2d*QALPHA[0]
        QGAMMA[1] = AUXCI2[1]*powpk.CM2d + powpk.Q2d*QALPHA[1]
        QGAMMA[2] = AUXCI2[2]*powpk.CM3d + powpk.Q2d*QALPHA[2]
        QGAMMA[3] = AUXCI2[3]*powpk.CM4d + powpk.Q2d*QALPHA[3]
        QGAMMA[4] = AUXCI2[4]*powpk.CM5d + powpk.Q2d*QALPHA[4]
        QGAMMA[5] = AUXCI2[5]*powpk.CM6d + powpk.Q2d*QALPHA[5]
        
        A22 = 0.0
        for M in range(6):
            A22 = A22 + QGAMMA[M]*QLAMBD[M]

        A2 = (DKK - QB*(1.0+DKK))/QL
        powpk.Q2 = -A22/(A2+A11)
        powpk.Q2d = powpk.Q2

        powpk.CM1  = powpk.Q2*QALPHA[0] + QGAMMA[0]
        powpk.CM2  = powpk.Q2*QALPHA[1] + QGAMMA[1]
        powpk.CM3  = powpk.Q2*QALPHA[2] + QGAMMA[2]
        powpk.CM4  = powpk.Q2*QALPHA[3] + QGAMMA[3]
        powpk.CM5  = powpk.Q2*QALPHA[4] + QGAMMA[4]
        powpk.CM6  = powpk.Q2*QALPHA[5] + QGAMMA[5]
        powpk.CM1d = powpk.CM1
        powpk.CM2d = powpk.CM2
        powpk.CM3d = powpk.CM3
        powpk.CM4d = powpk.CM4
        powpk.CM5d = powpk.CM5
        powpk.CM6d = powpk.CM6

    powpk.QMOY = series1(time) #0.076656
    powpk.PO   = powpk.Q2 + powpk.QMOY
    return powpk.PO

def fun3(*comps):
    return comps[0].DTST

geometry = opensd.Geometry([circuit1])
geometry.export_to_xml()

# action_setup.Action(None,None,rho_fb)
# post.Calculate(rho_fb)
# post.Calculate(powpk)
# post.Calculate(lambda: powpk.RTOT*1.E5)
# post.Calculate(lambda: powpk.Q2)
# post.Calculate(lambda: rho_fb.TRFL)
# post.Calculate(lambda: rho_fb.TRDOP)
# post.Calculate(lambda: rho_fb.RBMF)
# post.Calculate(lambda: rho_fb.TRCL)
# post.Calculate(lambda: rho_fb.TRNA)
# post.Calculate(lambda: rho_fb.RGEM)
# post.Calculate(lambda: rho_fb.RG)
# post.Calculate(lambda: rho_fb.RC)

# post.Calculate(fun3,"MassDTRC")
# post.Calculate(fun3,"MassDTRG")
# post.Calculate(fun3,"MassDTRV")
