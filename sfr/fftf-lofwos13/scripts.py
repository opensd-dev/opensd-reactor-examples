import math
import sys


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

def _read_series(filename):
    points = []
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [part.strip() for part in line.strip().split(",", 1)]
            if len(parts) != 2 or not parts[0] or not parts[1]:
                continue
            t, y = parts
            points.append((float(t), float(y)))
    return points

def _interp(points, time):
    if time <= points[0][0]:
        return points[0][1]
    if time >= points[-1][0]:
        return points[-1][1]

    for i in range(len(points) - 1):
        t0, y0 = points[i]
        t1, y1 = points[i + 1]
        if t0 <= time <= t1:
            if t1 == t0:
                return y1
            frac = (time - t0) / (t1 - t0)
            return y0 + frac * (y1 - y0)

    return points[-1][1]

_msource_points = None
_decay_points = None

time = 0.0
delt = 1.0

def mdot(time, delt):
    global _msource_points
    if _msource_points is None:
        _msource_points = _read_series("msource.csv")
    return _interp(_msource_points, time)

def power_trans(time, delt):
    po = getattr(powpk, "PO", 0.986)
    hilab = [9752.7, 54413, 70593, 102010, 88518, 109640]
    hiac = [3316500, 20275000, 28279000, 41711000, 43479000, 54461000]
    hiuab = [14643, 84355, 121930, 188410, 176920, 231170]

    values = []
    for i in range(6):
        values.append(po * hilab[i])
        values.append(po * hiac[i])
        values.append(po * hiuab[i])
    values.append(po * 7113740.0)
    return values


class LumpedMass:
    TOREF = 200.0
    registry = {}

    def __init__(self, identifier, tau, node_id):
        self.identifier = identifier
        self.tau = tau
        self.node_id = node_id
        self.DTST = 0.0
        self.DTST0 = 0.0
        self.DT0 = 0.0
        self.registry[identifier] = self

    def update(self, sim_time, sim_delt):
        node = get_comp(self.node_id)
        dt = node.stemp_gues - self.TOREF
        if not hasattr(self, "_initialized"):
            self.DTST = dt
            self.DTST0 = dt
            self.DT0 = dt
            self._initialized = True
        else:
            self.DTST = (0.5 * sim_delt / self.tau * (dt + self.DT0 - self.DTST0) + self.DT0) / (
                1.0 + 0.5 * sim_delt / self.tau
            )
            self.DTST0 = self.DTST
            self.DT0 = dt


LumpedMass("MassDTRG", 4.7572736, "node1")
LumpedMass("MassDTRV", 27.314107, "node1")
LumpedMass("MassDTRC", 22.0, "node7")


def get_comp(identifier):
    if identifier in LumpedMass.registry:
        return LumpedMass.registry[identifier]
    import bindings
    return bindings.get_comp(identifier)


def update_lumped_masses():
    for mass in LumpedMass.registry.values():
        mass.update(time, delt)
    return 0.0


def rho_fb():
    # reading components
    DHSNZLAB = []
    DHSNZAC = []
    DHSNZUAB = []
    PipeNZLAB = []
    PipeNZAC = []
    PipeNZUAB = []

    for zone in range(1, 7):
        DHSNZLAB.append(get_comp(f"hslab{zone}1"))
        DHSNZAC.append(get_comp(f"hslab{zone}2"))
        DHSNZUAB.append(get_comp(f"hslab{zone}3"))
        PipeNZLAB.append(get_comp(f"pipe{zone}1"))
        PipeNZAC.append(get_comp(f"pipe{zone}2"))
        PipeNZUAB.append(get_comp(f"pipe{zone}3"))

    pipe7 = get_comp("pipe7")

    from rpdat import TOREF
    from rpdat import (
        TC1AC, TC1LAB, TC1UAB, TC2AC, TC2LAB, TC2UAB, TC3AC, TC3LAB, TC3UAB,
        TCDAC, TCDLAB, TCDUAB, TCD2AC, TCD2LAB, TCD2UAB, gemW, gemh,
    )

    rho_fb.TRFL = 0.0
    rho_fb.TRCL = 0.0
    rho_fb.TRNA = 0.0
    rho_fb.TRDOP = 0.0

    NZNR = 6
    IFR = 0
    for I in range(NZNR):
        RFL = 0.0
        RCL = 0.0
        RNA = 0.0
        RDOP = 0.0

        for K in range(10):
            TNAAC = PipeNZAC[I].faces[K].dnode.stemp_gues
            RNA = RNA + TC3AC[I][K - 1] * (TNAAC - TOREF)

            TDA3 = DHSNZAC[I].layers[0].nodes[10 + K].temp_gues
            TDA4 = (DHSNZAC[I].layers[0].nodes[10 + K].temp_gues + DHSNZAC[I].layers[1].nodes[K].temp_gues) / 2.0
            TDA6 = 0.5 * (TDA3 + TDA4)
            QDNA = TDA6 / TOREF
            if TDA6 <= 1000.0 and TDA6 > 473.0:
                RDOP = RDOP + TCDAC[I][K - 1] * math.log(QDNA)
            elif TDA6 > 1000.0 and TDA6 <= 2000.0:
                RDOP = RDOP + TCD2AC[I][K - 1] * math.log(QDNA)
            else:
                print("warning: Doppler out of range. time = " + str(time))
                sys.exit(1)

            TAC1 = DHSNZAC[I].layers[2].nodes[K].temp_gues
            TAC2 = DHSNZAC[I].layers[2].nodes[K].eface.temp_gues
            TAC = 0.5 * (TAC1 + TAC2)
            RCL = RCL + TC2AC[I][K - 1] * (TAC - TOREF)

            if IFR == 1:
                RFL = RFL + TC1AC[I][K - 1] * (TAC - TOREF)
            else:
                RFL = RFL + TC1AC[I][K - 1] * (TDA6 - TOREF)

        TLAB1 = DHSNZLAB[I].layers[2].nodes[0].temp_gues
        TLAB2 = DHSNZLAB[I].layers[2].nodes[0].eface.temp_gues
        TLAB = 0.5 * (TLAB1 + TLAB2)
        RCL = RCL + TC2LAB[I] * (TLAB - TOREF)

        TNALAB = PipeNZLAB[I].dnode.stemp_gues
        RNA = RNA + TC3LAB[I] * (TNALAB - TOREF)

        TDL3 = DHSNZLAB[I].layers[0].nodes[1].temp_gues
        TDL4 = (DHSNZLAB[I].layers[0].nodes[1].temp_gues + DHSNZLAB[I].layers[1].nodes[0].temp_gues) / 2.0
        TDL6 = 0.5 * (TDL3 + TDL4)
        QDNL = TDL6 / TOREF
        if TDL6 <= 1000.0 and TDL6 > 473.0:
            RDOP = RDOP + TCDLAB[I] * math.log(QDNL)
        elif TDL6 > 1000.0 and TDL6 <= 2000.0:
            RDOP = RDOP + TCD2LAB[I] * math.log(QDNL)
        else:
            print("warning: Doppler out of range. time = " + str(time))
            sys.exit(1)

        if IFR == 1:
            RFL = RFL + TC1LAB[I] * (TLAB - TOREF)
        else:
            RFL = RFL + TC1LAB[I] * (TDL6 - TOREF)

        TUAB1 = DHSNZUAB[I].layers[2].nodes[0].temp_gues
        TUAB2 = DHSNZUAB[I].layers[2].nodes[0].eface.temp_gues
        TUAB = 0.5 * (TUAB1 + TUAB2)
        RCL = RCL + TC2UAB[I] * (TUAB - TOREF)

        TNAUAB = PipeNZUAB[I].dnode.stemp_gues
        RNA = RNA + TC3UAB[I] * (TNAUAB - TOREF)

        TDU3 = DHSNZUAB[I].layers[0].nodes[1].temp_gues
        TDU4 = (DHSNZUAB[I].layers[0].nodes[1].temp_gues + DHSNZUAB[I].layers[1].nodes[0].temp_gues) / 2.0
        TDU6 = 0.5 * (TDU3 + TDU4)
        QDNU = TDU6 / TOREF
        if TDU6 <= 1000.0 and TDU6 > 473.0:
            RDOP = RDOP + TCDUAB[I] * math.log(QDNU)
        elif TDU6 > 1000.0 and TDU6 <= 2000.0:
            RDOP = RDOP + TCD2UAB[I] * math.log(QDNU)
        else:
            print("warning: Doppler out of range. time = " + str(time))
            sys.exit(1)

        if IFR == 1:
            RFL = RFL + TC1UAB[I] * (TUAB - TOREF)
        else:
            RFL = RFL + TC1UAB[I] * (TDU6 - TOREF)

        rho_fb.TRCL = rho_fb.TRCL + RCL
        rho_fb.TRFL = rho_fb.TRFL + RFL
        rho_fb.TRNA = rho_fb.TRNA + RNA
        rho_fb.TRDOP = rho_fb.TRDOP + RDOP

    rho_fb.RBMF = 0.0
    for I in range(NZNR):
        SUM1 = 0.0
        for K in range(10):
            TDA3 = DHSNZAC[I].layers[0].nodes[10 + K].temp_gues
            TDA4 = (DHSNZAC[I].layers[0].nodes[10 + K].temp_gues + DHSNZAC[I].layers[1].nodes[K].temp_gues) / 2.0
            TDA6 = 0.5 * (TDA3 + TDA4)
            TDL3 = DHSNZLAB[I].layers[0].nodes[1].temp_gues
            TDL4 = (DHSNZLAB[I].layers[0].nodes[1].temp_gues + DHSNZLAB[I].layers[1].nodes[0].temp_gues) / 2.0
            TDL6 = 0.5 * (TDL3 + TDL4)
            TDU3 = DHSNZUAB[I].layers[0].nodes[1].temp_gues
            TDU4 = (DHSNZUAB[I].layers[0].nodes[1].temp_gues + DHSNZUAB[I].layers[1].nodes[0].temp_gues) / 2.0
            TDU6 = 0.5 * (TDU3 + TDU4)
            SUM1 = SUM1 + TDA6 - TOREF

        SUM1 = SUM1 + TDL6 - TOREF
        SUM1 = SUM1 + TDU6 - TOREF
        rho_fb.RBMF = rho_fb.RBMF + TC1AC[I][9] * SUM1 * (-1.0)

    QR = sum([pipe.mflow for pipe in PipeNZLAB]) + pipe7.mflow
    geml = 265.0 - 539504 / (2440.13 + QR * QR / (734.08 * 3 * 734.08 * 3) * 100 * 100)
    geml = geml - 18.59
    if geml <= 25.3 or geml >= 207.35:
        print("warning: geml out of range time = " + str(time))

    rho_fb.RGEM = 0.0
    for K in range(1, 22 + 1):
        if geml <= gemh[K - 1]:
            rho_fb.RGEM = (
                gemW[K - 1] * 450.0 / 497.59601
                - (gemW[K - 1] * 450.0 / 497.59601 - gemW[K - 1 - 1] * 450.0 / 497.59601)
                / (gemh[K - 1] - gemh[K - 1 - 1])
                * (gemh[K - 1] - geml)
            )
            break

    DTC = get_comp("MassDTRC").DTST
    DTG = get_comp("MassDTRG").DTST
    DTR = get_comp("MassDTRV").DTST

    TCGR = -0.952410817E-5
    WDCR = -7.82E-5 * 1000.0
    rho_fb.RG = TCGR * (DTG - 273.15)
    ALCR = 1.60e-5
    ALRV = 1.60e-5
    EFLCR = 4.89
    EFLRV = 7.341
    DLCR = ALCR * EFLCR * (DTC - 273.15)
    DLRV = ALRV * EFLRV * (DTR - 273.15)
    rho_fb.RC = WDCR * (DLCR - DLRV)

    rho_fb.TFR = (
        rho_fb.TRFL + rho_fb.TRCL + rho_fb.TRNA + rho_fb.TRDOP
        + rho_fb.RGEM + rho_fb.RBMF + rho_fb.RG + rho_fb.RC
    )
    return rho_fb.TFR * 1.0E5


def powpk():
    global _decay_points
    if _decay_points is None:
        _decay_points = _read_series("decay.csv")

    if not hasattr(rho_fb, "TFR"):
        rho_fb.TFR = 0.0

    if not hasattr(powpk, "REXTSS"):
        powpk.REXTSS = -rho_fb.TFR
        powpk.RTOT = 0.0
    else:
        powpk.REXT = powpk.REXTSS
        powpk.RTOT = powpk.REXT + rho_fb.TFR
    DKK = powpk.RTOT

    DT = delt
    QALPHA = [0, 0, 0, 0, 0, 0]
    QGAMMA = [0, 0, 0, 0, 0, 0]
    AUXCI1 = [0, 0, 0, 0, 0, 0]
    AUXCI2 = [0, 0, 0, 0, 0, 0]

    QBETA = (8.98184E-05, 7.51112E-04, 6.56706E-04, 1.22693E-03, 5.30909E-04, 1.69559E-04)
    QLAMBD = (1.29701E-02, 3.13017E-02, 1.34437E-01, 3.41507E-01, 1.35304E+00, 3.70998E+00)
    QL = 4.780164E-07

    QB = 0.0
    for M in range(6):
        QB = QB + QBETA[M]

    for M in range(6):
        AUXCI1[M] = 2.0 + DT * QLAMBD[M]
        AUXCI2[M] = (2.0 - DT * QLAMBD[M]) / AUXCI1[M]

    if not hasattr(powpk, "Q2d"):
        powpk.Q2SS = 0.936418
        powpk.Q2d = powpk.Q2SS
        powpk.Q2 = powpk.Q2SS
        powpk.CM1d = (QBETA[0] * powpk.Q2d) / (QL * QLAMBD[0])
        powpk.CM2d = (QBETA[1] * powpk.Q2d) / (QL * QLAMBD[1])
        powpk.CM3d = (QBETA[2] * powpk.Q2d) / (QL * QLAMBD[2])
        powpk.CM4d = (QBETA[3] * powpk.Q2d) / (QL * QLAMBD[3])
        powpk.CM5d = (QBETA[4] * powpk.Q2d) / (QL * QLAMBD[4])
        powpk.CM6d = (QBETA[5] * powpk.Q2d) / (QL * QLAMBD[5])
    else:
        A1 = DT * (1.0 + DKK) / QL
        A11 = 0.0
        for M in range(6):
            QALPHA[M] = QBETA[M] * A1 / AUXCI1[M]
            A11 = A11 + QALPHA[M] * QLAMBD[M]

        QGAMMA[0] = AUXCI2[0] * powpk.CM1d + powpk.Q2d * QALPHA[0]
        QGAMMA[1] = AUXCI2[1] * powpk.CM2d + powpk.Q2d * QALPHA[1]
        QGAMMA[2] = AUXCI2[2] * powpk.CM3d + powpk.Q2d * QALPHA[2]
        QGAMMA[3] = AUXCI2[3] * powpk.CM4d + powpk.Q2d * QALPHA[3]
        QGAMMA[4] = AUXCI2[4] * powpk.CM5d + powpk.Q2d * QALPHA[4]
        QGAMMA[5] = AUXCI2[5] * powpk.CM6d + powpk.Q2d * QALPHA[5]

        A22 = 0.0
        for M in range(6):
            A22 = A22 + QGAMMA[M] * QLAMBD[M]

        A2 = (DKK - QB * (1.0 + DKK)) / QL
        powpk.Q2 = -A22 / (A2 + A11)
        powpk.Q2d = powpk.Q2
        powpk.CM1d = powpk.Q2 * QALPHA[0] + QGAMMA[0]
        powpk.CM2d = powpk.Q2 * QALPHA[1] + QGAMMA[1]
        powpk.CM3d = powpk.Q2 * QALPHA[2] + QGAMMA[2]
        powpk.CM4d = powpk.Q2 * QALPHA[3] + QGAMMA[3]
        powpk.CM5d = powpk.Q2 * QALPHA[4] + QGAMMA[4]
        powpk.CM6d = powpk.Q2 * QALPHA[5] + QGAMMA[5]

    powpk.QMOY = _interp(_decay_points, time)
    powpk.PO = powpk.Q2 + powpk.QMOY
    return powpk.PO


powpk.PO = 0.986
