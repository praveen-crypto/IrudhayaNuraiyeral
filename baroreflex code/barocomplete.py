import numpy as np

Aav, Amv, Apv, Atv, Gpw = np.zeros(5)
Elaa, Elab, Elva, Elvb, Eraa, Erab, Erva, Ervb, Epua, Epuc, Epuv, Epwa, Epwc, Epwv = np.zeros(14)
yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv = np.zeros(10)
Ra, Raa, Rav, Rca, Rda, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwa, Rpwc, Rpwv, Rtv, Rv, Rvc, bav, bmv, bpv, btv = np.zeros(20)
Spua, Spuc, Spuv, Spwa, Spwc, Spwv = np.zeros(6)
Zpua, Zpuc, Zpuv, Zpwa, Zpwc, Zpwv = np.zeros(6)
Caor, Cart, Ccap, Cven, Cvca=np.zeros(5)
yaor, yart, ycap, yven, yvca=np.zeros(5)
Raor, Rart, Rcap, Rv, Rvc=np.zeros(5)
Saor, Sart, Scap, Sven, Svca =np.zeros(5)
P_0d =np.zeros(shape=(2,101))
dvq =np.zeros(2)
dv, v, dq, q=np.zeros(4)
elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv, Sra, Srv, ppp, ppc, pit, qco, FL, FR1, STR=np.zeros(19)
Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1, Rpv1, Rtv1, bav1, bmv1, bpv1, btv1=np.zeros(16)
Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1, ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2=np.zeros(20)
n_val, m_cvst, m_cvrg, n_vrg=np.zeros(4)
timestep, Tduration, ddt, tee, tcr, tac, tar, t, odic=np.zeros(9)
Rven,Rvca=np.zeros(2)
aoa,nstep,aao,aod=np.zeros(4)
gainElv,gainR,gainErv,gainTv,gainTs=np.zeros(5)


def barocomplete():

    def Integrated_ode():
        global Aav, Amv, Apv, Atv, Gpw #Valve Parameters
        global  Elaa, Elab, Elva, Elvb, Eraa, Erab, Erva, Ervb, Epua, Epuc, Epuv, Epwa, Epwc, Epwv #E_cardiopul
        global  yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv #yL_cardiopul
        global  Ra, Raa, Rav, Rca, Rda, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwa, Rpwc, Rpwv, Rtv, Rv, Rvc, bav, bmv, bpv, btv #R_cardiopul
        global  Spua, Spuc, Spuv, Spwa, Spwc, Spwv #S_cardiopul
        global  Zpua, Zpuc, Zpuv, Zpwa, Zpwc, Zpwv  #Z_cardiopul
        global  Caor, Cart, Ccap, Cven, Cvca #C_peripheral
        global  yaor, yart, ycap, yven, yvca #yL_peripheral
        global  Raor, Rart, Rcap, Rven, Rvca #R_peripheral
        global  Saor, Sart, Scap, Sven, Svca #S_peripheral
        global  dvq ,P_0d  ,nstep, gainR, gainElv, gainErv, gainTs, gainTv #sdvsdqdvdq
        global  dv, v, dq, q, aoa, aao, aod  #dvdq_cardiopul 
        global  elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv, Sra, Srv ,ppp, ppc, pit, qco, FL, FR1, STR  #cardiac_parameter
        global  Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1, Rpv1, Rtv1, bav1, bmv1, bpv1, btv1 #R_cardiopulc
        global Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1, ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2
        global n_val, m_cvst, m_cvrg, n_vrg
        global  Tduration, ddt, tcr, tee, tac, tar, t #timestep

        dvq[0, 0] = q[0, 14] - q[0, 0]  # Venous volume  dvq(1)= q(15) - q(1);

        P_0d[0, 0] = v[0, 0] / Cven + Sven * dv[0, 0]  # Venous  Pressure

        dvq[0, 1] = (v[0, 0] / Cven + Sven * dv[0, 0] - Rven * q[0, 0] - v[0, 1] / Cvca - Svca * dv[
                0, 1]) / yven  # Venous flow

        dvq[0, 2] = q[0, 0] - q[0, 1]  # VC volume

        P_0d[0, 1] = v[0, 1] / Cvca + Svca * dv[0, 1]  # VC Pressure

        dvq[0, 3] = (v[0, 1] / Cvca - era * v[0, 2] - Rvca * q[0, 1] + Svca * dv[0, 1] - Sra * dv[
               0, 2] - ppc - pit) / yvca  # VC Flow

        qco = 0.0

        dvq[0, 4] = q[0, 1] + qco - q[0, 2]  # RA volume

        P_0d[0, 3] = era * v[0, 2] + Sra * dv[0, 2] + ppc + pit  # RA pressure

        dvq[0, 5] = (era * v[0, 2] - prv - Rtv * q[0, 2] - btv * q[0, 2] * np.abs(q[0, 2]) + Sra * dv[0, 2] - Srv * dv[
                0, 3]) / ytv  # TV flow

        dvq[0, 6] = q[0, 2] - q[0, 3]  # RV volume

        P_0d[0, 5] = prv + Srv * dv[0, 3] + ppc + pit  # RV pressure

        dvq[0, 7] = (prv - Epua * Zpua - Rpv * q[0, 3] - bpv * q[0, 3] * (np.abs(q[0, 3])) + Srv * dv[0, 3] - Spua * dv[
                0, 4] + ppc) / ypv  # PV flow

        dvq[0, 8] = q[0, 3] - q[0, 4] - q[0, 7]  # Pulmonary Artery volume

        P_0d[0, 7] = Epua * Zpua + Spua * dv[0, 4] + pit  # Pulmonary Artery Pressure

        dvq[0, 9] = (Epua * Zpua - Epuc * Zpuc - Rpua * q[0, 4] + Spua * dv[0, 4] - Spuc * dv[
                0, 5]) / ypua  # Pulmonary Artery Flow

        dvq[0, 10] = q[0, 4] - q[0, 5]  # Pulmonary Capillary volume

        P_0d[0, 9] = Epuc * Zpuc + Spuc * dv[0, 5] + pit  # Pulmonary Capillary Pressure

        dvq[0, 11] = (Epuc * Zpuc - Epuv * Zpuv - Rpuc * q[0, 5] + Spuc * dv[0, 5] - Spuv * dv[
                0, 6]) / ypuc  # Pulmonary Capillary Flow

        dvq[0, 12] = q[0, 5] - q[0, 6]  # Pulmonary Vein volume

        P_0d[0, 11] = Epuv * Zpuv + Spuv * dv[0, 6] + pit  # Pulmonary Vein Pressure

        dvq[0, 13] = (Epuv * Zpuv - ela * v[0, 9] - Rpuv * q[0, 6] + Spuv * dv[0, 6] - Sla * dv[
                0, 9] - ppc) / ypuv  # Pulmonary Vein flow

        if Gpw > 0:
            dvq[0, 14] = (Epua * Zpua - Epwc * Zpwc - q[0, 7] / gpw + Spua * dv[0, 4] - Spwc * dv[
                    0, 7]) / ypwa  # Pulmonary Wedge ##Artery flow
        else:
            dvq[0, 14] = 0

        dvq[0, 15] = q[0, 7] - q[0, 8]  # Pulmonary Wedge capillary volume

        P_0d[0, 14] = Epwc * Zpwc + Spwc * dv[0, 7] + pit  # Pulmonary Wedge capillary Pressure

        dvq[0, 16] = (Epwc * Zpwc - Epwv * Zpwv - Rpwc * q[0, 8] + Spwc * dv[0, 7] - Spwv * dv[
            0, 8]) / ypwc  # Pulmonary Wedge capillary Flow

        dvq[0, 17] = q[0, 8] - q[0, 9]  # Pulmonary Wedge vein Volume

        P_0d[0, 16] = Epwv * Zpwv + Spwv * dv[0, 8] + pit  # Pulmonary Wedge vein Pressure

        dvq[0, 18] = (Epwv * Zpwv - ela * v[0, 9] - Rpwv * q[0, 9] + Spwv * dv[0, 8] - Sla * dv[
                0, 9] - ppc) / ypwv  # Pulmonary Wedge vein Flow

        dvq[0, 19] = q[0, 6] + q[0, 9] - q[0, 10]  # LA volume

        P_0d[0, 18] = ela * v[0, 9] + Sla * dv[0, 9] + ppc + pit  # LA Pressure

        dvq[0, 20] = (ela * v[0, 9] - plv - Rmv * q[0, 10] - bmv * q[0, 10] * np.abs(q[0, 10]) + Sla * dv[0, 9] - Slv *
                          dv[0, 10]) / ymv  # Mitral flow

        dvq[0, 21] = q[0, 10] - q[0, 11]  # LV volume

        P_0d[0, 20] = plv + Slv * dv[0, 10] + ppc + pit  # LV Pressure

        dvq[0, 22] = (plv - v[0, 11] / Caor - Rav * q[0, 11] - (bav * q[0, 11]) * np.abs(q[0, 11]) + Slv * dv[
                0, 10] - Saor * dv[0, 11] + ppc + pit) / yav  # Aortic Flow

        # Adjust the state of cardiac valve
        if Aav == 0.0 and q[0, 11] <= 0.000000001:
            dvq[0, 22] = 0.0

        if np.abs(tcr - Tduration) <= ddt * 0.5 or tcr < 0.1 or (Amv == 0.0 and q[0, 10] < 0.000000001):
            dvq[0, 20] = 0.0

        if Apv == 0.0 and q[0, 3] <= 0.00000001:
            dvq[0, 7] = 0.0

        if np.abs(tcr - Tduration) <= ddt * 0.5 or tcr < 0.1 or (Atv == 0.0 and q[0, 2] <= 0.000000001):
            dvq[0, 5] = 0.0

        dvq[0, 23] = q[0, 11] - q[0, 12]  # Aorta volume

        P_0d[0, 22] = v[0, 11] / Caor + Saor * dv[0, 11]  # Aorta Presuure

        dvq[0, 24] = (v[0, 11] / Caor + Saor * dv[0, 11] - q[0, 12] * Raor - v[0, 12] / Cart - Sart * dv[
                0, 12]) / yaor  # Aorta Flow

        dvq[0, 25] = q[0, 12] - q[0, 13]  # Artery Volume

        P_0d[0, 24] = v[0, 12] / Cart + Sart * dv[0, 12]  # Artery Pressure

        dvq[0, 26] = (v[0, 12] / Cart + Sart * dv[0, 12] - q[0, 13] * Rart - v[0, 13] / Ccap - Scap * dv[
            0, 13]) / yart  # Artery Flow

        dvq[0, 27] = q[0, 13] - q[0, 14]  # Capillarya Volume

        P_0d[0, 26] = v[0, 13] / Ccap + Scap * dv[0, 13]  # Capillary Pressure

        dvq[0, 28] = (v[0, 13] / Ccap + Scap * dv[0, 13] - q[0, 14] * Rcap - v[0, 0] / Cven - Sven * dv[
            0, 0]) / ycap  # Capillary Flow

        dvq[0, 29] = (aao+(6.37*aod)-q[0, 15])/2.076
        dvq[0, 30] = (-q[0, 16]+gainR)/6
        dvq[0, 31] = (-q[0, 17]+gainElv)/8
        dvq[0, 32] = (-q[0, 18]+gainErv)/8
        dvq[0, 33] = (-q[0, 19]+gainTs)/2
        dvq[0, 34] = (-q[0, 20]+gainTv)/1.5
       
    def Ecal(EEE, ZZZ, vol):
        EcalR = EEE * np.exp(vol / ZZZ)
        return EcalR

    def Lvecal():
        global Elva, Elvb
        global elv, FL
        global tee, tcr
        tcal = tcr
        if tcal <= tee:
            elv = FL * Elva * 0.5 * (1.0 - np.cos(3.1415926 * tcal / tee)) + Elvb / FL
        else:
            if tcal <= 1.5 * tee:
                elv = FL * Elva * 0.5 * (1.0 + np.cos(3.1415926 * (tcal - tee) / (0.5 * tee))) + Elvb / FL
            else:
                elv = Elvb / FL

    def Laecal():
        global tcr, ela, Elaa, tac, tar, Tduration, Elab
        tcal = tcr
        teec = tar - tac
        teer = teec
        tap = tar + teer - Tduration
        if (tcal >= 0.0 and tcal <= tap):
            ela = Elaa * 0.5 * (1.0 + np.cos(3.1415926 * (tcal + Tduration - tar) / teer)) + Elab
        
        if (tcal > tap and tcal <= tac):
            ela = Elab
            
        if (tcal > tac and tcal <= tar):
            ela = Elaa * 0.5 * (1.0 - np.cos(3.1415926 * (tcal - tac) / teec)) + Elab
        
        if (tcal > tar and tcal <= Tduration):
            ela = Elaa * 0.5 * (1.0 + np.cos(3.1415926 * (tcal - tar) / teer)) + Elab

    def Rvecal():  
        global tcr, FR1, Erva, tee, Ervb
        global erv
        tcal = tcr
        if tcal <= tee:
            erv = FR1 * Erva * 0.5 * (1.0 - np.cos(3.1415926 * tcal / tee)) + Ervb / FR1
        else:
            if tcal <= 1.5 * tee:
                erv = FR1 * Erva * 0.5 * (1.0 + np.cos(2.0 * 3.1415926 * (tcal - tee) / tee)) + Ervb / FR1
            else:
                erv = Ervb / FR1

    def Raecal():
        global Eraa, Erab
        global tar, tac, tcr
        global era, Tduration
        teec = tar - tac
        teer = teec
        tcal = tcr
        tap = tar + teer - Tduration
        if 0 <= tcal <= tap:
            era = Eraa * 0.5 * (1.0 + np.cos(3.1415926 * (tcal + Tduration - tar) / teer)) + Erab
        if tap < tcal <= tac:
            era = Erab
        if tcal > tac and tcal <= tar:
            era = Eraa * 0.5 * (1.0 - np.cos(3.1415926 * (tcal - tac) / teec)) + Erab
        if tcal > tar and tcal <= Tduration:
            era = Eraa * 0.5 * (1.0 + np.cos(3.1415926 * (tcal - tar) / teer)) + Erab

    def AAav():
        global Caor, dv, v
        global plv, Slv, ppc
        intee = plv + Slv * dv[0, 10] + ppc - v[0, 11] / Caor
        if intee > 0.0:
            AAav = 4.0
        else:
            AAav = 0.0
        return AAav

    def AAmv():
        global v, ela, plv
        intee = ela * v[0, 9] - plv
        if intee > 0.0:
            AAmv = 4.0
        else:
            AAmv = 0.0
        return AAmv

    def AApv():
        global Epua, Zpua, prv
        intee = prv - Epua * Zpua
        if intee > 0.0:
            AApv = 4.0
        else:
            AApv = 0.0
        return AApv

    def AAtv():
        global era, v, prv
        intee = era * v[0, 2] - prv
        if intee > 0.0:
            AAtv = 4.0
        else:
            AAtv = 0.0
        return AAtv

    def rkbaro(subresultcr):
        global Aav, Amv, Apv, Atv
        global P_0d, dvq
        global dv, v, dq, q 
        global timestep, Tduration, ddt
        
        # Declaration of variables
        dfl = np.zeros(shape=(2, 35))
        dq = np.zeros(shape=(2, 30))
        subrukuk = np.zeros(shape=(4, 35))
        inter = np.zeros(shape=(2, 35))
        newpara = np.zeros(shape=(2, 35))
        
        # ---------------------------------------------------------
        # call Integrated_ode 
        Integrated_ode()
        
        #---------------------------------------------------------
        dfl[0, 0:35]=dvq[0, 0:35]
        
        dv[0, 0:7] = dfl[0, 0:13:2]
        dv[0, 7:14] = dfl[0, 15:28:2]

        #dq[0, 0:7]   = dfl[0, 1:14:2]
        #dq[0, 7:15]  = dfl[0, 14:29:2]
        #dq[0, 15:21] = dfl[0, 29:35:1]

        for nrk in range(4):
            subrukuk[nrk, 0:35] = np.multiply(ddt, dfl[0, 0:35])
            if nrk == 0:               
                Aav = AAav()
                Amv = AAmv()
                Apv = AApv()
                Atv = AAtv() 
            if nrk<4:
                inter[0, 0:35] = np.multiply(0.5, (subrukuk[nrk, 0:35]))
            else:
                inter[0, 0:35]=subrukuk[nrk, 0:35]
        inter[0, 0:35]= np.multiply(0.5, (subrukuk[nrk, 0:35]))
        newpara[0, 0:35]= np.add(subresultcr[0, 0:35], inter[0, 0:35])
        
        q[0, 0:7]=newpara[0, 1:14:2]
        q[0, 7:15]=newpara[0, 14:29:2]
        q[0, 15:21]=newpara[0, 29:35:1]
        
        v[0, 0:7]=newpara[0, 0:13:2]
        v[0, 7:14]=newpara[0, 15:28:2]
        return subrukuk


    global  Aav, Amv, Apv, Atv, Gpw # valve,
    global  Elaa, Elab, Elva, Elvb, Eraa, Erab, Erva, Ervb, Epua, Epuc, Epuv, Epwa, Epwc, Epwv#E_cardiopul,
    global  yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv
    global  Ra, Raa, Rav, Rca, Rda, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwa, Rpwc, Rpwv, Rtv, Rv, Rvc, bav, bmv, bpv, btv
    global  Spua, Spuc, Spuv, Spwa, Spwc, Spwv
    global  Zpua, Zpuc, Zpuv, Zpwa, Zpwc, Zpwv
    global  Caor, Cart, Ccap, Cven, Cvca
    global  yaor, yart, ycap, yven, yvca
    global  Raor, Rart, Rcap, Rven, Rvca
    global  Saor, Sart, Scap, Sven, Svca
    global  dvq, P_0d, aoa, nstep, aao, aod
    global  dv, v, dq, q, gainElv, gainR, gainErv, gainTv, gainTs#dvdq_cardiopul,
    global  elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv ,Sra,Srv ,ppp ,ppc, pit, Pit, qco, FL, FR1, STR
    global  Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1 ,Rpv1, Rtv1 ,bav1, bmv1, bpv1, btv1
    global  Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1 ,ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2
    global  n_val, m_cvst, m_cvrg, n_vrg
    global  timestep, Tduration, ddt, tcr, tee, tac, tar, t
    global  odic
    
    v = np.zeros(shape=(2, 15))
    dv = np.zeros(shape=(2, 21))
    
    q = np.zeros(shape=(2, 36))
    dvq = np.zeros(shape=(2, 36))
    
    result = np.zeros(shape=(2, 101))
    P_0d = np.zeros(shape=(2, 101))
    diffv = np.zeros(shape=(2, 20))
    resultcr = np.zeros(shape=(2, 101))
    rukuk = np.zeros(shape=(4, 110))

    MyResult = np.zeros(shape=(100000, 39))
    MyResult1 = np.zeros(shape=(100000, 39))
    MyResult2 = np.zeros(shape=(100000, 39))
    Ervbaro, Rartbaro, Elvbaro, Tbaro=np.zeros(4)

    odic_new = np.zeros(35)
    
    odic=np.array([1068.2371, 52.4983,	181.7233,	-41.7618,	65.0625,	0,	122.6637,	
                   0,	67.0272,	-0.3118,	135.11,	-2.1737,	198.7568,	-64.1791,	
                   0,	2.7983,	0.1357,	2.7932,	1.1042,	68.8587,	0,	121.2539,	0,	
                   67.3641,	41.8262,	22.0472,	56.6627,	1.8539,	57.6473, 120,
                   0.8, 1.4, 0.5, -0.13, 0])

    Pit = -2.5
    pit = Pit
    jj = 0
    kk = 0

    # ======================================
    #Elva =  2.87     #!Peak-systolic elastance of left ventricle
    Elvb = 0.06      #!Basic diastolic elastance of left ventricle
    Elaa = 0.07      #!Peak-systolic elastance of left atrium
    Elab = 0.075     #!Basic diastolic elastance of left atrium
    #Erva = 0.52      #!Peak-systolic elastance of right ventricle
    Ervb = 0.043     #!Basic diastolic elastance of right ventricle
    Eraa = 0.055     #!Peak-systolic elastance of right atrium
    Erab = 0.06      #!Basic diastolic elastance of right atrium
    
    Vmax = 900.0    #Reference volume of Frank-Starling law
    Es = 45.9       #!Effective septal elastance
    Vpc0 = 380.0    #!Reference total pericardial and cardiac volume
    Vpe = 30.0      #!Pericardial volume of heart
    Vcon = 40.0     #!Volume constant
    Sva0 = 0.0005   #!Coefficient of cardiac viscoelasticity

    # ! Cardiac valve parameters
    #!(aortic valve(AV),mitral valve(MV), tricuspid valve(TV),pulmonary valve(PV))
    bav = 0.000025 #!Bernoulli's resistance of AV
    bmv = 0.000016 #!Bernoulli's resistance of MV
    btv = 0.000016#!Bernoulli's resistance of TV
    bpv = 0.000025 #!Bernoulli's resistance of PV
    Rav = 0.005    #!Viscous resistance of AV
    Rmv = 0.005    #!Viscous resistance of MV
    Rtv = 0.005    #!Viscous resistance of TV
    Rpv = 0.005   # !Viscous resistance of PV
    yav = 0.0005   #!Inertance of AV
    ymv = 0.0002   #!Inertance of MV
    ytv = 0.0002   #!Inertance of TV
    ypv = 0.0005   #!Inertance of PV

    #! Pulmonary circulation
    Epua0 = 0.02
    Epuc0 = 0.02
    Epuv0 = 0.02
    Epwc0 = 0.7
    Epwv0 = 0.7
    Rpua = 0.04
    Rpuc = 0.04
    Rpuv = 0.005
    Rpwa = 0.0005
    Rpwc = 0.4
    Rpwv = 0.4
    ypua = 0.0005
    ypuc = 0.0005
    ypuv = 0.0005
    ypwa = 0.0005
    ypwc = 0.0005
    ypwv = 0.0005
    Zpua = 20.0
    Zpuc = 60.0
    Zpuv = 200.0
    Zpwa = 1.0
    Zpwc = 1.0
    Zpwv = 1.0
    Spua = 0.01
    Spuc = 0.01
    Spuv = 0.01
    Spwa = 0.01
    Spwc = 0.01
    Spwv = 0.01

    #! Peripheral circulation
    Caor = 0.9
    Cart = 0.3
    Ccap = 0.06
    Cven = 100.0
    Cvca = 30.0

    yaor = 0.005
    yart = 0.001
    ycap = 0.0005
    yven = 0.0005
    yvca = 0.0005
    
    Raor = 0.03
    #Rart = 0.75
    Rcap = 0.35
    Rven = 0.07
    Rvca = 0.001
    
    Saor=0.01
    Sart=0.01
    Scap=0.01
    Sven=0.01
    Svca=0.01
    qco=0.0

    ddt=0.001
    
    for cycle in range(65):
        if cycle==0:
            Tduration=0.855
            Rart=0.75
            Elva=2.87
            Erva=0.52
        else:
            Tduration=Tbaro
            Rart=Rartbaro
            Elva=Elvbaro
            Erva=Ervbaro

        for nstep in np.arange(0 , Tduration, ddt):
            if nstep==0:
                ppc=0.0
                result=np.zeros(shape=(2, 101))
                result[0, 0:35] = odic

                for i in range(14):
                    if i<=6:
                        v[0, i] = result[0,2*i]
                    else:
                        v[0, i] = result[0,2*i+1]
    
                for i in range(15):
                    if i<=6:
                        q[0, i] = result[0,2*i+1]
                    else:
                        q[0, i] = result[0,2*i]
                q[0, 15:21] = result[0, 29:35]
                STR = 1.0
                FL = 1.0
                FR1 = 1.0
                Gpw = 0.0
                Aav = 0.0
                Amv = 0.0
                Apv = 0.0
                Atv = 0.0
                Pit = -2.5
                P_0d = np.zeros(shape=(2,101))
                aod = -41.7
                aao = 73
                
            tcr = nstep
            tee = 0.3*np.sqrt(Tduration)
            tac = Tduration - 0.5*tee - 0.02*(Tduration/0.855)    
            tar = Tduration - 0.02*(Tduration/0.855)    
            
            ncount=0
            ncountadd=ncount+1
            resultcr[0, 0:101]=result[0, 0:101]
            
            #%c.... Compute the elastances
            Epua = Ecal(Epua0, Zpua, v[0, 4])
            Epuc = Ecal(Epuc0, Zpuc, v[0, 5])
            Epuv = Ecal(Epuv0, Zpuv, v[0, 6])
            Epwc = Ecal(Epwc0, Zpwc, v[0, 7])
            Epwv = Ecal(Epwv0, Zpwv, v[0, 8])
            
            #%c.....Update nolinear cardiac parameters
            if tcr==0.0:
                FL=1.0-(result[0,21] / Vmax)
                FR1=1.0-(result[0,6] / Vmax)
               
            Lvecal()
            Laecal()
            Rvecal()
            Raecal()
   
            cklr = erv/(Es+erv)
            ckrl = elv/(Es+elv)
            plv = (ckrl*Es * v[0, 10] + ckrl * cklr * Es * v[0, 3]) / (1.0-cklr)
            prv = (cklr * Es * v[0, 3] + ckrl * cklr * Es* v[0, 10]) / (1.0-ckrl)
            Sla = Sva0 * v[0, 9] * ela
            Slv = Sva0 * plv
            Sra = Sva0 * v[0, 2] * era
            Srv = Sva0 * prv
            ppp = (v[0, 2] + v[0, 3] + v[0, 9] + v[0, 10] + Vpe-Vpc0) / Vcon
            ppc = np.exp(ppp)
    
            fmin = 2.52
            fmax = 47.78
            Pn = 92
            ka = 11.75
            fcs1 = (fmin + fmax * np.exp((q[0, 15]-Pn) / ka))/ (1 + np.exp((q[0, 15]-Pn) / ka))
            fcs = np.around(fcs1, 3)    # General expression is (fsc1, -3), numpy operates differently. 

            kes=0.0675
            fes0=16.11
            fes1=2.1
            es1=(-kes)*fcs
            es2=np.exp(es1)
            es3=(fes0-fes1)
            fes1=fes1+(es3*es2)
            fes=np.around(fes1, 3)
    
            fev0=3.2
            fev1=6.3
            kev=7.06
            fcs0=25
            fev1=(fev0+fev1*np.exp((q[0, 15]-fcs0)/kev))/(1+np.exp((q[0, 15]-fcs0)/kev))
            fev=np.around(fev1, 3)
           
            #================Arterial Resistance=====================
            if jj>2000:
                fesmin = 2.71
                if fes > fesmin:
                    gainR1= fes - fesmin+1
                    gainR2=np.log(gainR1)
                    GR=0.135
                    gainR=GR*gainR2
                else:
                    gainR=0
            else:
                gainR=0

            #================== LV Elastance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fesmin1=2.66
            Gelv=0.5
            if jj>2000:
                if fes>fesmin:
                    gainElv1 = fes-fesmin1+1
                    gainElv2 = np.log(gainElv1)
                    gainElv = Gelv*gainElv2
                else:
                    gainElv=0
            else:
                gainElv=0

            #================== RV Elastance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fesmin1=2.66
            Gerv=0.09
            if jj>2000:
                if fes>fesmin:
                    gainErv1 = fes-fesmin1+1
                    gainErv2 = np.log(gainErv1)
                    gainErv = Gerv*gainErv2
                else:
                    gainErv=0
            else:
                gainErv=0

            #====================Cardiac Cycle sympathetic%%%%%%%%%%%%%%%%%%%%
            fesmin1=2.66
            GTs= -0.13
            if jj>2000:
                if fes>fesmin:
                    gainTs1 = fes-fesmin1 + 1
                    gainTs2 = np.log(gainTs1)
                    gainTs=GTs*gainTs2
                else:
                    gainTs=0
            else:
                gainTs=0

            #====================Cardiac cycle-vaga;%%%%%%%%%%%%%%%%%
            Gtv=0.09
            if jj>200:
                gainTv=Gtv*fev
            else:
                gainTv=0
    
            Ervbaro1=q[0, 18]+0.3956
            Ervbaro=np.around(Ervbaro1, 3)
            
            Tbaro1=q[0, 19]+q[0, 20]+0.58
            Tbaro=np.around(Tbaro1, 3)
            
            Elvbaro1=q[0, 17]+2.0829
            Elvbaro=np.around(Elvbaro1, 3)
            
            Rartbaro1=0.534+q[0, 16]
            Rartbaro=np.around(Rartbaro1, 3)
            
            # Update dv
            Integrated_ode()
            
            diffv[0, 0:7] = dvq[0, 0:13:2]
            diffv[0, 7:14] = dvq[0, 15:28:2]
            dv[0, 0:14]=diffv[0, 0:14]
            
            rukuk=rkbaro(resultcr)
            
            #c.....Update variables with Runge-Kutta method
            for i in range(35):
                result[ncountadd, i] = result[ncount, i] + (rukuk[0, i] + 2.0*(rukuk[1, i] + rukuk[2, i]) + (rukuk[3, i])) / 6.0
             
            #%c.....Specify blood flows through cardiac valves
            delta = 0.00000001
            if Aav==0.0 and resultcr[0, 22]<=delta:
                resultcr[0, 22] = 0.0
                result[ncountadd,22] = 0.0
                q[0, 11]=0.0
    
            if(np.abs(tcr-Tduration) <= ddt*0.5 or tcr<0.1 or (Amv==0.0 or resultcr[0, 20]<=delta)):
                resultcr[0, 20]=0.0
                result[ncountadd, 20]=0.0
                q[0, 10]=0.0
    
            if(Apv==0.0 and resultcr[0, 7]<=delta):
                resultcr[0, 7]=0.0
                result[ncountadd, 7]=0.0
                q[0, 3]=0.0
    
            if(np.abs(tcr-Tduration) <= ddt*0.5 or tcr<0.1 or (Atv==0.0 or resultcr[0, 5]<=delta)):
                resultcr[0, 5]=0.0
                result[ncountadd,5] = 0.0
                q[0, 2] = 0.0

            v[0, 0:7] = result[ncountadd,0:13:2]
            v[0, 7:14] = result[ncountadd,15:28:2]
            
            q[0, 0:7] = result[ncountadd,1:14:2]
            q[0, 7:15] = result[ncountadd,14:29:2]
            q[0, 15:21] = result[ncountadd,29:35:1]
            result[0, 0:35] = result[1, 0:35]
            jj=jj+1
    
            MyResult[jj-1, 0:27] = P_0d[0, 0:27]
            
            MyResult1[jj-1, 0] = fcs
            MyResult1[jj-1, 1:36] = result[0, 0:35]###################
            
            MyResult2[jj-1, 0] = fcs
            MyResult2[jj-1, 1] = fes
            MyResult2[jj-1, 2] = fev
            MyResult2[jj-1, 4] = gainElv
            MyResult2[jj-1, 5] = elv
            #MyResult2[jj, 7] = Elvbaro
            #MyResult2[jj, 8] = Ervbaro
            #MyResult2[jj, 9] = Tbaro
            #MyResult2[jj, 10] = Rartbaro
            #MyResult2(jj,11)=nstep1;
            MyResult2[jj-1,12] = nstep
            MyResult2[jj-1,13] = tcr
    

            if jj>1:
                a = jj-2
                aod=(MyResult[jj-1, 22]-MyResult[a,22])*1000
            MyResult2[jj-1, 3] = aod
            
            if(nstep == Tduration):
                odic_new = result[0, 0:35]
                
        odic = np.around(odic_new, 3)
        
        np.savetxt('MyResult.txt', MyResult, delimiter=",")
        np.savetxt('MyResult1.txt', MyResult1, delimiter=",")
        np.savetxt('MyResult2.txt', MyResult2, delimiter=",")
          
        
        
if __name__ == "__main__":
    import time
    st = time.time()
    barocomplete()
    end = time.time()
    print("Total time is: ", (end-st))
        