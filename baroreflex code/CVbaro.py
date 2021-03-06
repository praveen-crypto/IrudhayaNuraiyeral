import numpy as np


Aav, Amv, Apv, Atv, Gpw = 0
Elaa, Elab, Elva, Elvb, Eraa, Erab, Erva, Ervb, Epua, Epuc, Epuv, Epwa, Epwc, Epwv = 0
yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv = 0
Ra, Raa, Rav, Rca, Rda, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwa, Rpwc, Rpwv, Rtv, Rv, Rvc, bav, bmv, bpv, btv = 0
Spua, Spuc, Spuv, Spwa, Spwc, Spwv = 0
Zpua, Zpuc, Zpuv, Zpwa, Zpwc, Zpwv = 0
Caor, Cart, Ccap, Cven, Cvca=0
yaor, yart, ycap, yven, yvca=0
Raor, Rart, Rcap, Rv, Rvc=0
Saor, Sart, Scap, Sven, Svca =0
dvq, P_0d =0
dv, v, dq, q=0
elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv, Sra, Srv, ppp, ppc, pit, qco, FL, FR1, STR=0
Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1, Rpv1, Rtv1, bav1, bmv1, bpv1, btv1=0
Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1, ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2=0
n_val, m_cvst, m_cvrg, n_vrg=0
timestep, Tduration, ddt, tee, tcr, tac, tar, t, odic=0
Rven,Rvca=0
aoa,nstep,aao,aod=0
gainElv,gainR,gainErv,gainTv,gainTs=0


def cvbaro():

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
        global  dvq ,P_0d , Pod ,nstep, gainR, gainElv, gainErv, gainTs, gainTv #sdvsdqdvdq
        global  dv, v, dq, q, aoa, aao, aod  #dvdq_cardiopul 
        global  elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv, Sra, Srv ,ppp, ppc, pit, qco, FL, FR1, STR  #cardiac_parameter
        global  Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1, Rpv1, Rtv1, bav1, bmv1, bpv1, btv1 #R_cardiopulc
        global Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1, ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2
        global n_val, m_cvst, m_cvrg, n_vrg
        global  Tduration, ddt, tcr, tee, tac, tar, t #timestep

        dvq[0, 0] = q[0, 14] - q[0, 0]  # Venous volume  dvq(1)= q(15) - q(1);

        P_0d[0, 0] = v[0, 0] / Cven + Sven * dv[0, 0]  # Venous  Pressure

        dvq[0, 1] = (v[0, 0] / Cven + Sven * dv[0, 0] - Rv * q[0, 0] - v[0, 1] / Cvca - Svca * dv[
                0, 1]) / yven  # Venous flow

        dvq[0, 2] = q[0, 0] - q[0, 1]  # VC volume

        P_0d[0, 1] = v[0, 1] / Cvca + Svca * dv[0, 1]  # VC Pressure

        dvq[0, 3] = (v[0, 1] / Cvca - era * v[0, 2] - Rvc * q[0, 1] + Svca * dv[0, 1] - Sra * dv[
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
            0, 0]) / ycap  # Capillary Pressure

        #% -------------------------------------------------------------------------
        #%y(25) Nbr:t q(16)
        #% dvq(30) = q(17);
        #% if nstep==0
        #%     Pod=-41.76;
        #%     'chek'
        #% else
        #% Pod=(diff(P_0d(23)))*1000;
        #% end

        dvq[30] = (aao+(6.37*aod)-q(16))/2.076;
        #%y(26) Nbr_t:t
        #%modeqn

        dvq[31] = (-q(17)+gainR)/6;
        #% dvq(31)   =(1.05*P_0d(23) + 0.036*1.05*((q(12)-q(13))/0.9)-(0.0018+0.001)*q(17) - q(16))/(0.0018*0.001); 
        #%modeqn
        dvq[32] = (-q(18)+gainElv)/8;
        dvq[33] = (-q(19)+gainErv)/8;
        dvq[34] = (-q(20)+gainTs)/2;
        dvq[35] = (-q(21)+gainTv)/1.5;
        #%y(27) N_hrv:t
        ##%modeqn
        #% if 0+t>0.2
        #% dvq(32)  =  (-q(18) + (1.2*q(16)))/1.8;
        #% else 
        #% dvq(32) = 0;
        #% end
        #%modeqn

        #%y(28) N_hrs:t
        #%modeqn
        #% if (0+t>3.0)
        #% dvq(33)  = (-q(19) + (1.15*q(16)))/10.0;
        #% else
        #% dvq(33) = 0;
        #% end
        #% %modeqn
        #% 
        #% %y(29) N_con:t
        #% %modeqn
        #% if (0+t>3.0)
        #% dvq(34)  =  (-q(20) + (1.35*q(16)))/10.0;
        #% else
        #% dvq(34) =0;
        #% end
        #% %modeqn
        #% 
        #% %q(21) N_vaso:t
        #% %modeqn
        #% if (0+t>3.0)
        #% dvq(35) =  (-q(21) + (1.35*q(16)))/6.0;
        #% else 
        #% dvq(35)  =0; 
        #% end
        #% %modeqn

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
        # c if (tcal > tar. and.tcal <= (tar+teer)) then
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
        
        # call Integrated_ode % == == == == == = make a note........................
        Integrated_ode()
        # Declaration of variables
        dfl = np.zeros(shape=(2, 30))
        dq = np.zeros(shape=(2, 15))
        subrukuk = np.zeros(shape=(4, 29))
        inter = np.zeros(shape=(2, 29))
        newpara = np.zeros(shape=(2, 29))
        # -------------------------------------------------------------------------
        dfl[1:35]=dvq[1:35]

        dq[1:7]   = dfl[2:2:14]
        dq[8:15]  = dfl[15:2:29]
        dq[16:21] = dfl[30:1:35]

        dv[1:7] = dfl[1:2:13]
        dv[8:14] = dfl[16:2:28]

        for nrk in range(4):
            subrukuk[nrk,1:35] = ddt*dfl[1:35]
            if nrk == 0:               
                Aav = AAav()
                Amv = AAmv()
                Apv = AApv()
                Atv = AAtv() #%=========>> make a note.........................&&
            if nrk<4:
                inter[1:35] = 0.5 * (subrukuk[nrk, 1:35])
            else:
                inter[1:35]=(subrukuk[nrk,1:35])
        inter[1:35]=0.5*(subrukuk[nrk,1:35])

        newpara[1:35]=subresultcr[1:35]+inter[1:35]
        q[1:7]=newpara[2:2:14]
        q[8:15]=newpara[15:2:29]
        q[16:21]=newpara[30:1:35]
        v[1:7]=newpara[1:2:13]
        v[8:14]=newpara[16:2:28]
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

    v = np.zeros(shape=(2, 14))
    q = np.zeros(shape=(2, 21))
    dvq = np.zeros(shape=(2, 29))
    result = np.zeros(shape=(2, 101))
    P_0d = np.zeros(shape=(2, 101))
    diffv = np.zeros(shape=(2, 14))
    resultcr = np.zeros(shape=(2, 101))
    dv = np.zeros(shape=(2, 21))
    MyResult = np.zeros(shape=(100000, 28))
    MyResult1 = np.zeros(shape=(100000, 29))

    odic = np.array(
            [1068.2371, 52.4983, 181.7233, -41.7618, 65.0625, 0, 122.6637, 0, 67.0272, -0.3118, 135.1100, -2.1737, 198.7568,
             -64.1791, 0, 2.7983, 0.1357, 2.7932, 1.1042, 68.8587, 0, 121.2539, 0, 67.3641, 41.8262, 22.0472, 56.6627,
             1.8539, 57.6473])

    dv = np.zeros(1, 21)
    Pit = -2.5
    pit = Pit
    jj = 0
    kk = 0

    # ======================================
    Elva =  2.87     #!Peak-systolic elastance of left ventricle
    Elvb = 0.06      #!Basic diastolic elastance of left ventricle
    Elaa = 0.07      #!Peak-systolic elastance of left atrium
    Elab = 0.075     #!Basic diastolic elastance of left atrium
    Erva = 0.52      #!Peak-systolic elastance of right ventricle
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
    Rart = 0.75
    Rcap = 0.35
    Rven = 0.07
    Rvca = 0.001
    
    Saor = 0.01
    Sart = 0.01
    Scap = 0.01
    Sven = 0.01
    Svca = 0.01
    qco =  0.0

    #%---------------------------------------
    #%Tduration=input('Please specify cardiac duration(s)');
    dt= input('Please specify time step(s)')
    #% tee=0.3*sqrt(Tduration); %!Moment when ventricular contractility reaches the peak
    #% tac=Tduration - 0.5*tee -  0.02* (Tduration/0.855);% !Moment when atrium begins to contract
    #% tar=Tduration - 0.02* (Tduration/0.855);% !Moment when atrium begins to relax
    ddt = dt
    #% ncycle=input('how many cardiac cycles to run ?');
    ntotal =30

    #---------------------------------------
    for nstep in range(0 , ddt, ntotal):
        if nstep==0:
            tcr=0.0
            ppc=0.0
            result=np.zeros[1,101]
            result[1, 1:35] = odic

            #       Rartbaro=0.75
            #       Ervbaro=0.58
            #       Elvbaro=2.87
            for i in range(14):
                if i<=7:
                    v[i] = result[1,2*i-1]
                else:
                    v[i] = result[1,2*i]

            for i in range(15):
                if i<=7:
                    q[i] = result[1,2*i]
                else:
                    q[i] = result[1,2*i-1]

            q[16:21] = result[30:1:35]
            STR = 1.0
            FL = 1.0
            FR1 = 1.0
            Gpw = 0.0
            Aav = 0.0
            Amv = 0.0
            Apv = 0.0
            Atv = 0.0
            Pit = -2.5
            P_0d = np.zeros(1,101)
            aod = -41.7
            aao = 73

            if nstep<0.855:
                Tduration=0.855
                tcr1=np.rem(nstep,Tduration)
                tcr= np.roundn(tcr1,-3)
                tee=0.3* np.sqrt(Tduration);                            #%!Moment when ventricular contractility reaches the peak
                tac=Tduration - 0.5*tee -  0.02* (Tduration/0.855); #% !Moment when atrium begins to contract
                tar=Tduration - 0.02* (Tduration/0.855);            #% !Moment when atrium begins to relax
            else:
                Tduration = Tbaro
                tcr1 = np.rem(nstep, Tduration)
                tcr = np.roundn(tcr1, -3)
                tee = 0.3*np.sqrt(Tduration);                               #%!Moment when ventricular contractility reaches the peak
                tac = Tduration - 0.5*tee -  0.02;  #%* (Tduration/0.855); #% !Moment when atrium begins to contract
                tar = Tduration - 0.02;             #%* (Tduration/0.855);            #% !Moment when atrium begins to relax
    
        ncount=1
        ncountadd=ncount+1
        resultcr[1:101]= result[1, 1:101]

        #%c.... Compute the elastances
        Epua = Ecal(Epua0,Zpua, v[5])
        Epuc = Ecal(Epuc0,Zpuc, v[6])
        Epuv = Ecal(Epuv0,Zpuv, v[7])
        Epwc = Ecal(Epwc0,Zpwc, v[8])
        Epwv = Ecal(Epwv0,Zpwv, v[9])

        #%%%%%%%%%%%%*******************Order changed**************%%%%%%%%%%%
        #%c.....Update nolinear cardiac parameters
        if tcr==0.0:
            FL=1.0-(result[1,22] / Vmax)
            FR1=1.0-(result[1,7] / Vmax)
        Lvecal()
        Laecal()
        Rvecal()
        Raecal()

        cklr = erv / (Es+erv)
        ckrl = elv/(Es+elv)
        plv = (ckrl * Es * v[11] + ckrl * cklr * Es * v[4]) / (1.0-cklr)
        prv = (cklr * Es * v[4] + ckrl * cklr * Es* v[11]) / (1.0-ckrl)
        Sla = Sva0 * v[10] * ela
        Slv = Sva0 * plv
        Sra = Sva0 * v[3] * era
        Srv = Sva0 * prv
        ppp = (v[3] + v[4] + v[10] + v[11] + Vpe-Vpc0) / Vcon
        ppc = np.exp(ppp)

        fmin = 2.52
        fmax = 47.78
        Pn = 92
        ka = 11.75
        fcs1 = (fmin + fmax * np.exp((q[16]-Pn) / ka))/ (1 + np.exp((q[16]-Pn) / ka))
        fcs = np.roundn(fcs1, -3)

        kes=0.0675
        fes0=16.11
        fes1=2.1
        es1=(-kes)*fcs
        es2=np.exp(es1)
        es3=(fes0-fes1)
        fes1=fes1+(es3*es2)
        fes=np.roundn(fes1,-3)

        fev0=3.2
        fev1=6.3
        kev=7.06
        fcs0=25
        fev1=(fev0+fev1*np.exp((q(16)-fcs0)/kev))/(1+np.exp((q(16)-fcs0)/kev));
        fev=np.roundn(fev1,-3)

    #================Arterial Resistance=====================

        if nstep>2:
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


        Rartbaro1 = 0.534 + q(17)
        Rartbaro = np.roundn(Rartbaro1,-3)


    #================== RV Elastance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fesmin1=2.66
        Gerv=0.09
        if nstep>2:
            if fes>fesmin:
                gainErv1 = fes-fesmin1+1
                gainErv2 = np.log(gainErv1)
                gainErv = Gerv*gainErv2
            else:
                gainErv=0
        else:
            gainErv=0

        Ervbaro1 = q(19)+0.3956
        Ervbaro = np.roundn(Ervbaro1,-3)

    #====================Cardiac Cycle sympathetic%%%%%%%%%%%%%%%%%%%%
        fesmin1=2.66
        GTs=-0.13
        if nstep>2:
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
        if nstep>0.200:
            gainTv=Gtv*fev
        else:
            gainTv=0

        Tbaro1 = q[20] +q[21]+0.58
        Tbaro = np.roundn(Tbaro1, -3)

        #%c.... Update dv
        Integrated_ode() #%============================>> make a note.............................
        diffv[0:7] = dvq[1:2:13]
        diffv[8:14] = dvq[16:2:28]
        dv[1:14]=diffv[1:14]

        rukuk=rkbaro(resultcr)

        #c.....Update variables with Runge-Kutta method
        for i in range(35):
            result[ncountadd,i] = result(ncount,i)+(rukuk[1,i] + 2.0 * (rukuk[2,i] + rukuk[3,i]) + (rukuk[4,i])) / 6.0;

        #%c.....Specify blood flows through cardiac valves
        delta = 0.00000001

        if Aav==0.0 and resultcr[23]<=delta:
            resultcr[23] = 0.0
            result[ncountadd,23] = 0.0
            q[12]=0.0

        if(abs(tcr-Tduration) <= ddt*0.5 or tcr<0.1 or (Amv==0.0 or resultcr[21]<=delta)):

            resultcr[21]=0.0
            result[ncountadd, 21]=0.0
            q[11]=0.0

        if(Apv==0.0 and resultcr[8]<=delta):
            resultcr[8]=0.0
            result[ncountadd,8]=0.0
            q[4]=0.0

        if(abs(tcr-Tduration) <= ddt*0.5 or tcr<0.1 or (Atv==0.0 or resultcr[6]<=delta)):
            resultcr[6]=0.0
            result[ncountadd,6] = 0.0
            q[3] = 0.0

        v[1:7] = result[ncountadd,1:2:13]
        v[8:14] = result[ncountadd,16:2:28]
        q[1:7] = result[ncountadd,2:2:14]
        q[8:15] = result[ncountadd,15:2:29]
        q[16:21] = result[ncountadd,30:1:35]
        result[1,1:35] = result[2,1:35]
        jj=jj+1

        #%           MyResult(jj,1)=t;

        MyResult[jj, 1:27] = P_0d[1:27]
        MyResult1[jj, 1] = fcs
        MyResult1[jj, 2:36] = result[1, 1:35]
        MyResult2[jj, 1] = fcs
        MyResult2[jj, 2] = fes
        MyResult2[jj, 3] = fev
        MyResult2[jj, 5] = gainElv
        MyResult2[jj, 6] = elv
        MyResult2[jj, 7] = Elvbaro
        MyResult2[jj, 8] = Ervbaro
        MyResult2[jj, 9] = Tbaro
        MyResult2[jj, 10] = Rartbaro
    #%     MyResult2(jj,11)=nstep1;
        MyResult2[jj,12] = tcr
        MyResult2[jj,13] = Tduration
    #%     MyResult2(jj,14)=tcheck;

        if jj>1:

            a = jj-1
            aod=(MyResult[jj, 23]-MyResult[a,23])*10000

        MyResult2[jj,4] = aod

        if(nstep == ntotal):
            odic_new = result[1, 1:35]

    #plot(Myresult2(:,1));

