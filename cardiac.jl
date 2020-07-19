function Integrated_ode()
    global Aav, Amv, Apv, Atv, Gpw = zeros(1,5)
    global Epua, Epuc, Epuv, Epwc, Epwv =zeros(5)
    global yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv = zeros(10)
    global Rav, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwc, Rpwv, Rtv, bav, bmv, bpv, btv = zeros(13)
    global Spua, Spuc, Spuv, Spwc, Spwv = zeros(5)
    global Zpua, Zpuc, Zpuv, Zpwc, Zpwv = zeros(5)
    global Caor, Cart, Ccap, Cven, Cvca = zeros(5)
    global yaor, yart, ycap, yven, yvca = zeros(5)
    global Raor, Rart, Rcap, Rven, Rvca = zeros(5)
    global Saor, Sart, Scap, Sven, Svca = zeros(5)
    global dvq, P_0d = zeros(2)
    global dv, v, q = zeros(3)
    global ela, era = zeros(2)
    global plv, prv, Sla, Slv, Sra, Srv, ppc, pit, qco, Tduration, ddt, tcr = zeros(12)
    
    dvq[1] = q[15] - q[1]  # Venous volume  dvq(1)= q(15) - q(1);
    
    P_0d[1] = v[1]/Cven + Sven*dv[1]  # Venous  Pressure
    
    dvq[2] = (v[1]/Cven + Sven*dv[1] - Rven*q[1]-v[2] / Cvca-Svca*dv[2]) / yven  # Venous flow
    
    dvq[3]=q[1]-q[2] #%VC volume

    P_0d[2]=v[2]/Cvca + Svca*dv[2]#% VC Pressure

    dvq[4]=(v[2]/Cvca-era*v[3]-Rvca*q[2]+Svca*dv[2] -Sra*dv[3]- ppc - pit)/yvca; #% VC Flow

    qco = 0.0

    dvq[5]=q[2]+qco-q[3] #% RA volume


    P_0d[4] = era * v[3] + Sra * dv[3] + ppc + pit  # RA pressure
    
    dvq[6] = (era * v[3] - prv - Rtv * q[3] - btv * q[3] * abs(q[3]) + Sra * dv[3] - Srv * dv[4]) / ytv  # TV flow
    
    dvq[7] = q[3] - q[4]  # RV volume

    P_0d[6] = prv + Srv * dv[4] + ppc + pit  # RV pressure

    dvq[8] = (prv - Epua * Zpua - Rpv * q[4] - bpv * q[4] * (abs(q[4])) + Srv * dv[4] - Spua * dv[5] + ppc) / ypv  # PV flow

    dvq[9] = q[4] - q[5] - q[8]  # Pulmonary Artery volume
    
    P_0d[8] = Epua * Zpua + Spua * dv[5] + pit  # Pulmonary Artery Pressure
    
    dvq[10] = (Epua * Zpua - Epuc * Zpuc - Rpua * q[5] + Spua * dv[5] - Spuc * dv[6]) / ypua  # Pulmonary Artery Flow
    
    dvq[11] = q[5] - q[6]  # Pulmonary Capillary volume
    
    P_0d[10] = Epuc * Zpuc + Spuc * dv[6] + pit  # Pulmonary Capillary Pressure
    
    dvq[12] = (Epuc * Zpuc - Epuv * Zpuv - Rpuc * q[6] + Spuc * dv[6] - Spuv * dv[7]) / ypuc  # Pulmonary Capillary Flow
    
    dvq[13] = q[6] - q[7]  # Pulmonary Vein volume
    
    P_0d[12] = Epuv * Zpuv + Spuv * dv[7] + pit  # Pulmonary Vein Pressure
    
    dvq[14] = (Epuv * Zpuv - ela * v[8] - Rpuv * q[7] + Spuv * dv[7] - Sla * dv[10] - ppc) / ypuv  # Pulmonary Vein flow
    
    if Gpw > 0
        dvq[15] = (Epua * Zpua - Epwc * Zpwc - q[8] / gpw + Spua * dv[5] - Spwc * dv[8]) / ypwa  # Pulmonary Wedge ##Artery flow
    else
        dvq[15] = 0
    end

    dvq[16] = q[8] - q[9]  # Pulmonary Wedge capillary volume
    
    P_0d[15] = Epwc * Zpwc + Spwc * dv[8] + pit  # Pulmonary Wedge capillary Pressure
    
    dvq[17] = (Epwc * Zpwc - Epwv * Zpwv - Rpwc * q[9] + Spwc * dv[8] - Spwv * dv[9]) / ypwc  # Pulmonary Wedge capillary Flow
    
    dvq[18] = q[9] - q[10]  # Pulmonary Wedge vein Volume
    
    P_0d[17] = Epwv * Zpwv + Spwv * dv[9] + pit  # Pulmonary Wedge vein Pressure
    
    dvq[19] = (Epwv * Zpwv - ela * v[10] - Rpwv * q[10] + Spwv * dv[9] - Sla * dv[10] - ppc) / ypwv  # Pulmonary Wedge vein Flow
    
    dvq[20] = q[7] + q[8] - q[11]  # LA volume
    
    P_0d[19] = ela * v[10] + Sla * dv[10] + ppc + pit  # LA Pressure
    
    dvq[21] = (ela*v[10] - plv-Rmv * q[11]-bmv * q[11]*abs(q[11]) + Sla*dv[10] - Slv *dv[11]) / ymv  # Mitral flow
    
    dvq[22] = q[11] - q[12]  # LV volume

    P_0d[21] = plv + Slv * dv[11] + ppc + pit  # LV Pressure
    
    dvq[23] = (plv - v[12] / Caor - Rav * q[12] - (bav * q[12]) * abs.(q[12]) + Slv * dv[11] - Saor * dv[12] + ppc + pit) / yav  # Aortic Flow
    
    # Adjust the state of cardiac valve
    if Aav == 0.0  && q[12] <= 0.000000001
        dvq[23] = 0.0
    end

    if abs.(tcr - Tduration) <= ddt * 0.5  || tcr < 0.1  || (Amv == 0.0  && q[11] < 0.000000001)
        dvq[21] = 0.0
    end

    if Apv == 0.0  && q[4] <= 0.00000001
        dvq[8] = 0.0
    end

    if abs.(tcr-Tduration) <= ddt*0.5 ||tcr<0.1  || (Atv == 0.0  && q[3] <= 0.000000001)
        dvq[6] = 0.0
    end

    dvq[24] = q[12] - q[13]  # Aorta volume

    P_0d[23] = v[12] / Caor + Saor * dv[12]  # Aorta Presuure
   
    dvq[25] = (v[12] / Caor + Saor * dv[12] - q[13] * Raor - v[13] / Cart - Sart * dv[13]) / yaor  # Aorta Flow
    
    dvq[26] = q[13] - q[14]  # Artery Volume
    
    P_0d[25] = v[13] / Cart + Sart * dv[13]  # Artery Pressure
    
    dvq[27] = (v[13] / Cart + Sart * dv[13] - q[14] * Rart - v[14] / Ccap - Scap * dv[14]) / yart  # Artery Flow
    
    dvq[28] = q[14] - q[15]  # Capillarya Volume
    
    P_0d[27] = v[14] / Ccap + Scap * dv[14]  # Capillary Pressure
    
    dvq[29] = (v[14] / Ccap + Scap * dv[14] - q[15] * Rcap - v[1] / Cven - Sven * dv[1]) / ycap  # Capillary Pressure
    
end

function Ecal(EEE, ZZZ, vol)
    EcalR = EEE * exp(vol / ZZZ)
    return EcalR
end

function Lvecal()
    global Elva, Elvb = zeros(2)
    global elv, FL = zeros(2)
    global tee, tcr = zeros(2)

    tcal = tcr
    if tcal <= tee
        elv = FL * Elva * 0.5 * (1.0 - cos(3.1415926 * tcal / tee)) + Elvb / FL
    else
        if tcal <= 1.5 * tee
            elv = FL * Elva * 0.5 * (1.0 + cos(3.1415926 * (tcal - tee) / (0.5 * tee))) + Elvb / FL
        else
            elv = Elvb / FL
        end
    end
end

function Laecal()
    global tcr, ela, Elaa, tac, tar, Tduration, Elab= zeros(7)
    tcal = tcr
    teec = tar - tac
    teer = teec
    tap = tar + teer - Tduration

    if (tcal >= 0.0  && tcal <= tap)
        ela = Elaa * 0.5 * (1.0 + cos(3.1415926 * (tcal + Tduration - tar) / teer)) + Elab
    end

    if (tcal > tap && tcal <= tac)
        ela = Elab
    end

    if (tcal > tac && tcal <= tar)
        ela = Elaa * 0.5 * (1.0 - cos(3.1415926 * (tcal - tac) / teec)) + Elab
    # c if (tcal > tar. and.tcal <= (tar+teer)) then
    end

    if (tcal > tar && tcal <= Tduration)
        ela = Elaa * 0.5 * (1.0 + cos(3.1415926 * (tcal - tar) / teer)) + Elab
    end

end

function Rvecal()
    global tcr, FR1, Erva, tee, Ervb, erv= zeros(6)
    tcal = tcr
    if tcal <= tee
        erv = FR1 * Erva * 0.5 * (1.0 - cos(3.1415926 * tcal / tee)) + Ervb / FR1
    else
        if tcal <= 1.5 * tee
            erv = FR1 * Erva * 0.5 * (1.0 + cos(2.0 * 3.1415926 * (tcal - tee) / tee)) + Ervb / FR1
   else
             erv = Ervb / FR1
         end
     end
end


function Raecal() 
    global Eraa, Erab= zeros(2)
    global tar, tac, tcr= zeros(3)
    global era, Tduration= zeros(2)
    teec = tar - tac
    teer = teec
    tcal = tcr
    tap = tar + teer - Tduration

    if 0 <= tcal <= tap 
        era = Eraa * 0.5 * (1.0 + cos(3.1415926 * (tcal + Tduration - tar) / teer)) + Erab
    end

    if tap < tcal <= tac 
        era = Erab
    end

    if tcal > tac  && tcal <= tar 
        era = Eraa * 0.5 * (1.0 - cos(3.1415926 * (tcal - tac) / teec)) + Erab
    end

    if tcal > tar  && tcal <= Tduration 
        era = Eraa * 0.5 * (1.0 + cos(3.1415926 * (tcal - tar) / teer)) + Erab
    end
end

function AAav() 
    global Caor, dv, v= zeros(3)
    global plv, Slv, ppc= zeros(3)
    intee = plv + Slv * dv[11] + ppc - v[12] / Caor
    
    if intee > 0.0 
        AAav = 4.0
    else 
        AAav = 0.0
    end
        
    return AAav
end

function AAmv() 
    global v, ela, plv = zeros(3)
    intee = ela * v[10] - plv
    if intee > 0.0 
        AAmv = 4.0
    else 
        AAmv = 0.0
    end
    
    return AAmv
end

function AApv() 
    global Epua, Zpua, prv = zeros(3)
    intee = prv - Epua * Zpua
    if intee > 0.0 
        AApv = 4.0
    else 
        AApv = 0.0
    
    end
    return AApv
end

function AAtv() 
    global era, v, prv = zeros(3)
    intee = era * v[3] - prv
    if intee > 0.0 
        AAtv = 4.0
    else 
        AAtv = 0.0
    end    
    return AAtv
end

function cardiac_state(subresultcr)
    
    global Aav, Amv, Apv, Atv, dvq = zeros(5)
    global dv, v, dq, q, ddt = zeros(5)
    # call Integrated_ode % == == == == == = make a note........................
    Integrated_ode()
    # Declaration of variables
    global c = 0
    dfl = zeros(2, 30)
    dq = zeros(2, 15)
    global subrukuk = zeros(4, 29)
    global inter = zeros(2, 29)
    newpara = zeros(2, 29)
    # -------------------------------------------------------------------------
    dfl[1:29] = dvq[1:29]
    dq[1:7] = dfl[1:2:14]
    dq[8:15] = dfl[14:2:29]

    dv[1:7] = dfl[1:2:13]
    dv[8:14] = dfl[15:2:28]

    for nrk in 4
        global subrukuk
        global inter
        subrukuk[nrk, 1:29] = (ddt* dfl[1,1:29])

        if nrk == 1 
            Aav = AAav()
            Amv = AAmv()
            Apv = AApv()
            Atv = AAtv()
        end
        if nrk < 4
            inter[1, 1:29] = (0.5 * subrukuk[nrk, 1:29])
        else
            inter[1, 1:29] = subrukuk[nrk, 1:29]
        end
        c = nrk
    end
    
    inter[1:29] = (0.5 * subrukuk[c, 1:29])
    newpara[1:29] = (subresultcr[1:29] + inter[1:29])
    q[1:7] = newpara[1:2:14]
    q[8:15] = newpara[14:2:29]
    v[1:7] = newpara[1:2:13]
    v[8:14] = newpara[15:2:28]

    return subrukuk
end
# End of function definition ----------------------------------

function cardiac()
    global Aav, Amv, Apv, Atv, Gpw    =zeros(1,5)                    # valve, 
    global Elaa, Elab, Elva, Elvb, Eraa, Erab, Erva, Ervb, Epua, Epuc, Epuv, Epwa, Epwc, Epwv= zeros(14) # E_cardiopul,
    global yav, ymv, ypv, ytv, ypua, ypuc, ypuv, ypwa, ypwc, ypwv = zeros(10) # yL_cardiopul,
    global Ra, Raa, Rav, Rca, Rda, Rmv, Rpua, Rpuc, Rpuv, Rpv, Rpwa, Rpwc, Rpwv, Rtv, Rv, Rvc = zeros(16)
    global bav, bmv, bpv, btv = zeros(4) # R_cardiopul
    global Spua, Spuc, Spuv, Spwa, Spwc, Spwv = zeros(6)            # S_cardiopul,
    global  Zpua, Zpuc, Zpuv, Zpwa, Zpwc, Zpwv = zeros(6)           # Z_cardiopul,
    global  Caor, Cart, Ccap, Cven, Cvca = zeros(5)                 #
    global  yaor, yart, ycap, yven, yvca = zeros(5)                 # yL_peripheral,
    global  Raor, Rart, Rcap, Rven, Rvca = zeros(5)                    # R_peripheral,
    global  Saor, Sart, Scap, Sven, Svca = zeros(5)                 # S_peripheral,
    global  dvq, P_0d = zeros(2)                                    # sdvsdqdvdq
    global  dv, v, dq, q = zeros(4)                                 # dvdq_cardiopul,
    global  elv, ela, erv, era, cklr, ckrl, plv, prv, Sla, Slv, Sra, Srv, ppp, ppc, pit, qco, FL, FR1, S0TR = zeros(19) # cardiac_parameter,
    global  Rav0, Rmv0, Rpv0, Rtv0, bav0, bmv0, bpv0, btv0, Rav1, Rmv1, Rpv1, Rtv1, bav1, bmv1, bpv1, btv1 = zeros(16)# R_cardiopulc,
    global Rav2, Rmv2, Rpv2, Rtv2, bav2, bmv2, bpv2, btv2, yav0, ymv0, ypv0, ytv0, yav1, ymv1, ypv1, ytv1, yav2, ymv2, ypv2, ytv2= zeros(20)#
    global n_val, m_cvst, m_cvrg, n_vrg = zeros(4)
    global timestep, Tduration, ddt, tee, tcr, tac, tar, t, odic = zeros(9)


    #%-----------------------------------------
    Elva=2.87          #   !Peak-systolic elastance of left ventricle
    Elvb=0.06          #   !Basic diastolic elastance of left ventricle
    Elaa=0.07          #   !Peak-systolic elastance of left atrium
    Elab=0.075         #   !Basic diastolic elastance of left atrium
    Erva=0.52          #   !Peak-systolic elastance of right ventricle
    Ervb=0.043         #   !Basic diastolic elastance of right ventricle
    Eraa=0.055         #   !Peak-systolic elastance of right atrium
    Erab=0.06          #   !Basic diastolic elastance of right atrium

    Vmax=900.0         #   !Reference volume of Frank-Starling law
    Es=45.9            #   !Effective septal elastance
    Vpc0=380.0         #   !Reference total pericardial and cardiac volume
    Vpe=30.0           #   !Pericardial volume of heart
    Vcon=40.0          #   !Volume constant
    Sva0=0.0005        #   !Coefficient of cardiac viscoelasticity

    # ! Cardiac valve parameters
    # !(aortic valve(AV),mitral valve(MV), tricuspid valve(TV),pulmonary valve(PV))
    bav=0.000025        #%!Bernoulli's resistance of AV
    bmv=0.000016        #%!Bernoulli's resistance of MV
    btv=0.000016        #%!Bernoulli's resistance of TV
    bpv=0.000025        #%!Bernoulli's resistance of PV
    Rav=0.005           #%!Viscous resistance of AV
    Rmv=0.005           #%!Viscous resistance of MV
    Rtv=0.005           #%!Viscous resistance of TV
    Rpv=0.005           #% !Viscous resistance of PV
    yav=0.0005          #%!Inertance of AV
    ymv=0.0002          #%!Inertance of MV
    ytv=0.0002          #%!Inertance of TV
    ypv=0.0005          #%!Inertance of PV

    #! Pulmonary circulation
    Epua0=0.02
    Epuc0=0.02
    Epuv0=0.02
    Epwc0=0.7
    Epwv0=0.7
    Rpua=0.04
    Rpuc=0.04
    Rpuv=0.005
    Rpwa=0.0005
    Rpwc=0.4
    Rpwv=0.4
    ypua=0.0005
    ypuc=0.0005
    ypuv=0.0005
    ypwa=0.0005
    ypwc=0.0005
    ypwv=0.0005
    Zpua=20.0
    Zpuc=60.0
    Zpuv=200.0
    Zpwa=1.0
    Zpwc=1.0
    Zpwv=1.0
    Spua=0.01
    Spuc=0.01
    Spuv=0.01
    Spwa=0.01
    Spwc=0.01
    Spwv=0.01

    #%! Peripheral circulation
    Caor=0.9
    Cart=0.3
    Ccap=0.06
    Cven=100.0
    Cvca=30.0
    yaor=0.005
    yart=0.001
    ycap=0.0005
    yven=0.0005
    yvca=0.0005

    Raor=0.03
    Rart=0.75
    Rcap=0.35
    Rven=0.07
    Rvca=0.001

    Saor=0.01
    Sart=0.01
    Scap=0.01
    Sven=0.01
    Svca=0.01
    qco = 0.0
    #--------------------------------------------------------------------
    v = zeros(1, 14)
    q =  zeros(1, 21)
    dvq =  zeros(1, 29)
    result =  zeros(2, 101)
    P_0d =  zeros(2, 101)
    diffv =  zeros(1, 14)
    resultcr =  zeros(2, 101)
    dv =  zeros(1, 21)
    MyResult =  zeros(100000, 28)
    MyResult1 =  zeros(100000, 29)
    rukuk = zeros(4,29)

    # initial values of all state equations (should be equal to number of equations in dvdqsdvdq.m file)
    odic = [1068.2371, 52.4983, 181.7233, -41.7618, 65.0625, 0, 122.6637, 0, 67.0272, -0.3118, 135.1100, -2.1737, 198.7568,
         -64.1791, 0, 2.7983, 0.1357, 2.7932, 1.1042, 68.8587, 0, 121.2539, 0, 67.3641, 41.8262, 22.0472, 56.6627,
         1.8539, 57.6473]

    Pit = -2.5
    pit = Pit
    jj = 0

    # ------------------------------------------------------------------------------------------------
    HR = 70
    Tduration = 60 / HR                                     # input('Please specify cardiac duration(s)')
    dt = 0.001                                              # input'Please specify time step(s)')
    ncycle = 8
    tee = 0.3 * sqrt(Tduration)                             # !Moment when ventricular contractility reaches the peak
    tac = Tduration-0.5 * tee-0.02 * (Tduration/0.855)      # !Moment when atrium begins to contract
    tar = Tduration - 0.02 * (Tduration / 0.855)            # !Moment when atrium begins to relax
    ddt = dt
    ntotal = (ncycle * Tduration / dt)

    # ------------------------------------------------------------------------
    for nstep in 1:ntotal
        if nstep == 1
            tcr = 0.0
            ppc = 0.0
            cn = 0

            for i in 29
                result[1, i] = odic[i]
            end

            for i in 14
                if i <= 7
                    v[i] = result[2 * i - 1]
                else 
                    v[i] = result[2 * i ]
                end
            end

            for i in 15
                if i <= 7 
                    q[i] = result[2 * i]
                else 
                    q[i] = result[2 * i - 1]
                end
            end

            STR = 1.0
            FL = 1.0
            FR1 = 1.0
            Gpw = 0.0
            Aav = 0.0
            Amv = 0.0
            Apv = 0.0
            Atv = 0.0
            Pit = -2.5
        end

        # ----------------------Start computation---------------------------
        ncount = 1
        ncountadd = ncount + 1

        tcr = nstep * dt % Tduration
        t = nstep * dt

        # c.... Compute the pulmonary elastances
        Epua = Ecal(Epua0, Zpua, v[5])
        Epuc = Ecal(Epuc0, Zpuc, v[6])
        Epuv = Ecal(Epuv0, Zpuv, v[7])
        Epwc = Ecal(Epwc0, Zpwc, v[8])
        Epwv = Ecal(Epwv0, Zpwv, v[9])


        # c.....Update nolinear cardiac parameters
        if tcr == 0.0 
            FL = 1.0 - (result[1, 22] / Vmax)        # Left ventricle scaling factor
            FR1 = 1.0 - (result[1, 7] / Vmax)        # Right ventricle scaling factor
        end
    
        Lvecal()                                     # LV elastance function calling
        Laecal()                                     # LA elastance function calling
        Rvecal()                                     # RV elastance function calling
        Raecal()                                     # RA elastance function calling

        # Spetum cross talk pressure calculations
        cklr = erv / (Es+erv)
        ckrl = elv / (Es+elv)
        plv = (ckrl*Es*v[11]+ckrl*cklr*Es*v[4])/(1.0 - cklr)
        prv = (cklr*Es*v[4]+ckrl*cklr*Es*v[11])/(1.0 - ckrl)

        # All cardiac chambers viscoelastance calculation
        Sla = Sva0 * v[10] * ela
        Slv = Sva0 * plv
        Sra = Sva0 * v[3] * era
        Srv = Sva0 * prv

        # Pericardium pressure calculations
        ppp = (v[3] + v[4] + v[10] + v[11] + Vpe - Vpc0) / Vcon
        ppc = exp(ppp)

        # c....Update dv and state equation function calling
        Integrated_ode()  

        diffv[1:7] = dvq[1:2:13]
        diffv[8:14] = dvq[16:2:28]
        dv[1:14] = diffv[1:14]

        # c.....Implement fourth - order Runge - Kutta method
        rukuk = cardiac_state(resultcr)

        #println(size(rukuk))

        # c.....Update variables with Runge-Kutta method
        for j in 29
            result[ncountadd, j] = result[ncount, j] + (rukuk[1, j] + 2.0 * (rukuk[2, j] + rukuk[3, j]) + rukuk[4, j]) / 6.0
        end

        # %c.....Specify blood flows through cardiac valves
        delta = 0.00000001

        # update all four cardiac valve flow as zero when they were closed
        if Aav == 0.0  && resultcr[23] <= delta
        resultcr[23] = 0.0
        result[ncountadd, 23] = 0.0
        q[12] = 0.0
        end

        if abs.(tcr - Tduration) <= ddt * 0.5 || tcr < 0.1 || (Amv == 0.0  && resultcr[21]) <= delta 
        resultcr[21] = 0.0
        resultcr[ncountadd, 21] = 0.0
        q[11] = 0.0
        end

        if (Apv == 0.0)  && (resultcr[8] <= delta)
        resultcr[8] = 0.0
        resultcr[ncountadd, 8] = 0.0
        q[4] = 0.0
        end

        if abs.(tcr - Tduration) <= ddt * 0.5 || tcr < 0.1 || (Atv==0.0  && resultcr[6]<=delta) 
        resultcr[6] = 0.0
        result[ncountadd, 6] = 0.0
        q[3] = 0.0
        end

        # c.....Update q() and v() to be used at next time
        v[1:7] = result[ncountadd, 1:2:13]
        v[8:14] = result[ncountadd, 16:2:28]

        q[1:7] = result[ncountadd, 2:2:14]
        q[8:15] = result[ncountadd, 15:2:29]
        q[16:21] = result[ncountadd, 30:1:35]

        result[1, 1:29] = result[2, 1:29]

        # load calculated hemodynamic variables in an array
        jj += 1
        MyResult[jj, 1:27] = P_0d[1:27]
        MyResult1[jj, 1:29] = result[1, 1:29]

        # Update the initial conditions for next cardiac cycle
        if nstep == ntotal
            odic_new = result[1, 1:29]
        end

    end
    println("completed")
    return MyResult1
end



