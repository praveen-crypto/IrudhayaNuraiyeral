function rungekuttabaro(subresultcr)
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        global odic,jj, result, fesmin, nrk

        dfl = zeros(1,35)
        dq = zeros(1, 30)
        subrukuk = zeros(4, 35)
        inter = zeros(1,35)
        newpara = zeros(1,35)



        dvdqsdvsdq1()
        dfl[1:35]=dvq[1:35]

        dq[1:7]=dfl[2:2:14]
        dq[8:15]=dfl[15:2:29]
        dq[16:21]=dfl[30:1:35]

        dv[1:7]=dfl[1:2:13]
        dv[8:14]=dfl[16:2:28]

        for i=1:1:4
                global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
                global nrk
                subrukuk[i,1:35]=ddt*dfl[1:35]

                if (i==1)
                        Aav=AAav()
                        Amv=AAmv()
                        Apv=AApv()
                        Atv=AAtv()
                end
                if i<4
                        inter[1:35]=0.5*(subrukuk[i,1:35])
                else
                         inter[1:35]=(subrukuk[i,1:35])
                end
                nrk = i
        end
        inter[1:35]=0.5*(subrukuk[nrk,1:35])
        newpara[1:35]=subresultcr[1:35]+inter[1:35]

        q[1:7]=newpara[2:2:14]
        q[8:15]=newpara[15:2:29]
        q[16:21]=newpara[30:1:35]

        v[1:7]=newpara[1:2:13]
        v[8:14]=newpara[16:2:28]

        return subrukuk
end

function AAav()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        intee=plv + Slv*dv[11] + ppc - v[12]/Caor
        if intee>0.0
                AAav=4.0
        else
                AAav=0.0
        end
        return AAav
end

function AAmv()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        intee=ela*v[10]-plv
        if (intee>0.0)
                AAmv=4.0
        else
                AAmv=0.0
        end
        return AAmv
end

function AApv()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        intee=prv-Epua*Zpua
        if (intee>0.0)
            AApv=4.0
        else
            AApv=0.0
        end
        return AApv
end

function AAtv()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        intee=era*v[3]-prv
        if (intee>0.0)
                AAtv=4.0
        else
                AAtv=0.0
        end
        return AAtv
end

function Ecal(EEE,ZZZ,vol)
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        EcalR=EEE*exp(vol/ZZZ)
        return EcalR
end

function Laecal()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        tcal=tcr
        teec=tar-tac
        teer=teec
        tap=tar+teer-Tduration

        if (tcal>=0.0 && tcal<=tap)
                ela=Elaa*0.5*(1.0+cos(3.1415926*(tcal+Tduration-tar)/teer))+Elab     ;
        end
        if (tcal>tap && tcal<=tac)
                ela=Elab;
        end
        if (tcal>tac && tcal<=tar)
                ela=Elaa*0.5*(1.0-cos(3.1415926*(tcal-tac)/teec))+Elab;
        end
        # c      if(tcal>tar.and.tcal<=(tar+teer)) then
        if (tcal>tar && tcal<=Tduration)
                ela=Elaa*0.5*(1.0+cos(3.1415926*(tcal-tar)/teer))+Elab;
        end
end

function Lvecal()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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

        tcal = tcr
        if (tcal<=tee)
            elv=FL*Elva*0.5*(1.0-cos(3.1415926*tcal/tee))+Elvb/FL
        else
            if (tcal<=1.5*tee)
                elv=FL*Elva*0.5*(1.0+cos(3.1415926*(tcal-tee)/(0.5*tee))) +Elvb/FL
            else
                elv=Elvb/FL
            end
        end
end

function raecal()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        tcal=tcr
        teec=tar-tac
        teer=teec
        tap=tar+teer-Tduration

        if (tcal>=0 && tcal<=tap)
            era=Eraa*0.5*(1.0+cos(3.1415926*(tcal+Tduration-tar)/teer))+Erab
        end

        if (tcal>tap && tcal<=tac)
            era=Erab
        end

        if (tcal>tac && tcal<=tar)
            era=Eraa*0.5*(1.0-cos(3.1415926*(tcal-tac)/teec))+Erab
        end

        if (tcal>tar && tcal<=Tduration)
            era=Eraa*0.5*(1.0+cos(3.1415926*(tcal-tar)/teer))+Erab
        end
end

function rvecal()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        tcal=tcr;

        if (tcal<=tee)
            erv=FR1*Erva*0.5*(1.0-cos(3.1415926*tcal/tee))+Ervb/FR1;
        else
            if (tcal<=1.5*tee)
                erv=FR1*Erva*0.5*(1.0+cos(2.0*3.1415926*(tcal-tee)/tee))+Ervb/FR1
            else
                erv=Ervb/FR1
            end
        end
end

function dvdqsdvsdq1()
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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

        dvq[1]= q[15] - q[1] # Veneous volume

        P_0d[1]=v[1]/Cven + Sven*dv[1] #Venous  Pressure

        dvq[2]=(v[1]/Cven + Sven*dv[1] - Rven*q[1] - v[2]/Cvca - Svca*dv[2])/yven # Venous flow

        dvq[3]=q[1]-q[2] #VC volume

        P_0d[2]=v[2]/Cvca + Svca*dv[2]# VC Pressure

        dvq[4]=(v[2]/Cvca-era*v[3]-Rvca*q[2]+Svca*dv[2] -Sra*dv[3]- ppc - pit)/yvca # VC Flow

        qco = 0.0

        dvq[5]=q[2]+qco-q[3] # RA volume

        P_0d[4]=era*v[3] + Sra*dv[3]  + ppc + pit #RA pressure

        dvq[6]=(era*v[3]-prv-Rtv*q[3]-btv*q[3]*abs(q[3])+ Sra*dv[3]-Srv*dv[4])/ytv # TV flow

        dvq[7]=q[3]-q[4] #RV volume

        P_0d[6]=prv + Srv*dv[4] + ppc + pit #RV pressure

        dvq[8]=(prv-Epua*Zpua-Rpv*q[4]-bpv*q[4]*(abs(q[4]))+Srv*dv[4]-  Spua*dv[5]+ppc)/ypv #PV flow

        dvq[9]=q[4]-q[5]-q[8] #Pulmonary Artery volume

        P_0d[8]=Epua*Zpua + Spua*dv[5] + pit #Pulmonary Artery Pressure

        dvq[10]=(Epua*Zpua-Epuc*Zpuc-Rpua*q[5]+Spua*dv[5] -Spuc*dv[6])/ypua #Pulmonary Artery Flow

        dvq[11]=q[5]-q[6]#Pulmonary Capillary volume

        P_0d[10]=Epuc*Zpuc+Spuc*dv[6]+pit#Pulmonary Capillary Pressure

        dvq[12]=(Epuc*Zpuc-Epuv*Zpuv- Rpuc*q[6]+Spuc*dv[6]-Spuv*dv[7])/ypuc#Pulmonary Capillary Flow

        dvq[13]=q[6]-q[7] #Pulmonary Vein volume

        P_0d[12]=Epuv*Zpuv+Spuv*dv[7]+pit#Pulmonary Vein Pressure

        dvq[14]=(Epuv*Zpuv-ela*v[10]-  Rpuv*q[7]+Spuv*dv[7]-Sla*dv[10]-ppc)/ypuv #Pulmonary Vein flow

        if (Gpw>0)
                dvq[15]=(Epua*Zpua-Epwc*Zpwc- q[8]/gpw+Spua*dv[5]-Spwc*dv[8])/ypwa #Pulmonary Wedge Artery flow
        else
                dvq[15]=0
        end

        dvq[16]=q[8]-q[9] #Pulmonary Wedge capillary volume

        P_0d[15]=Epwc*Zpwc+Spwc*dv[8]+pit #Pulmonary Wedge capillary Pressure

        dvq[17]=(Epwc*Zpwc-Epwv*Zpwv-  Rpwc*q[9]+Spwc*dv[8]-Spwv*dv[9])/ypwc #Pulmonary Wedge capillary Flow

        dvq[18]=q[9]-q[10]  #Pulmonary Wedge vein Volume

        P_0d[17]=Epwv*Zpwv+Spwv*dv[9]+pit#Pulmonary Wedge vein Pressure

        dvq[19]=(Epwv*Zpwv-ela*v[10]-  Rpwv*q[10]+Spwv*dv[9]-Sla*dv[10]-ppc)/ypwv #Pulmonary Wedge vein Flow

        dvq[20]=q[7]+q[10]-q[11] #LA volume

        P_0d[19]=ela*v[10]+Sla*dv[10]+ppc+pit #LA Pressure


        dvq[21]=(ela*v[10]-plv-Rmv*q[11]- bmv*q[11]*abs(q[11])+Sla*dv[10]-Slv*dv[11])/ymv  #Mitral flow

        dvq[22]=q[11]-q[12] #LV volume

        P_0d[21]=plv+Slv*dv[11]+ppc+pit#LV Pressure

        dvq[23]=(plv-v[12]/Caor-Rav*q[12]-(bav*q[12])*abs(q[12])+ Slv*dv[11]-Saor*dv[12]+ ppc + pit)/yav#Aortic Flow

        # Adjust the state of cardiac valve
        if ((Aav==0.0 && q[12]<=0.000000001))
                dvq[23]=0.0
        end

        if ((abs(tcr-Tduration)<=ddt*0.5 || tcr<0.1 || (Amv==0.0 && q[11]<0.000000001)))
                dvq[21]=0.0
        end

        if ((Apv==0.0 && q[4]<=0.00000001))
                dvq[8] = 0.0
        end

        if ((abs(tcr-Tduration)<=ddt*0.5 || tcr<0.1 || (Atv==0.0 && q[3]<=0.000000001)))
                dvq[6] = 0.0
        end

        dvq[24]=q[12] - q[13] #Aorta volume

        P_0d[23]=v[12]/Caor + Saor*dv[12] #Aorta Presuure


        dvq[25]= (v[12]/Caor + Saor*dv[12] - q[13]*Raor - v[13]/Cart - Sart*dv[13])/yaor#Aorta Flow
        dvq[26]=q[13] - q[14] #Artery Volume

        P_0d[25]=v[13]/Cart + Sart*dv[13]#Artery Pressure

        dvq[27]= (v[13]/Cart + Sart*dv[13] - q[14]*Rart - v[14]/Ccap -  Scap*dv[14])/yart #Artery Flow

        dvq[28]=q[14] - q[15] #Capillarya Volume

        P_0d[27]=v[14]/Ccap + Scap*dv[14]#Capillary Pressure

        dvq[29]= (v[14]/Ccap + Scap*dv[14] - q[15]*Rcap - v[1]/Cven - Sven*dv[1])/ycap#Capillary flow

        # -------------------------------------------------------------------------
        #y[25] Nbr:t q[16]
        # dvq[30] = q[17]
        # if nstep==0
        #     Pod=-41.76
        #     'chek'
                        # else
        # Pod=(diff(P_0d[23]))*1000
        # end

        dvq[30] =(aao+(6.37*aod)-q[16])/2.076
        #y[26] Nbr_t:t
        #modeqn

        dvq[31]=(-q[17]+gainR)/6
        # dvq[31]   =[1.05*P_0d[23] + 0.036*1.05*((q[12]-q[13])/0.9]-[0.0018+0.001]*q[17] - q[16])/[0.0018*0.001]
        #modeqn
        dvq[32]=(-q[18]+gainElv)/8
        dvq[33]=(-q[19]+gainErv)/8
        dvq[34]=(-q[20]+gainTs)/2
        dvq[35]=(-q[21]+gainTv)/1.5

end


function main_baro(HR)
        global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
        global odic,jj, result, fesmin

        global odic=[1068.2371, 52.4983, 181.7233, -41.7618, 65.0625, 0, 122.6637, 0,
        67.0272, -0.3118, 135.11, -2.1737, 198.7568, -64.1791, 0, 2.7983, 0.1357, 2.7932, 1.1042, 68.8587, 0, 121.2539, 0, 67.3641, 41.8262, 22.0472, 56.6627, 1.8539, 57.6473, 120, 0.8, 1.4, 0.5, -0.13, 0]

        v = zeros(1, 15)
        dv = zeros(1, 21)
        q = zeros(1, 36)
        dvq = zeros(1, 36)

        result = zeros(4, 101)
        P_0d = zeros(1, 101)
        diffv = zeros(1, 20)
        resultcr = zeros(1, 101)
        rukuk = zeros(4, 110)

        MyResult = zeros(100000, 28)
        MyResult1 = zeros(100000, 37)
        MyResult2 = zeros(100000, 15)
        Ervbaro, Rartbaro, Elvbaro, Tbaro = zeros(1,4)

        odic_new = zeros(1,35)
        print(typeof(HR))
        hr = copy(HR)

        Pit=-2.5
        pit=Pit
        jj=0
        kk=0

        #-----------------------------------------
        # Elva=2.87      #!Peak-systolic elastance of left ventricle
        Elvb=0.06      #!Basic diastolic elastance of left ventricle
        Elaa=0.07      #!Peak-systolic elastance of left atrium
        Elab=0.075     #!Basic diastolic elastance of left atrium
        # Erva=0.52      #!Peak-systolic elastance of right ventricle
        Ervb=0.043     #!Basic diastolic elastance of right ventricle
        Eraa=0.055     #!Peak-systolic elastance of right atrium
        Erab=0.06      #!Basic diastolic elastance of right atrium

        Vmax=900.0#Reference volume of Frank-Starling law
        Es=45.9     #!Effective septal elastance
        Vpc0=380.0  #!Reference total pericardial and cardiac volume
        Vpe=30.0    #!Pericardial volume of heart
        Vcon=40.0   #!Volume constant
        Sva0=0.0005 #!Coefficient of cardiac viscoelasticity

        # ! Cardiac valve parameters
        #!(aortic valve(AV),mitral valve(MV), tricuspid valve(TV),pulmonary valve(PV))
        bav=0.000025 #!Bernoulli's resistance of AV
        bmv=0.000016 #!Bernoulli's resistance of MV
        btv=0.000016#!Bernoulli's resistance of TV
        bpv=0.000025 #!Bernoulli's resistance of PV
        Rav=0.005    #!Viscous resistance of AV
        Rmv=0.005    #!Viscous resistance of MV
        Rtv=0.005    #!Viscous resistance of TV
        Rpv=0.005   # !Viscous resistance of PV
        yav=0.0005   #!Inertance of AV
        ymv=0.0002   #!Inertance of MV
        ytv=0.0002   #!Inertance of TV
        ypv=0.0005   #!Inertance of PV

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

        #! Peripheral circulation
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
        # Rart=0.75
        Rcap=0.35
        Rven=0.07
        Rvca=0.001

        Saor=0.01
        Sart=0.01
        Scap=0.01
        Sven=0.01
        Svca=0.01
        qco = 0.0
        
        #---------------------------------------
        # Tduration=input('Please specify cardiac duration(s)')
        # dt= input('Please specify time step(s)')
        # tee=0.3*sqrt(Tduration) #!Moment when ventricular contractility reaches the peak
        # tac=Tduration - 0.5*tee -  0.02* (Tduration/0.855)# !Moment when atrium begins to contract
        # tar=Tduration - 0.02* (Tduration/0.855)# !Moment when atrium begins to relax
        
        # ncycle=input('how many cardiac cycles to run ?')
        # ntotal =5
        ddt = 0.001

        #---------------------------------------------
        for cycle in 1:1:65
                global  Aav, Amv, Apv, Atv, Gpw #Valve Parameters
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
                global odic,jj, result, fesmin
                if cycle == 1
                        Tduration=0.855
                        Rart=0.75
                        Elva=2.87
                        Erva=0.52
                else
                        Tduration=Tbaro
                        Rart=Rartbaro
                        Elva=Elvbaro
                        Erva=Ervbaro
                end
                for nstep in 0:ddt:Tduration
                        if nstep == 0
                                tcr = 0
                                ppc = 0
                                result = zeros(4,101)
                                result[1,1:35] = odic
                                for i in 1:14
                                        if i<8
                                                v[i] = result[1,2*i-1]
                                        else
                                                v[i] = result[1, 2*i]
                                        end
                                end
                                for i in 1:15
                                        if i<8
                                                q[i] = result[1,2*i]
                                        else
                                                q[i] = result[1,2*i-1]
                                        end
                                end
                                q[16:21] = odic[30:35]
                                STR=1.0
                                FL=1.0
                                FR1=1.0
                                Gpw=0.0

                                Aav=0.0
                                Amv=0.0
                                Apv=0.0
                                Atv=0.0

                                Pit=-2.5
                                P_0d=zeros(1,101)

                                aod=-41.7
                                aao=73
                        end
                        tcr = nstep
                        tee=0.3*sqrt(Tduration)
                        tac=Tduration - 0.5*tee -  0.02* (Tduration/0.855)
                        tar=Tduration - 0.02* (Tduration/0.855)

                        ncount=1
                        ncountadd=ncount+1
                        resultcr[1:101]= result[1,1:101]

                        Epua=Ecal(Epua0,Zpua,v[5])
                        Epuc=Ecal(Epuc0,Zpuc,v[6])
                        Epuv=Ecal(Epuv0,Zpuv,v[7])
                        Epwc=Ecal(Epwc0,Zpwc,v[8])
                        Epwv=Ecal(Epwv0,Zpwv,v[9])

                        if tcr == 0
                                FL = 1.0-(result[1,22]/Vmax)
                                FR1=1.0-(result[1,7]/Vmax)
                        end
                        Lvecal()
                        Laecal()
                        rvecal()
                        raecal()

                        cklr=erv/(Es+erv)
                        ckrl=elv/(Es+elv)
                        plv=(ckrl*Es*v[11]+ckrl*cklr*Es*v[4])/(1.0-cklr)
                        prv=(cklr*Es*v[4]+ckrl*cklr*Es*v[11])/(1.0-ckrl)
                        Sla=Sva0*v[10]*ela
                        Slv=Sva0*plv
                        Sra=Sva0*v[10]*era
                        Srv=Sva0*prv
                        ppp=(v[3]+v[4]+v[10]+v[11]+Vpe-Vpc0)/Vcon
                        ppc=exp(ppp)
                        fmin=2.52
                        fmax=47.78
                        Pn=92
                        ka=11.75
                        fcs1=(fmin+fmax*exp((q[16]-Pn)/ka))/(1+exp((q[16]-Pn)/ka))
                        fcs=round(fcs1,digits=3)

                        kes=0.0675
                        fes0=16.11
                        fes1=2.1
                        es1=(-kes)*fcs
                        es2=exp(es1)
                        es3=(fes0-fes1)
                        fes1=fes1+(es3*es2)
                        fes=round(fes1,digits=3)

                        fev0=3.2
                        fev1=6.3
                        kev=7.06
                        fcs0=25
                        fev1=(fev0+fev1*exp((q[16]-fcs0)/kev))/(1+exp((q[16]-fcs0)/kev))
                        fev=round(fev1,digits=3)

                        ############################ ARTERIAL RESISTANCE##########################
                        if jj>2000
                                fesmin = 2.71
                                if fes>fesmin
                                        gainR1=fes-fesmin+1
                                        gainR2=log(gainR1)
                                        GR=0.135
                                        gainR=GR*gainR2
                                else
                                        gainR=0
                                end
                        else
                                gainR=0
                        end
                        ############################ LV elastance ##########################
                        fesmin1 = 2.66
                        Gelv=0.5
                        if jj>0
                                if fes>fesmin
                                        gainElv1=fes-fesmin1+1
                                        gainElv2=log(gainElv1)
                                        gainElv=Gelv*gainElv2
                                else
                                        gainElv=0
                                end
                        else
                                gainElv=0
                        end
                        ############################ RV elastance #######################################
                        fesmin1=2.66
                        Gerv=0.09
                        if jj>2000
                                if fes>fesmin
                                        gainErv1=fes-fesmin1+1
                                        gainErv2=log(gainErv1)
                                        gainErv=Gerv*gainErv2
                                else
                                        gainErv=0
                                end
                        else
                                        gainErv=0
                        end

                        ############################ Cardiac cycle sympathetic ##########################
                        fesmin=2.66
                        GTs=-0.13
                        if jj>2000
                                if fes>fmin
                                        gainTs1=fes-fesmin1+1
                                        gainTs2=log(gainTs1)
                                        gainTs=GTs*gainTs2
                                else
                                        gainTs=0
                                end
                        else
                                gainTs=0
                        end
                        ############################ Cardiac cycle-vaga ##########################
                        Gtv = 0.09
                        if jj>200
                                gainTv=Gtv*fev
                        else
                                gainTv=0
                        end
                        Ervbaro1=q[19]+0.3956
                        Ervbaro=round(Ervbaro1, digits=3)
                        Tbaro1=q[20]+q[21]+0.58
                        Tbaro=round(Tbaro1,digits=3)
                        Elvbaro1=q[18]+2.0829
                        Elvbaro=round(Elvbaro1,digits=3)
                        Rartbaro1=0.534+q[17]
                        Rartbaro=round(Rartbaro1,digits=3)

                        dvdqsdvsdq1()
                        diffv[1:7]=dvq[1:2:13]
                        diffv[8:14] = dvq[16:2:28]
                        dv[1:14]=diffv[1:14]

                        rukuk = rungekuttabaro(resultcr)

                        for i=1:35
                                result[ncountadd, i] = result[ncount,i]+(rukuk[1,i]+2.0*(rukuk[2,i]+rukuk[3,i])+(rukuk[4,i]))/6.0
                        end

                        delta = 0.00000001

                        if (Aav==0.0 && resultcr[23]<=delta)
                                resultcr[23]=0.0
                                result[ncountadd,23]=0.0
                                q[12]=0.0
                        end

                        if (abs(tcr-Tduration)<=ddt*0.5 || tcr<0.1 || (Amv==0.0 && resultcr[21]<=delta))
                                resultcr[21]=0.0
                                result[ncountadd,21]=0.0
                                q[11]=0.0
                        end

                        if (Apv==0.0 && resultcr[8]<=delta)
                                resultcr[8]=0.0
                                result[ncountadd,8]=0.0
                                q[4]=0.0
                        end

                        if (abs(tcr-Tduration)<=ddt*0.5 || tcr<0.1 || (Atv==0.0 && resultcr[6]<=delta))
                                resultcr[6]=0.0
                                result[ncountadd,6]=0.0
                                q[3]=0.0
                        end

                        v[1:7]=result[ncountadd,1:2:13]
                        v[8:14]=result[ncountadd,16:2:28]
                        q[1:7]=result[ncountadd,2:2:14]
                        q[8:15]=result[ncountadd,15:2:29]
                        q[16:21]=result[ncountadd,30:1:35]
                        result[1,1:35]=result[2,1:35]

                        jj = jj+1

                        MyResult[jj,1:27]=P_0d[1:27]
                        MyResult1[jj,1]=fcs
                        MyResult1[jj,2:36]=result[1,1:35]
                        MyResult2[jj,1]=fcs
                        MyResult2[jj,2]=fes
                        MyResult2[jj,3]=fev
                        MyResult2[jj,5]=gainElv
                        MyResult2[jj,6]=elv

                        MyResult2[jj,13]=nstep
                        MyResult2[jj,14]=tcr

                        if jj>1
                                a = jj-1
                                a0d = (MyResult[jj,23]-MyResult[a,23])*1000
                        end

                        MyResult2[jj,4]=aod

                        if nstep==Tduration
                                odic_new = result[1,1:35]
                        end
                end
                for i in size(odic_new)
                        odic[i] = round(odic_new[i], digits=3)
                end
        end

        return MyResult1[16000:25999, 24]
end

#@time result = main_baro()
#println(size(result))
