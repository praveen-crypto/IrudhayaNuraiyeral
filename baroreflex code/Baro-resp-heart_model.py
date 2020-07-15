import numpy as np

odic = np.array([80,3565,3300,-80,-5.6,3565,134.08,125.95,106.26,13.21,19.967,36.26,96.1,34.92,19.486,36.499,66.54,67,1068.2371,52.4983,181.7233,-41.7618,65.0625,0,122.6637,0,67.0272,-0.3118,135.11,-2.1737,198.7568,-64.1791,0,2.7983,0.1357,2.7932,1.1042,68.8587,0,121.2539,0,67.3641,41.8262,22.0472,56.6627,1.8539,57.6473])

dv = np.zeros(shape=(2,21))

# Pit = -2.5

jj = 0
kk = 0

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


FRC =  3.9
VT =  0.5
RV =  1650
TLC =  5550
VD =  185
RespR =  10
Ac =  7.09/1.36
Al = 0.2/1.36
As =  2.2/(1.36*1000)
Au =  0.34/(1.36*1000)
Bc =  37.3/1.36
Bcp =  3.73/1.36
Bl =  -0.5/1.36
Bs =  0.02/(1.36*1000)
Cve =  0.5*(1.36*1000)
Kc_air =  0.21/(1.36*1000)
Kl =  1/1000
Ks =  -10.9
Ku =  0.46/(1.36*10)
Ku = pow(Ku, 6)
Rve =  1/(1.36*1000)
Vstar =  5.3*1000
Vcmax =  0.185*1000
A = -2.17
B = 5.6
f = 15
tau = 0.1
f1 = f/60


####################Gas Exchange##############
Tbody = 300 
Pstp = 760  
Tstp = 273  
nH = 3.5  
P50_O2 = 26.5  
CHb = 0.0000204  
Hcrit = 0.45  
alphaO2 = 1.36E-09  
Vpcmax = 0.07125*1000  
Visf = 48000  
PS1 = 166.67  
Vcytox = 0.21 
Kcytox = 1.00E-12  
P50_CO2 = 250  
alphaCO2 = 3.26E-08 
RQ = 0.8  
Pao = 760  
PH2O = 47  
tweight = 60
r_Pao_O2 = 2.10E-01
r_Pao_CO2 = 0.0003


#  CO2flux = 0
#  O2flux = 0
 
#---------------------------------------
# Tduration = input('Please specify cardiac duration(s)')
dt =  float(input('Please specify time step(s)'))
# tee = 0.3*sqrt(Tduration) #!Moment when ventricular contractility reaches the peak
# tac = Tduration - 0.5*tee -  0.02* (Tduration/0.855)# !Moment when atrium begins to contract
# tar = Tduration - 0.02* (Tduration/0.855)# !Moment when atrium begins to relax
ddt =  dt
# ncycle = input('how many cardiac cycles to run ?')
ntotal = 45


# Tduration = 0.855#input('Please specify cardiac duration(s)')
# dt = 0.0001#input('Please specify time step(s)')

# ddt =  dt
# ncycle = 20#input('how many cardiac cycles to run ?')
# ntotal = ((ncycle)*Tduration)

for nstep in np.arange(0,ntotal,ddt):
    if nstep == 0:

        #tcr = 0
        ppc = 0.0
        result = np.zeros(shape=(2,101))
        result[1, :47] = odic[i]

        for i in np.range(18):
            #if i<6:
                q[0, i] = result[0, i+1]
            #else:
            #    q[i] = result[0, 2*i]
        q[0, 18:25:1] = result[0, 19:32:2]
        q[0, 25:33:1] = result[0, 32:47:2]

        v[0, 0:7:1] = result[0, 18:31:2]
        q[0, 18:25:1] = result[0, 19:32:1]
        STR = 1.0
        FL=1.0
        FR1=1.0
        Gpw=0.0

        Aav=0.0
        Amv=0.0
        Apv=0.0
        Atv=0.0

        Pit=0
        P_0d=np.zeros(1,101)


    #--------------------------------------
    #%c*-----------------------Start computation----------------------*c
    print(nstep)
    if nstep>=0 and nstep<=0.855:
        Tduration = 0.855
        tcr=np.rem(nstep,Tduration)
        tee=0.3*np.sqrt(Tduration);                         #!Moment when ventricular contractility reaches the peak
        tac=Tduration - 0.5*tee -  0.02* (Tduration/0.855); # !Moment when atrium begins to contract
        tar=Tduration - 0.02* (Tduration/0.855);            # !Moment when atrium begins to relax
        
    elif aa1==0:
        gbaro=0.005
        Tduration= gbaro +0.855;                            # *(P1-130)+0.855
        tcr=np.rem(nstep,Tduration)
        tee=0.3*np.sqrt(Tduration);                         # !Moment when ventricular contractility reaches the peak
        tac=Tduration - 0.5*tee -  0.02;                    # *(Tduration/0.855);% !Moment when atrium begins to contract
        tar=Tduration - 0.02;                               #  (Tduration/0.855);% !Moment when atrium begins to relax
        # else 
        # tcr=rem(nstep,Tduration)
    else:
         tcr=np.rem(nstep,Tduration)
 
    aa1=np.roundn(tcr,-3)
    ncount=1
    ncountadd=ncount+1
    resultcr[0:47]= result[0, 0:47]
#     tcr=rem((nstep),Tduration)
    tresp=np.rem((nstep),4)
#     Pit=Ppl
#     Pit = (-1*sin(2*pi*f1*(tresp)) - 1.7)
#     pit=Pit

    a  =  np.roundn(tee,-3)
    b = np.roundn(tcr,-3)
    if a == b:
        P1=P_0d[23]

    Epua = Ecal(Epua0,Zpua,v(5))
    Epuc = Ecal(Epuc0,Zpuc,v(6))
    Epuv = Ecal(Epuv0,Zpuv,v(7))
    Epwc = Ecal(Epwc0,Zpwc,v(8))
    Epwv = Ecal(Epwv0,Zpwv,v(9))
    
    if tcr==0.0:
        FL  = 1.0-(result[1,30] / Vmax)
        FR1 = 1.0-(result[1,25] / Vmax)
    
    Lvecal()
    Laecal()
    rvecal()
    raecal()

    cklr=np.divide(erv, (Es+erv))
    ckrl=np.divide(elv, (Es+elv))

    plv= (ckrl * Es * v(11) + ckrl * cklr * Es * v(4)) / (1.0 - cklr)
    prv= (cklr * Es * v(4) + ckrl * cklr * Es * v(11)) / (1.0 - ckrl)
    Sla= Sva0 * v(10) * ela
    Slv= Sva0 * plv
    Sra= Sva0 * v(3) * era
    Srv= Sva0 * prv
    ppp= (v(3)+v(4)+v(10)+v(11)+Vpe-Vpc0)/Vcon
    ppc= exp(ppp)


Vcw = q(3) + VD
aa = Vcmax/q(1)
bb = 1/aa
Rco = 0
if q(1) > Vcmax:
    Rc = Kc_air + Rco
else:
    Rc = Kc_air * (Vcmax / q(1))^2 + Rco



a1 = q(3) - RV

a2 = Ks * a1
a3 = Vstar-RV
a4 = exp(a2/a3)
Rs = As * a4+Bs
Ru = Au+Ku * np.abs(q(2))

#% b1=TLC-RV;
#% b2=Vcw-RV;
#% b3=b1/b2;
#% b4=b3-0.999;
#% b5=log(b4);
#% Pcw=-Acw+Bcw*b5;
#%  Pcw = (A*sin(2*pi*f1*(tresp)) - B); %Mechanical Ventilation%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Valsalva/Mueller%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pin=30;
Pfinal=-2;
Tcons=0.8;
te=2;
% for Tcons=0.5:0.1:2.5
aa=1-exp(-te/Tcons);
a=(Pin-Pfinal)/aa;
b=Pin-a;
% count=0;
% for t1=0:0.01:60
%     count=count+1;
    if nstep<=8
     Pcw=-2.5;
    else
    nstep1=nstep-8;
if nstep1<=2
    Pin1=-2.5;
    Pfinal1=30;
    Tcons1=0.5;
    te1=2;
% for Tcons=0.5:0.1:2.5
aa1=1-exp(-te1/Tcons1);
a1=(Pin1-Pfinal1)/aa1;
b1=Pin1-a1;
    Pcw=a1*exp(-nstep1/Tcons1)+b1;
elseif nstep1>=12
    Pcw=-2.5;
else
    nstep=nstep1-2;
    PitZ=rem(nstep1,10);
    if PitZ>=0 && PitZ<=8
         Pcw=30;
    else
    PitZ1=PitZ-8;
    Pcw=a*exp(-PitZ1/Tcons)+b;
%     if Tcons==0.5
    end
end








