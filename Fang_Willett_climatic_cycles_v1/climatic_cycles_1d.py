# This py.file is the main codes of 1d model
# Author @X. Fang Aug. 2026
##
import numpy as np
from scipy.optimize import minimize_scalar
# Function nuclide_calculation computes the present-day surface 10Be concentration for a given erosion rate history
def nuclide_calculation(tc,erosion_rates,start_mear_index,P1,P2,P3,attL1,attL2,attL3,density,L_10Be):
    timeStep = tc[1] - tc[0]
    # Spallation, slow muons, fast muons respectively
    L1, L2, L3 = attL1 / density, attL2 / density, attL3 / density

    # surface concentration at time t for Spallation, slow muons, fast muons respectively
    Css_final1, Css_final2, Css_final3 = np.zeros(tc.shape), np.zeros(tc.shape), np.zeros(tc.shape)

    # concentration start from 0
    # erosion depth for every time step
    ero_h = erosion_rates * timeStep*100  # cm
    # erosion depth for the whole history
    ero_h_acc = np.cumsum(ero_h)  # cm

    for i_index in range(start_mear_index, len(tc) - 1, 1):
        ero_h_acc_tem = ero_h_acc[:i_index + 1]
        ero_h_acc_tem = np.insert(ero_h_acc_tem, 0, 0)
        # erosion depth
        # depth
        ero_h_acc_tem = ero_h_acc_tem[-1] - ero_h_acc_tem
        ero_h_acc_tem = ero_h_acc_tem[:-1]
        #
        tem_coff = L_10Be * timeStep - 1
        powers_of_temcoff = tem_coff ** np.arange(i_index + 1)
        signs = (-1) ** np.arange(i_index + 1)

        temP1 = P1 * np.exp(-ero_h_acc_tem / L1) * timeStep
        temP2 = P2 * np.exp(-ero_h_acc_tem / L2) * timeStep
        temP3 = P3 * np.exp(-ero_h_acc_tem / L3) * timeStep

        temP1_rever = np.flip(temP1)
        temP2_rever = np.flip(temP2)
        temP3_rever = np.flip(temP3)

        Css_final1[i_index + 1] = np.dot(signs * powers_of_temcoff, temP1_rever)
        Css_final2[i_index + 1] = np.dot(signs * powers_of_temcoff, temP2_rever)
        Css_final3[i_index + 1] = np.dot(signs * powers_of_temcoff, temP3_rever)
    Css_total = Css_final1 + Css_final2 + Css_final3
    return Css_total

# Function apparent_erosion_rates computes the Ecos when the cosmogenic nuclide contentration is known
def apparent_erosion_rates(P1,P2,P3,attL1,attL2,attL3,density,L_10Be,Css_total):
    L1, L2, L3 = attL1 / density, attL2 / density, attL3 / density
    apparentE=np.full(len(Css_total), np.nan)

    for i_index in range(len(Css_total)):
        print(i_index,'/',len(Css_total))
        lower_limit = 1e-8 # cm/y
        upper_limit = 1e1 # cm/y

        low_lim = P1 * L1 / (lower_limit + L1 * L_10Be) + P2 * L2 / (
                    lower_limit + L2 * L_10Be) + P3 * L3 / (lower_limit + L3 * L_10Be) - Css_total[i_index]
        up_lim  = P1 * L1 / (upper_limit + L1 * L_10Be) + P2 * L2 / (
                    upper_limit + L2 * L_10Be) + P3 * L3 / (upper_limit + L3 * L_10Be) - Css_total[i_index]

        # make sure upper_limit >Ereal, lower_limit<Ereal
        n_limit = 10
        n1 = 1
        while (up_lim > 0 and n1 < n_limit):
            upper_limit = upper_limit * 10
            up_lim = P1 * L1 / (upper_limit + L1 * L_10Be) + P2 * L2 / (
                    upper_limit + L2 * L_10Be) + P3 * L3 / (upper_limit + L3 * L_10Be) - Css_total[i_index]
            n1 = n1 + 1

        n2 = 1
        while (low_lim < 0 and n2 < n_limit):
            lower_limit = lower_limit / 10
            low_lim = P1 * L1 / (lower_limit + L1 * L_10Be) + P2 * L2 / (
                    lower_limit + L2 * L_10Be) + P3 * L3 / (lower_limit + L3 * L_10Be) - Css_total[i_index]
            n2 = n2 + 1

        #
        n = 50
        i = 0
        a = 1e8

        if (np.sign(low_lim) != np.sign(up_lim)):
            while (abs(a) > 1e-3 and i < n):
                i = i + 1
                mid = 10 **((np.log10(lower_limit) + np.log10(upper_limit)) / 2)
                if i > 1:
                    if ((mid - upper_limit) == 0):
                        break
                    elif((mid - lower_limit) == 0):
                        break

                a =P1 * L1 / (mid + L1 * L_10Be) + P2 * L2 / (
                    mid + L2 * L_10Be) + P3 * L3 / (mid + L3 * L_10Be) - Css_total[i_index]

                if (np.sign(low_lim) == np.sign(a)):
                    lower_limit = mid
                else:
                    upper_limit = mid

                if (a == 0):
                    break

        apparentE[i_index]=mid
    apparentE=apparentE*1e-2 # change from cm/y to m/y
    return apparentE
##
# ###############################################################
# ###############################################################
# Section One: Paramter setting
# 1 catchment input setting
# uplift: m/yr; Pm: m/yr; Pa: m/yr
Point1=[0.5e-3,0.4,0.2]
Point2=[0.5e-3,0.6,0.3]

Point3=[0.5e-3,0.4,0.1]
Point4=[0.5e-3,1.2,0.3]

Point5=[0.5e-3,0.4,0.05]
Point6=[0.5e-3,2.4,0.3]

PointList=[Point1,Point2,Point3,Point4,Point5,Point6]


# 2  1d river profile setting
# river length
L=10e3 # m
disStep=25 # m

# precipitation parameters
timeTotal=10000e3 # 10 Myr unit:yr

timeStep=1e2 # unit: yr
Period_prec = 20e3 # 20kyr unit: yr
timecycle = timeTotal / Period_prec  # cycle number

# Stream power law and Hack's law
m=0.5
n=1
Kp=1e-6 # m-1 * y-1
ka=6.69
h=1.67
Ac=1e6 # m2

# cosmogenic parameters
rock_density=2.65 # g/cm3
Half_10Be=1.387e6 # y
L_10Be=np.log(2)/Half_10Be # namuda
attL_spal_10Be=160 # g/cm2 Spallation attenuation length (Braucher et al., 2011) [g/cm2]
attL_ms_10Be = 1500  # Slow muons attenuation length (Braucher et al., 2011) [g/cm2]
attL_mf_10Be = 4320 # Fast muons attenuation length (Braucher et al., 2011) [g/cm2]
P_spal_SLHL=4 # at/g/yr
P_ms_SLHL = 0.012 # at/g/yr
P_mf_SLHL = 0.039 # at/g/yr

A0,P0=1,1
t=np.linspace(0,timeTotal-timeStep,int(timeTotal/timeStep))
x=np.linspace(0,L,int(L/disStep)+1) # upstream distance
A=ka*np.power(L-x,h)+Ac
A=np.append(A,0)
# Ein, Ecos, ks, ksp,kspm
Ein_tList,Ecos_tList,Ks_tList,Ksp_tList,Kspm_tList=[],[],[],[],[]
# Ecos calcualtion
start_cos=int(len(t)-1-100*(Period_prec/timeStep)) # the last 100 cycles
end_cos=len(t)-1
t_cos=t[start_cos:end_cos]
# precipitation history
Prec_tList=[]

##
# ###############################################################
# ###############################################################
# Section two: 1d model running
# for each catchment
for point_index,point_item in enumerate(PointList):
    U=point_item[0] # U
    Prec_0=point_item[1] # Pm
    Prec_amplitude=point_item[2] # Pa
    print('Point ',point_index,'/',str(len(PointList)))
    # ###############################################################
    # 1 Model running: calculating elevation for every x and every t
    Prec_t = Prec_amplitude * np.sin(2 * np.pi * t / Period_prec) + Prec_0
    Prec_tList.append(Prec_t)
    Z=np.zeros((t.size,x.size))
    coffi=np.power(U/Kp,1/n)

    # 1 calculate Z during time
    # simplified equation
    for j in range(0,len(t)-1,1):
        for i in range(1,len(x),1):
            temE_dz_dx = Kp * np.power(A[i] * Prec_t[j + 1], m)
            dz_dx_before=(np.power(Z[j][i]-Z[j][i-1],n-1))/(np.power(disStep,n))
            Z[j + 1][i]=(U*timeStep+temE_dz_dx*dz_dx_before*timeStep*Z[j+1][i-1]+Z[j][i])/(1+temE_dz_dx*dz_dx_before*timeStep)
    print(str(point_index),'Z calculation finished')
    # ###############################################################
    # 2 calculate
    # E & Ks & Ks_p
    # x-averaged E, Ks, and Ks_p
    # E & Ks & Ks_p
    E,ks,ks_p,chi_p=np.zeros((t.size,x.size)),np.zeros((t.size,x.size)),np.zeros((t.size,x.size)),np.zeros((t.size,x.size))
    ks_pm=np.zeros((t.size,x.size))
    for j in range(start_cos,len(t),1):
        tem_tao_sum=0
        tem_AP=np.power(A * Prec_t[j], m)
        tem_Pmn=np.power(Prec_t[j], m / n)
        tem_APmm = np.power(A0 * P0 / (A * Prec_t[j]), m / n)
        for i in range(1,len(x),1):
            dz_dx=(Z[j][i]-Z[j][i-1])/disStep
            E[j][i]=Kp*tem_AP[i]*np.power(dz_dx,n)
            ks[j][i] = dz_dx * np.power(A[i], m / n)
            ks_p[j][i]= ks[j][i] * tem_Pmn
            ks_pm[j][i]=ks[j][i]*np.power(Prec_0, m / n)
            chi_p[j][i] = chi_p[j][i - 1] + ((tem_APmm[i - 1] + tem_APmm[i]) / 2) * disStep
    # calculate x-averaged erosion rates, Ks, and Ks_p
    E_avr=np.zeros(t.size)
    Ks_p_avr=np.zeros(t.size)
    Ks_avr=np.zeros(t.size)
    Ks_pm_avr=np.zeros(t.size)
    for j in range(start_cos,len(t),1):
        Esum,Ksnsum,Ksnqsum,Ksnqmsum=0,0,0,0
        for i in range(1,len(x),1):
            temSlope=(Z[j][i]-Z[j][i-1])/disStep
            Esum += Kp * np.power(Prec_t[j]*A[i], m) * np.power(temSlope, n) * (A[i] - A[i + 1])
            Ksnqsum += temSlope * np.power(A[i] * Prec_t[j], m / n) * (A[i] - A[i + 1])
            Ksnsum += temSlope * np.power(A[i], m / n) * (A[i] - A[i + 1])
            Ksnqmsum+=temSlope * np.power(A[i] * Prec_0, m / n) * (A[i] - A[i + 1])
        temEavr=Esum/A[0]
        temKsnqavr=Ksnqsum/A[0]
        temKsnavr=Ksnsum/A[0]
        temKsnqmavr=Ksnqmsum/A[0]
        E_avr[j]=temEavr
        Ks_p_avr[j]=temKsnqavr
        Ks_avr[j]=temKsnavr
        Ks_pm_avr[j]=temKsnqmavr
    Ein_tList.append(E_avr)
    Ks_tList.append(Ks_avr)
    Ksp_tList.append(Ks_p_avr)
    Kspm_tList.append(Ks_pm_avr)
    print(str(point_index),'river steepness finished')
##
# ###############################################################
# ###############################################################
# Section three: Ecos and Ecos-qm calculation
# 1 Ecos calculation
t_start_cos = t[start_cos]
t_end_cos=t[end_cos]

cos_timeStep=1e1
ratio_timeStep=int(timeStep/cos_timeStep)
# time period for cosmogenic calculation
t_cos = np.arange(t_start_cos, t_end_cos + timeStep, cos_timeStep)
# time for cosmogenic measurement the last 5 cycles
start_mear_index=int(len(t_cos)-5*(Period_prec/timeStep)*ratio_timeStep)
Ein_t_repList,Ecos_tList,C10Be_tList=[]

for point_index in range(0,len(Prec_tList),1):
    print('Point ', point_index, '/', str(len(Prec_tList)))
    temEin_t = Ein_tList[point_index][start_cos:end_cos+1]
    Ein_rep =np.array([num for num in temEin_t for _ in range(ratio_timeStep)])
    Ein_t_repList.append(Ein_rep)
    C10Be = nuclide_calculation(t_cos, Ein_rep, start_mear_index-1,P_spal_SLHL,P_mf_SLHL,P_ms_SLHL,attL_spal_10Be, attL_mf_10Be,attL_ms_10Be,rock_density,L_10Be)
    E_cos = apparent_erosion_rates(P_spal_SLHL, P_mf_SLHL, P_ms_SLHL, attL_spal_10Be,attL_mf_10Be,attL_ms_10Be, rock_density, L_10Be, C10Be[start_mear_index:])
    C10Be_tList.append(C10Be[start_mear_index:])
    Ecos_tList.append(E_cos)
    print('Point ', point_index, 'cos erosion finished')

startp = int(len(t)-5*(Period_prec/timeStep))
endp=len(t)-1

Ecos_tList_real=[ np.zeros(t.size) for i in range(len(Prec_tList))]
C10Be_tList_real=[ np.zeros(t.size) for i in range(len(Prec_tList))]

t_index=np.linspace(start=startp, stop=endp, num=endp-startp+1, dtype=int)
tem_cos_index=(t_index-start_cos)*ratio_timeStep-start_mear_index
for i in range(len(Prec_tList)):
    Ecos_tList_real[i][t_index]=Ecos_tList[i][tem_cos_index]
    C10Be_tList_real[i][t_index] = C10Be_tList[i][tem_cos_index]

 # ###############################################################
# 2 Ecos-qm calculation
Ecosmodi_tList=[np.zeros(t.size) for i in range(len(Prec_tList))]
bestc_tList=[np.zeros(t.size) for i in range(len(Prec_tList))]
bestb_tList=[np.zeros(t.size) for i in range(len(Prec_tList))]

for point_index in range(0,len(Prec_tList),1):
    for j in range(startp,len(t)-1,1):
        temEcos = Ecos_tList_real[point_index][j]
        temC10Be=C10Be_tList_real[point_index][j]
        temt=t[int(j-(Period_prec/timeStep)*5):j+1]
        temtemP=Prec_tList[point_index][int(j-(Period_prec/timeStep)*5):j+1]
        temtemEin=np.power(temtemP,0.5)/1e3 # change unit to cm/yr

        def error_func(c):
            temEin = temtemEin * c
            C10Be_last = nuclide_calculation(temt, temEin, len(temEin) - 2,
                                             P_spal_SLHL, P_mf_SLHL, P_ms_SLHL,
                                             attL_spal_10Be, attL_mf_10Be, attL_ms_10Be,
                                             rock_density, L_10Be)
            return abs(C10Be_last[-1] - temC10Be)

        res = minimize_scalar(error_func, bounds=(1e-4, 10), method='bounded')
        best_c = res.x
        best_error = res.fun
        temEcosmodi=np.mean(temtemEin*best_c)
        Ecosmodi_tList[point_index][j] = temEcosmodi
        bestc_tList[point_index][j] = best_c
        bestb_tList[point_index][j] = best_error


## Figure setting
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter

startp=int((Period_prec / timeStep) * (timecycle - 1))
endp=len(t)-1
lenp=int(endp+1-startp)

chooseinterval=1/20
chooseinterval_index=int(chooseinterval*(Period_prec/timeStep))
chooseIndex=[startp+temindex*chooseinterval_index for temindex in range(int(1/chooseinterval))]
chooseIndex=np.array(chooseIndex)
# color
# gradual color for markers
cmap_Reds = plt.get_cmap('Reds')
cmap_Reds_color=[cmap_Reds(i) for i in np.linspace(0, 1, len(chooseIndex)+1)]
cmap_Reds_color=cmap_Reds_color[1:]

cmap_Wistia = plt.get_cmap('Wistia')
cmap_Wistia_color=[cmap_Wistia(i) for i in np.linspace(0, 1, len(chooseIndex)+1)]
cmap_Wistia_color=cmap_Wistia_color[1:]

cmap_PuBu = plt.get_cmap('Blues')
cmap_PuBu_color=[cmap_PuBu(i) for i in np.linspace(0, 1, len(chooseIndex)+1)]
cmap_PuBu_color=cmap_PuBu_color[1:]

mycolor0 = [cmap_Reds_color[-10],cmap_Reds_color[-2],cmap_Wistia_color[-10],cmap_Wistia_color[-10],cmap_PuBu_color[-10],cmap_PuBu_color[-10]]  # red, orange, blue
mycolor=[cmap_Reds_color,cmap_Reds_color,cmap_Wistia_color,cmap_Wistia_color,cmap_PuBu_color,cmap_PuBu_color]

## Figure 1-1 1cycle, precipiation, Ein and Ecos, Ks, ksp, kspm
chooseinterval=1/4
chooseinterval_index=int(chooseinterval*(Period_prec/timeStep))
chooseIndex=[startp+temindex*chooseinterval_index for temindex in range(int(1/chooseinterval))]
chooseIndex=np.array(chooseIndex)
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
mymarker=['o','*']
mylinewidth=[2,1,0.5]
mymarkersize=[20,20,20]
plt.figure(figsize=(16,4))
# 1 precipitation
plt.subplot(1,3,1)
for point_index in range(0,len(Prec_tList),1):
    plt.plot(t[startp:endp+1]/ 1e3, Prec_tList[point_index][startp:endp+1], c=mycolor0[point_index],linewidth=mylinewidth[0],zorder=1)
    plt.scatter(t[chooseIndex] / 1e3, Prec_tList[point_index][chooseIndex], 20,c=mycolor0[point_index][:3],marker=mymarker[int(point_index%2)],edgecolors='k',linewidths=0.2)
plt.xticks([9980,9985,9990,9995,10000],[0,5,10,15,20])
plt.ylim(0,2.8)
plt.yticks([0,0.5,1,1.5,2,2.5])
# 2 Ein, Ecos, Ecos-qm
plt.subplot(1,3,2) # Ecos and Ein
for point_index in range(0,len(Prec_tList),1):
    plt.plot(t[startp:endp+1] / 1e3, Ein_tList[point_index][startp:endp+1]*1e3, c=mycolor0[point_index], linewidth=mylinewidth[0], linestyle='solid',zorder=1)
    plt.scatter(t[chooseIndex] / 1e3, Ein_tList[point_index][chooseIndex] * 1e3, mymarkersize[0],c=mycolor0[point_index][:3],marker=mymarker[int(point_index%2)],edgecolors='k',linewidths=0.2)
    plt.plot(t[startp:endp+1] / 1e3, Ecos_tList_real[point_index][startp:endp+1]*1e3, c=mycolor0[point_index],linewidth=mylinewidth[1], linestyle='dashed',zorder=1)
    plt.scatter(t[chooseIndex] / 1e3, Ecos_tList_real[point_index][chooseIndex] * 1e3,mymarkersize[1], c=mycolor0[point_index][:3],
                marker=mymarker[int(point_index % 2)],edgecolors='k',linewidths=0.2)
    plt.plot(t[startp:endp] / 1e3, Ecosmodi_tList[point_index][startp:endp]*1e3, c=mycolor0[point_index],
             linewidth=mylinewidth[2], linestyle='dotted', zorder=1)
    plt.scatter(t[chooseIndex] / 1e3, Ecosmodi_tList[point_index][chooseIndex] * 1e3,mymarkersize[2], c=mycolor0[point_index][:3],
                marker=mymarker[int(point_index % 2)],edgecolors='k',linewidths=0.2)

plt.xticks([9980,9985,9990,9995,10000],[0,5,10,15,20])
plt.ylim(0.3,0.7)
plt.yticks([0.3,0.4,0.5,0.6,0.7])
# 3 ks,ksp,kspm
plt.subplot(1,3,3)
for point_index in range(0,len(Prec_tList),1):
    #plt.plot(t[startp:endp + 1] / 1e3, Ks_tList[point_index][startp:endp + 1], c=mycolor0[point_index], zorder=1, linestyle='dashed',linewidth=mylinewidth[1])
    plt.plot(t[startp:endp + 1] / 1e3, Ksp_tList[point_index][startp:endp + 1], c=mycolor0[point_index], zorder=1,
             linestyle='solid',linewidth=mylinewidth[0])
    plt.scatter(t[chooseIndex] / 1e3, Ksp_tList[point_index][chooseIndex], mymarkersize[1],
                c=mycolor0[point_index][:3],
                marker=mymarker[int(point_index % 2)], edgecolors='k', linewidths=0.2)
    plt.plot(t[startp:endp + 1] / 1e3, Kspm_tList[point_index][startp:endp + 1], c=mycolor0[point_index], zorder=1,
             linestyle='dotted',linewidth=mylinewidth[2])
    plt.scatter(t[chooseIndex] / 1e3, Kspm_tList[point_index][chooseIndex], mymarkersize[1],
                c=mycolor0[point_index][:3],
                marker=mymarker[int(point_index % 2)], edgecolors='k', linewidths=0.2)
plt.xticks([9980,9985,9990,9995,10000],[0,5,10,15,20])
plt.ylim(100,900)
plt.yticks([100,300,500,700,900])
plt.tight_layout()
plt.savefig('F 1-1 PaPm PEKs.jpg',dpi=400,bbox_inches='tight')
plt.show()
## Figure 1-2 Ecos vs ksq, Ecos-qm vs ksom(ksqeff)
chooseinterval=1/100
chooseinterval_index=int(chooseinterval*(Period_prec/timeStep))
chooseIndex=[startp+temindex*chooseinterval_index for temindex in range(int(1/chooseinterval))]
chooseIndex=np.array(chooseIndex)
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
# Ecos vs ksq
for point_index in range(0,len(Prec_tList),1):
    axs[0].scatter(Ecos_tList_real[point_index][chooseIndex] * 1e3, Ksp_tList[point_index][chooseIndex], c=mycolor[point_index],edgecolors=mycolor[point_index][-10],s=15,linewidths=0.05,marker=mymarker[int(point_index%2)],zorder=2, alpha=1)

axs[0].set_ylim(300,700)
axs[0].set_xlim(0.3,0.7)
axs[0].set_xticks([0.3,0.4,0.5,0.6,0.7])
axs[0].set_yticks([300,400,500,600,700])
axs[0].xaxis.set_major_formatter(ScalarFormatter())
axs[0].yaxis.set_major_formatter(ScalarFormatter())

# Ecos-qm vs ksom(ksqeff)
for point_index in range(0,len(Prec_tList),1):
    axs[1].scatter(Ecosmodi_tList[point_index][chooseIndex] * 1e3, Kspm_tList[point_index][chooseIndex], c=mycolor[point_index],edgecolors=mycolor[point_index],s=15,linewidths=0.5,marker=mymarker[int(point_index%2)],zorder=2, alpha=1)

axs[1].set_ylim(300,700)
axs[1].set_xlim(0.3,0.7)
axs[1].set_xticks([0.3,0.4,0.5,0.6,0.7])
axs[1].set_yticks([300,400,500,600,700])
axs[1].xaxis.set_major_formatter(ScalarFormatter())
axs[1].yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.savefig('F 1-2 PaPm ks vs E.jpg',dpi=400,bbox_inches='tight')
plt.show()
