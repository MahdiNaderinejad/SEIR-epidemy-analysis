import matplotlib.pyplot as plt
import numpy


def f(y,**cte):
    mio=cte['mio']
    beta=cte['beta']
    sigma=cte['sigma']
    delta=cte['delta']
    omega=cte['omega']
    nuo=cte['nuo']
    gama=(1-delta)/14
    alpha=delta/14
    S,E,I,R=y[0],y[1],y[2],y[3]
    N=S+E+I+R
    f_arr=[mio*N-beta*I*S/N+omega*R-mio*S-nuo*S,
           -sigma*E-mio*E+beta*I*S/N,
           sigma*E-gama*I-(mio+alpha)*I,
           gama*I-omega*R-mio*R+nuo*S]
    return f_arr

def k1(y,**cte):
    return f(y,**cte)

def k2(y,**cte):
    y_k1=[yi+k1i*h/2 for yi,k1i in zip(y,k1(y,**cte))]
    return f(y_k1,**cte)

def k3(y,**cte):
    y_k2=[yi+k2i*h/2 for yi,k2i in zip(y,k2(y,**cte))]
    return f(y_k2,**cte)

def k4(y,**cte):
    y_k3=[yi+k3i*h for yi,k3i in zip(y,k3(y,**cte))]
    return f(y_k3,**cte)
    
def yNew(y,**cte):
    y_new=[yi+(k1i+2*k2i+2*k3i+k4i)*h/6 for yi,k1i,k2i,k3i,k4i in zip(y,k1(y,**cte),k2(y,**cte),k3(y,**cte),k4(y,**cte))]
    return y_new
    
days=365
N0=100
E0=0.1/100*N0

mio=1/(76*365)
beta=0.5
sigma=1/7
delta_list=[0,0.02]
omega=1/365
nuo_list=[0,1/365]
h=10

err_stop=0.5*10**(-2)
err=[1,1,1,1]
old_calc=[100,100,100,100]

while(err[0]>err_stop or err[1]>err_stop or err[2]>err_stop or err[3]>err_stop):
    h=h/10
    last_N,S_list,E_list,I_list,R_list,N_list=[],[],[],[],[],[]
    
    for nuo in nuo_list:
        for delta in delta_list:
            if nuo==nuo_list[1] and delta==delta_list[0]:
                continue
            
            y_list=[[N0-E0,E0,0,0]]
            died_of_illness,vaccinated=[],[]
            for day in range(1,int(days/h)):
                y_list.append(yNew(y_list[-1],h=h,mio=mio,beta=beta,sigma=sigma,delta=delta,omega=omega,nuo=nuo))
            
            days_list=h*numpy.array(list(range(len(y_list))))
            S_list.append([y[0] for y in y_list])
            E_list.append([y[1] for y in y_list])
            I_list.append([y[2] for y in y_list])
            R_list.append([y[3] for y in y_list])
            N_list.append([sum(y) for y in y_list])
    
    probability_noVac_notInfec=[(1-beta*S*I/N**2)*(1-nuo_list[1]) for S,I,N in zip(S_list[-1],I_list[-1],N_list[-1])]
    probability=numpy.prod(numpy.array(probability_noVac_notInfec))**(h)
    last_calc=[N_list[0][-1],N_list[1][-1],N_list[2][-1],probability]        
    err=[]
    for new,old in zip(last_calc,old_calc):
        err.append(abs(new-old))

    old_calc=last_calc


ttl=[' - (Death=N, Vac.=N)',' - (Death=Y, Vac.=N)',' - (Death=Y, Vac.=Y)']
for i in range(3):
    plt.subplot(3,1,i+1)
    plt.plot(days_list, S_list[i],color='#303030',label='Susceptible')
    plt.plot(days_list, E_list[i],color='green',label='Exposed')
    plt.plot(days_list, I_list[i],color='red',label='Infectious')
    plt.plot(days_list, R_list[i],color='blue',label='Recovered')
    plt.plot(days_list, N_list[i],color='orange',label='Society')
    plt.ylabel('Groups\' population(%)'+ttl[i])
    plt.xlabel('Days passed')
    if i==0:
        plt.legend(loc=1, bbox_to_anchor=(1.32,1))
plt.subplots_adjust(top=3.1)


plt.savefig('P-810698317-results.png',bbox_inches = 'tight', dpi=400)
res_file = open("P-810698317-results.txt", "w")
res_file.write(f'Note: Calculations are done using steps of length h= {h} day.\n')
res_file.write(f'>>>Ans-A1: Assuming no deaths due to disease, the population is {last_calc[0]:.2f}% after 365 days that means no change in population. (As expected)\n')
res_file.write(f'>>>Ans-A2: About {100-N_list[1][-1]:.2f}% of society has died of this illness in 365 days with no vaccination plan available. The population has dcreased to {N_list[1][-1]:.2f}%.\n')    
res_file.write(f'>>>Ans-B:  About {100-N_list[-1][-1]:.2f}% of society has died of this illness in 365 days while a vaccination procedure was ongoing. The population has decreased to {N_list[2][-1]:.2f}%.\n')  
res_file.write(f'>>>Ans-C:  This plan can vaccinate {sum(numpy.array(S_list[-1])*nuo)*h:.2f}% of the people which means there are {2*sum(numpy.array(S_list[-1])*nuo)*h:.2f} shots for every 100 people while everybody is intended to get two dosages of vaccine.\n')
res_file.write(f'>>>Ans-D:  The probability of not getting exposed to and infected by the diseases while not getting vaccinated for one year(365 days) is equal to {probability*100:.2f}%. ')
res_file.close()


