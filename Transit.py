import numpy as np 
import matplotlib.pyplot as plt 


solar_mass = 1.989*10**30
solar_radius = 6.957*10**8
G = 6.67408*10**(-11)
m_star = 1*solar_mass
jup_radius = 71492000
# R_star = 1*solar_radius
# R_planet = R_star*delta**2
time = np.linspace(0,80,1000)
flux = np.ones(len(time))
# t_T = 40
# t_F = 30
f = []


def light_curve(R_star,time,delta,flux,t_T,i,P,m_star):
    R_planet = (R_star*delta**0.5)/jup_radius
    print('Radius of Planet= %0.3f Jupiter_radius '%R_planet)             # in terms of jupiter mass
    a = ((P**2*G*m_star)/(4*np.pi**2))**(1/3)
    r = a/R_star    
    rat = R_planet/R_star                # ratio of radius of planet with respect to star
    # t_F = t_T*(((1-rat)**2-(r*np.cos(i))**2)/((1+rat)**2-(r*np.cos(i))**2))**(1/2)
    t_F = t_T-10
    t_0 = time[0]                   # initial observing time
    t_1 = t_0 + 30                  # transit starting time
    t_2 = t_1+(t_T-t_F)/2           #flat transit starting time
    t_3 = t_2+t_F                   # flat transit ending time
    t_4 = t_1+t_T                   # transit over time
    t_5 = time[-1]                  #last observing time

    for i in range(len(time)):
        if t_0<=time[i]<=t_1:
            f_1 = flux[i]
            f.append(f_1)
    
    for i in range(len(time)):
        if t_1<=time[i]<=t_2:
            m_1 = -delta/((t_T-t_F)/2)
            c_1 = 1-m_1*t_1
            f_2 = m_1*time[i]+c_1
            f.append(f_2)
    
    for i in range(len(time)):
        if t_2<=time[i]<=t_3:
            f.append(1-delta)
    
    for i in range(len(time)):
        if t_3<=time[i]<=t_4:
            m_2 = delta/((t_T-t_F)/2)
            c_2 = 1-m_2*t_4
            f_3 = m_2*time[i]+c_2
            f.append(f_3)   
    
    for i in range(len(time)):
        if t_4<=time[i]<=t_5:
            f.append(flux[i])
    return f

flux = light_curve(R_star=1*solar_radius,time=time,delta=0.005,flux=flux,t_T=20,i=90,P=60,m_star=1*solar_mass)
fig = plt.figure()
plt.plot(time,flux)
plt.xlabel("Time")
plt.ylabel("Flux")
plt.show()
fig.savefig("Transit_light_curve.png")
