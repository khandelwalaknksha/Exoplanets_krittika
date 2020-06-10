import numpy as np 
import matplotlib.pyplot as plt 


solar_mass = 1.989*10**30
solar_radius = 6.957*10**8
G = 6.67408*10**(-11)
m_star = 1*solar_mass
jup_radius = 71492000
time = np.linspace(2000,3000,100000)
flux = np.ones(len(time))



f = []

def light_curve(R_star,time,R_planet,flux,i,P,m_star):
    delta = (R_planet/R_star)**2
    a = ((P**2*G*m_star)/(4*np.pi**2))**(1/3)
    r = a/R_star    
    rat = R_planet/R_star                # ratio of radius of planet with respect to star
    x1 = ((1-rat)**2-(r*np.cos(i))**2)**(1/2)
    x2 = ((1+rat)**2-(r*np.cos(i))**2)**(1/2)
    x12 = (P*R_star)/(np.pi*a)
    t_T = x2*x12
    t_F = x1*x12
    t_0 = time[0]                   # initial observing time
    t_1 = t_0 + 300                 # transit starting time
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
    return f,t_F,t_T,delta

flux,t_F,t_T,delta  = light_curve(R_star=1*solar_radius,time=time,R_planet=jup_radius,flux=flux,i=90,P=10,m_star=1*solar_mass)
fig = plt.figure()
plt.plot(time,flux)
plt.xlabel("Time")
plt.ylabel("Flux")
plt.show()
fig.savefig("Transit_light_curve.png")

# Limb Darkening Effect

# gamma is the angle b/w the rays emitted radially outward direction from centre of the stellar disk and 
# and else where of the disk



gamma = np.linspace(0,90,1000)       
mu = np.cos(np.deg2rad(gamma))

# linear Model
def linear(mu):
    u = 0.4
    l = 1-u*(1-mu)
    return l

# Quardetic Model

def quard(mu):
    u1,u2 = 0.2,0.2
    l = 1-u1*(1-mu)-u2*(1-mu)**2
    return l

#  Model

def three_par(mu):
    u1,u2,u3 = 0.1,0.1,0.2
    l = 1-u1*(1-mu)-u2*(1-mu**(3/2))-u3*(1-mu**2)
    return l

def non_lin(mu):
    u1,u2,u3,u4 = 0.05,0.15,0.15,0.08
    l = 1-u1*(1-mu**(1/2))-u2*(1-mu)-u3*(1-mu**(3/2))-u4*(1-mu**2)
    return l


plt.plot(gamma,linear(mu),label='Linear Model')
plt.plot(gamma,quard(mu),label='Quardatic Model')
plt.plot(gamma,three_par(mu),label='3 param')
plt.plot(gamma,non_lin(mu),label='4 param')
plt.xlim(0,100)
plt.xlabel("theta [deg]")
plt.ylabel("Relative intensity")
plt.title("Stellar limb darkening")
plt.legend()
plt.show()
plt.savefig("Stellar Limb Darkening"+'.png')

