import numpy as np
import matplotlib.pyplot as plt
import pdb

data = np.loadtxt("D:\\Onedrive\\Studium\\090\\Monte Carlo\\Day3\\ksca_kabs_lambda300.dat")
#google how to import a fuction to another
N_total = 100
tau = 100 #assume tau is 10
theta_0 = (60/180)*np.pi
delta_z = 10 #default = 1 cloud depth
beta_ext = tau/delta_z #beta_ext can be calculated from tau and delta_z
s = tau/beta_ext #get the distance s = delta_z
p_0 =  np.array([0,0,delta_z]) #initial position
g=0.85 #cloud factor
dummy = np.array([0,1,0])
phi = 0 #for phi = 0
dummy = np.array([0,1,0])

def generate_random_theta():
    theta = np.random.uniform(0,1)
    return theta
theta = generate_random_theta()


def direction_vector(theta, phi):
    n = np.array([np.sin(theta)*np.sin(phi), np.sin(theta)*np.cos(phi), np.cos(theta)])
    n_dir = n
    print(f"direction vector is {n_dir}")

    return n_dir
n_inc = direction_vector(theta=theta_0, phi=0) 

#3.generate ramdom tau // from task3 day 
def generate_tau():
    r = np.random.uniform(low=0, high=1)
    tau = -np.log(1-r)
    return tau

def position_vector(p, n, s):
    #p for position, n for direction, s for distance
    p_new = p+ s*n
    return p_new
n_inc = direction_vector(theta=theta_0, phi=0)
#z_p =  position_vector(p_0, n_inc, -s) #traversed position
#z_position = float(position_vector(p_0, n_inc, -s)[2]) #only delta z

def generate_coordinate(n_inc):
    n_0 = p_0
    e1 = n_inc
    e1 = np.sqrt(1/np.dot(e1,e1))*e1 #normalize
    e2 = dummy-(np.dot(e1,dummy)/np.dot(e1,e1))*e1
    e2 = np.sqrt(1/np.dot(e2,e2))*e2
    e3 = np.cross(e1,e2)
    e3 = np.sqrt(1/np.dot(e3,e3))*e3
    #make sure e1 is still orthogonal
    c1 = np.dot(e1,e2)
    c2 = np.dot(e2,e3)
    c3 = np.dot(e1,e3)
    print(c1,c2,c3)
    return e1,e2,e3
#e1,e2,e3 = generate_coordinate(n_inc)
#pdb.set_trace() #debugger

def calculate_n_phi(phi, e2, e3): #need e2 and e3
    n_phi = np.cos(phi)*e2 + np.sin(phi)*e3 #direction vecotr n_phi
    return n_phi
#n_phi = calculate_n_phi(phi, e2, e3)

def calculate_n_scatter(n_phi, theta, n_inc):
    n_scatter = np.cos(theta)*n_inc + np.sin(theta)*n_phi
    n_scatter = np.sqrt(1/np.dot(n_scatter,n_scatter))*n_scatter #normierung
    print(n_scatter)
    return n_scatter
#n_scatter = calculate_n_scatter(n_phi, theta, n_inc)

def random_mu(g):
    r = np.random.uniform(0,1)
    mu_r=1./(2.*g)*(1.+g*g-((1.-g*g)/(1.-g+2.*r*g))**2)
    return mu_r
mu = random_mu(g)
theta = np.arccos(mu)

### Main Fuction
N = N_total
N_t = 0
N_b = 0
theta_all = np.array([])
n_inc = direction_vector(theta=theta_0, phi=0)
while N-(N_t+N_b)>0:
    p_0 = p_0
    tau = generate_tau() 
    s = tau/beta_ext #get the random distance s
    position = position_vector(p_0, n_inc, -s)
    #print(position)
    n_scatter = n_inc
    while position[2]>0: #where the scattering starts
        tau = generate_tau()
        s = tau/beta_ext #get the distance s
        #position = position_vector(p_0, n_inc, s) #get new position vector
        
        theta = np.arccos(random_mu(g)) #random theta
        phi = generate_random_theta() * np.pi*2 #random phi
        
        e1,e2,e3 = generate_coordinate(n_scatter)
        n_phi = calculate_n_phi(phi, e2, e3)
        n_phi = np.sqrt(1/np.dot(n_phi,n_phi))*n_phi #normierung
        n_in = n_scatter
        n_scatter = calculate_n_scatter(n_phi, theta, n_scatter) #Theta according to Greenstein #normierung
        check = np.dot(n_scatter,n_in)
        print(np.cos(theta),check)
        #pdb.set_trace() #debugger
        
        p_scatter =  position_vector(position, n_scatter,-s) #traversed position
        position = p_scatter
        print(f"{position[2]}")
        theta_all = np.append(theta_all, theta)
        #pdb.set_trace() #debugger
        if position[2]>=delta_z:
            N_t = N_t +1
            break
        elif position[2] <= 0:
            N_b = N_b +1   
            break       
    else: print(f"position is: {position}")
else: 
    reflectance = N_t/N
    transmittance = N_b/N
    #calculate standard derivtion
    if N_b == 0:
        print("N_b is 0")
        sigma =0
    else: 
        sigma = np.sqrt((N-N_b)/(N*N_b))
    print(f"total photon number is {N}")
    print(f"{N_t} photon reflected, Reflectance is {reflectance}")
    print(f"{N_b} photon transmitted, Transmittance is {transmittance}")
    print(f"standard deviation is {sigma}")
    #print(theta_all)
    #plt.hist(theta_all)
    #plt.show()
#pdb.set_trace() #debugger
