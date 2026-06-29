import numpy as np
from scipy import signal 
import matplotlib.pyplot as plt
#Define variables
ms=77.5
mu=10
ks=11000
ku=150000
cu=49
#Initialise Time vector
T=np.linspace(0,5,10000,endpoint=True)
x0=np.array([0,0,0,0])
#change U value
U=np.zeros((10000,2))
#change U (input) to simulate a bump
H=0.08
v=10
L=0.8
t_start=0.5
t_end=t_start+L/v
for i, t in enumerate(T):
    if t_start <= t <= t_end:
        xg=H/2*(1-np.cos((2*np.pi*v)/L*(t-t_start)))
        xg_dot=(np.pi*v*H)/L*np.sin((2*np.pi*v)/L*(t-t_start))
        U[i,0]=xg
        U[i,1]=xg_dot
#store values for later 
zeta_list = []
settling_times = []
peak_accelerations_g = []
max_spring_forces = []
max_damper_forces = []        
#add zeta value and loop for solver
zeta_range=np.linspace(0.05, 1.5, 100)
for zeta in zeta_range:
    cs=2*zeta*np.sqrt(ks*ms)
#Define matrices
    A=np.array([               
        [0,0,1,0],
        [0,0,0,1],
        [-ks/ms,ks/ms,-cs/ms,cs/ms],
        [ks/mu,-(ks+ku)/mu,cs/mu,-(cs+cu)/mu]
    ])
    B=np.array([
        [0,0],
        [0,0],
        [0,0],
        [ku/mu,cu/mu]
    ])
    C=np.array([
        [1,-1,0,0],
        [0,0,1,0]
    ])
    D=np.array([
        [0,0],
        [0,0]
    ])
    #set up system
    sys=signal.StateSpace(A,B,C,D)
    #run system
    tout,yout,xout= signal.lsim(sys, U, T, x0, interp=True)
    #extract xs,xu,xs_dot,xu_dot
    xs=xout[:,0]
    xu=xout[:,1]
    xs_dot=xout[:,2]
    xu_dot=xout[:,3]
    # calculate forces
    F_spring= ks*(xs-xu)
    F_damper=cs*(xs_dot-xu_dot)
    #calculate max forces
    F_spring_max=np.max(np.abs(F_spring))
    F_damper_max=np.max(np.abs(F_damper))
    #peak body acceleration
    xs_dotdot=(-F_spring - F_damper)/ms
    #settling time
    threshold=0.001 #from matlab code
    outside_threshold=np.where(np.abs(xs) > threshold)[0]
    if len(outside_threshold) > 0: 
       settling_time=tout[outside_threshold[-1]]
    else:
       settling_time=0
    #save values into the tracking arrays
    zeta_list.append(zeta)
    settling_times.append(settling_time)
    peak_accelerations_g.append(np.max(np.abs(xs_dotdot))/9.81) #in g's to compare to datasheet
    max_spring_forces.append(F_spring_max)
    max_damper_forces.append(F_damper_max)
#convert lists to array
settling_times=np.array(settling_times)
peak_accelerations_g=np.array(peak_accelerations_g)
max_spring_forces = np.array(max_spring_forces)  
max_damper_forces = np.array(max_damper_forces)
zeta_list=np.array(zeta_list)
#Multiple scenario responses
#Scenario 1: Check if both criteria can be met simultaneously (<0.3g AND <1.0s)
ideal_indices = np.where((peak_accelerations_g < 0.3) & (settling_times < 1.0))[0]

#Scenario 2: Settling time balanced between 0.9s and 1.0s
balanced_indices = np.where((settling_times >= 0.9) & (settling_times <= 1.0))[0]

if len(balanced_indices) > 0:
# From the configurations that settle around 1 second, pick the one with the lowest acceleration
    idx_balanced = balanced_indices[np.argmin(peak_accelerations_g[balanced_indices])]
else:
# Fallback: just find the one closest to 1.0s overall
    idx_balanced = np.argmin(np.abs(settling_times - 1.0))

#Scenario 3: Minimum Settling Time
idx_min_settling = np.argmin(settling_times)

#Scenario 4: Minimum Acceleration
idx_min_accel = np.argmin(peak_accelerations_g)    

#Print the Results
#Output 1 - Check to see if best case possible
print("OUTPUT 1: IDEAL SPECIFICATION CHECK (Comfort < 0.3g AND Settling < 1s)")
if len(ideal_indices) > 0:
    print(f"Success! Ideal configuration found at Zeta: {zeta_list[ideal_indices[0]]:.3f}")
else:
    print("RESULT: SIMULTANEOUS TARGETS MATHEMATICALLY IMPOSSIBLE.")
    print("The system cannot physically achieve <0.3g and settle in <1.0s simultaneously.")

#Output 2: Balanced case (Used for optimisation)
print("OUTPUT 2: Balanced (Settling Time ~1.0s, Optimized Accel)")
print(f"Optimal Zeta:         {zeta_list[idx_balanced]:.3f}")
print(f"Settling Time:        {settling_times[idx_balanced]:.3f} s  <-- [Target Met: 0.9s - 1.0s]")
print(f"Peak Acceleration:    {peak_accelerations_g[idx_balanced]:.3f} g <-- [Minimized Compromise]")
print(f"Max Spring Force:     {max_spring_forces[idx_balanced]:.3f} N")
print(f"Max Damper Force:     {max_damper_forces[idx_balanced]:.3f} N")

#Output 3: Minimum Settling time
print("OUTPUT 3: Minimum Settling Time")
print(f"Optimal Zeta:         {zeta_list[idx_min_settling]:.3f}")
print(f"Settling Time:        {settling_times[idx_min_settling]:.3f} s")
print(f"Peak Acceleration:    {peak_accelerations_g[idx_min_settling]:.3f} g")
print(f"Max Spring Force:     {max_spring_forces[idx_min_settling]:.3f} N  ")
print(f"Max Damper Force:     {max_damper_forces[idx_min_settling]:.3f} N  ")

#Output 4: Minimum peak body acceleration
print("OUTPUT 4: Minimum Peak Body Acceleration")
print(f"Optimal Zeta:         {zeta_list[idx_min_accel]:.3f}")
print(f"Settling Time:        {settling_times[idx_min_accel]:.3f} s")
print(f"Peak Acceleration:    {peak_accelerations_g[idx_min_accel]:.3f} g")
print(f"Max Spring Force:     {max_spring_forces[idx_min_accel]:.3f} N  ")
print(f"Max Damper Force:     {max_damper_forces[idx_min_accel]:.3f} N  ")

#Plotting x values, chassis displacement (xs) and wheel displacment (xu)
cs_balanced = 2 * zeta_list[idx_balanced] * np.sqrt(ks * ms)
A_balanced = np.array([               
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [-ks/ms, ks/ms, -cs_balanced/ms, cs_balanced/ms],
    [ks/mu, -(ks+ku)/mu, cs_balanced/mu, -(cs_balanced+cu)/mu]
])
sys_balanced = signal.StateSpace(A_balanced, B, C, D)
tout, yout, xout = signal.lsim(sys_balanced, U, T, x0, interp=True)
# Plotting x values, chassis displacement (xs) and wheel displacment (xu)
chassis_displacement=xout[:,0]
wheel_displacement=xout[:,1]
suspension_travel=yout[:,0]
chassis_velocity=yout[:,1]
fig, ax = plt.subplots(2,1)
#Graph 1 for x values
ax[0].plot(tout, chassis_displacement, label="Chassis Displacement ($x_s$)", color="blue", linewidth=2)
ax[0].plot(tout, wheel_displacement, label="Wheel Displacement ($x_u$)", color="orange", linestyle="--")
ax[0].set_xlabel("Time (Seconds)")
ax[0].set_ylabel("Displacement (Metres)")
ax[0].legend()
ax[0].set_title("Balanced Dynamic Response")
ax[0].grid(True)
#graph 2.1
line1=ax[1].plot(tout, suspension_travel, label="Suspension Travel ($x_s$-$x_u$)", color="purple", linewidth=2)
ax[1].set_xlabel("Time (Seconds)")
ax[1].set_ylabel("Displacement (Metres)")
#graph 2.2
ax_twin=ax[1].twinx()
line2= ax_twin.plot(tout, chassis_velocity, label=r"Chassis Velocity ($\dot{x}_s$)", color="green", linestyle="--") #if look weird add r"chassis
ax_twin.set_ylabel("Chassis Velocity (Metres per Second)")
#combine lines
lines=line1+line2
labels=[l.get_label() for l in lines]
ax[1].legend(lines, labels)
ax[1].grid(True)

#plots for zeta against settling time and acceleration
fig, ax1 = plt.subplots()

ax1.plot(zeta_list, peak_accelerations_g, color="red")
ax1.set_xlabel("Damping ratio")
ax1.set_ylabel("Peak acceleration (g)")
ax1.set_title("Peak Body Acceleration and Settling Time vs. Damping Ratio")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(zeta_list, settling_times, color="blue")
ax2.set_ylabel("Settling time (s)")

plt.show()

