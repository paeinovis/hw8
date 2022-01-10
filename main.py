import numpy as np
import math


K_ = .001        # / s
T_out = 40       # C
T_init = 0       # C
time = 10 * 60   # s


def NLOW(k, T_i, O, t):
    # Newton's Law of Warming(/Cooling) which I derived but didn't know if I was supposed to type the proof out here
    return (T_i - O) * np.exp(-k * t) + O


def Num(K, T, O, dt):
    # Regular old dT
    return -K * (T - O) * dt


a_sol = NLOW(K_, T_init, T_out, time)
print("Analytical solution:", a_sol)


t_steps = [1, .5, .25, .125]
results = []
c_temp = T_init
c_time = 0

# Warming
for current_step in t_steps:
    while c_time <= time:
        c_temp += Num(K_, c_temp, T_out, current_step)
        c_time += current_step
    results.append(c_temp)
    c_temp = T_init
    c_time = 0

print("Numerical solutions for each time step:", results)
print("\nConverge to analytical solution", a_sol)
print("Accuracy:", a_sol/results)






























