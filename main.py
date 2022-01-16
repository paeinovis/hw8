import numpy as np
import matplotlib.pyplot as plt
import math


K_ = .001        # / s
T_out = 40       # C
T_init = 0       # C
time = 10 * 60   # s


def NLOW(k, T_i, O, t):
    # Newton's Law of Warming(/Cooling) which I derived but didn't know if I was supposed to type the proof out here
    return (T_i - O) * np.exp(-k * t) + O


def Num(K, T, O):
    # Regular old dT/dt
    return -K * (T - O)


a_sol = NLOW(K_, T_init, T_out, time)
print("Analytical solution:", a_sol)


t_steps = [1, .5, .25, .125, .05]
results = []
c_temp = T_init
c_time = 0

# Warming
for current_step in t_steps:
    while c_time <= time:
        c_temp += Num(K_, c_temp, T_out) * current_step
        c_time += current_step
    results.append(c_temp)
    c_temp = T_init
    c_time = 0

print("\nNumerical solutions for each time step (Euler):", results)
print("Converge to analytical solution", a_sol)
print("Accuracy:", a_sol/results)


# Runge-Kutta

results_RK = []
c_temp = T_init
c_time = 0

# RK second order
for current_step in t_steps:
    while c_time <= time:
        k1 = current_step * Num(K_, c_temp, T_out)
        k2 = current_step * Num(K_, c_temp + .5 * k1, T_out)
        c_temp += k2 + (current_step ** 3)
        c_time += current_step
    results_RK.append(c_temp)
    c_temp = T_init
    c_time = 0

print("\nNumerical solutions for each time step (Runge-Kutta, 2nd order):", results_RK)
print("Converge to analytical solution", a_sol)
print("Accuracy:", a_sol/results_RK)



# Q2
def dxdt(a, x, B, y):
    return (a * x) - (B * x * y)
def dydt(g, x, d, y):
    return (g * x * y) - (d * y)

a = 1
B = .5
g = .5
d = 2

x = 2
y = 2
c_time = 0
time = 30
current_step = .25
pop_r = [2]
pop_f = [2]

# RK fourth order
while c_time <= time:
    k1_r = current_step * dxdt(a, x, B, y)
    k1_f = current_step * dydt(g, x, d, y)

    k2_r = current_step * dxdt(a, x + .5 * k1_r, B, y + .5 * k1_f)
    k2_f = current_step * dydt(g, x + .5 * k1_r, d, y + .5 * k1_f)

    k3_r = current_step * dxdt(a, x + .5 * k2_r, B, y + .5 * k2_f)
    k3_f = current_step * dydt(g, x + .5 * k2_r, d, y + .5 * k2_f)

    k4_r = current_step * dxdt(a, x + k3_r, B, y + k3_f)
    k4_f = current_step * dydt(g, x + k3_r, d, y + k3_f)

    x += (1 / 6) * (k1_r + (2 * k2_r) + (2 * k3_r) + k4_r)
    pop_r.append(x)

    y += (1 / 6) * (k1_f + (2 * k2_f) + (2 * k3_f) + k4_f)
    pop_f.append(y)

    c_time += current_step

time_stepped = np.arange(0, c_time + current_step, current_step)
plt.plot(time_stepped, pop_r, label = "Rabbits", color="brown")
plt.plot(time_stepped, pop_f, label = "Foxes", color="orange")
plt.legend()
plt.show()

print("\nReasonable graph considering the lag on foxes increasing after the rabbits increase and decreasing when the "
      "food supply decreases.")
