from modelo_novo import e_nut2, hx22
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  # biblioteca para construção de gráficos

prof = 45
compri = 300
altura = 205
t_ar = 35
vel = np.arange(3, 22, 1)
rpm = np.arange(3000, 12001, 1600)
temp = 95
result = []
x = []
y = []
for i in rpm:
    q = hx22(prof, compri, altura, i, temp, t_ar, vel)
    result.append(q[0])
    x.append(q[1])
    y.append(q[2])

mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20
plt.plot(vel, result[0], '-', label='0,06 m/s')
plt.plot(vel, result[1], '--', label='0,16 m/s')
plt.plot(vel, result[2], 'o-', label='0,25 m/s')
plt.plot(vel, result[3], 'x-', label='0,35 m/s')
plt.plot(vel, result[5], '.-', label='0,54 m/s')
plt.ylabel('Potência Dissipada [kW]')
plt.xlabel('Velocidade do ar [m/s]')
plt.legend()
plt.grid()
plt.show()