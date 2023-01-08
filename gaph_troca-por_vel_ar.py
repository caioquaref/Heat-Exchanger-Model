from modelo_novo import e_nut2, hx12
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  # biblioteca para construção de gráficos

prof = 45
compri = 300
altura = 205
t_ar = 35
vel = np.arange(3, 27, 5)
rpm = np.arange(3000, 12000, 500)
temp = 95
result = []
x = []
for i in vel:
    q = hx12(prof, compri, altura, rpm, temp, t_ar, i)
    result.append(q[0])
    x.append(q[1])

mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20
plt.plot(x[0], result[0], '-', label='3 m/s')
plt.plot(x[1], result[1], '--', label='8 m/s')
plt.plot(x[2], result[2], 'o-', label='13 m/s')
plt.plot(x[3], result[3], 'x-', label='18 m/s')
plt.plot(x[4], result[4], '.-', label='23 m/s')
plt.ylabel('Potência Dissipada [kW]')
plt.xlabel('Vazão da bomba de água [kg/s]')
plt.legend()
plt.grid()
plt.show()