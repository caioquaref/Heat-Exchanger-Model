from modelo_novo import e_nut2
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  # biblioteca para construção de gráficos

prof = 45
compri = 300
altura = 205
prof_g = 45
t_ar = 35
compri_g = 270
altura_g = 200
vel = 7
rpm = [3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000]
temp_init = 95
q = []
e = []
for i in rpm:
    rad1 = e_nut2(prof, compri, altura, i, temp_init, t_ar, vel)
    rad2 = e_nut2(prof, compri, altura, i, rad1[1], t_ar, vel)
    q_tot = (rad1[5] + rad2[5])/1000
    e_tot = (rad1[2] + rad2[2])/2
    q.append(q_tot)
    e.append(e_tot)
q_gerin = []
e_g = []
for i in rpm:
    rad3 = e_nut2(prof_g, compri_g, altura_g, i, temp_init, t_ar, vel)
    rad4 = e_nut2(prof_g, compri_g, altura_g, i, rad3[1], t_ar, vel)
    q_tot = (rad3[5] + rad4[5])/1000
    e_tot = (rad3[2] + rad4[2])/2
    q_gerin.append(q_tot)
    e_g.append(e_tot)
q_17 = [3.67, 5.07, 6.34, 8.57, 10.36, 11.93, 12.39, 12.99, 13.96, 14.92]
q_30 = [6.16, 8.96, 11.2, 15.12, 18.3, 21.05, 21.86, 22.92, 24.64, 26.33]
q_mean = np.average(q)
q_gm = np.average(q_gerin)
q_30m = np.average(q_30)
q_17m = np.average(q_17)
print(q_30[-1], q[-1], q_gerin[-1], q_17[-1])
pct1 = np.interp(q_mean, [q_17m, q_30m], [17, 35])
pct2 = np.interp(q_gm, [q_17m, q_30m], [17, 35])
print('30', int(pct1), int(pct2), '17')
std = np.mean(e)
print(std)
print(rad2[2], (rad2[1]+rad1[1])/2)

mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20
plt.plot(rpm, q, '-', label='34% - Radiador real')
plt.plot(rpm, q_gerin, 'o-',label='30% - Radiador ideal')
plt.plot(rpm, q_17, '--', label='17% $P_{rad}$')
plt.plot(rpm, q_30, 'x-', label='30% $P_{rad}$')
plt.title('')
plt.ylim(0, 40)
plt.ylabel('Potência Dissipada [kW]')
plt.xlabel('Rotação [rpm]')
plt.legend()
plt.grid()
plt.show()

plt.plot(rpm, e, '-', label='34% - Dimensões reais')
plt.plot(rpm, e_g, 'o-',label='30% - Dimensões ideais')
plt.title('')
plt.ylabel('Efetividade')
plt.xlabel('Rotação [rpm]')
plt.legend()
plt.grid()
plt.show()