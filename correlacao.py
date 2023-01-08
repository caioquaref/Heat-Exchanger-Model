from modelo_novo import e_nut2, temp_comp2
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  # biblioteca para construção de gráficos
from sklearn.metrics import mean_absolute_percentage_error, mean_absolute_error
mpl.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20

# dados dos testes; ordem dos arquivos: werneck -> afonso, cima -> baixo, esquerda -> direita (abas)
temp_ar_ent = [42, 42, 42, 41, 40, 40, 39, 38, 37] # input modelo  ||| 27,26,26,26,32,31,31,29,31,30,29,33,31,31,30
temp_ar_real = [43, 43, 43, 43, 42, 42, 42, 40, 40] #  ||| 43, 43, 43, 43, 42, 42, 42, 40, 40
vel = [5.7,7.2,10.9,10.1,12.6,14.5,17.9,19,11.9,18.1,17,8.6,9.4,14.1,14.9] # input modelo
temp_sai_motor = [87.2,87.2, 90.9, 90, 91.7, 93.2, 96.2, 93.9, 93.5] # input modelo  ||| 87.3,85.3,90.2,85.3,87.8,84.5,89.8,84.3,85.2,90.9,85.1,87.7,85.5,98.8,85.7
temp_rad = [85.2,85.2, 89.3, 88.4, 90.1, 91.9, 94.9, 91.2, 92.5] #  ||| 85.2,85.2, 89.3, 88.4, 90.1, 91.9, 94.9, 91.2, 92.5
temp_ent_motor = [82.8,82.8, 87.2, 86.4, 88.5, 90.1, 93.2, 89.7, 91]  #  ||| 82.8,82.8, 87.2, 86.4, 88.5, 90.1, 93.2, 89.7, 91
rot = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000] # input modelo  ||| 8500,8500,12000,12000,8500,8500,12300,12300,8500,12000,12000,8500,8500,11800,11800
vel_sim = 4.3

# utilizando os inputs dos testes para gerar outputs no modelo ;; (prof, compri, altura, rpm, temp, vel):
inputs = [rot, temp_sai_motor, temp_ar_ent, vel]
prof = 45
compri = 300
altura = 205
temp_ar_out1 = []
temp_ar_out2 = []
temp_ag_out1 = []
temp_ag_out2 = []
rad1 = []
rad2 = []
for i in range(0, 9):
    result = temp_comp2(prof, compri, altura, rot[i], temp_sai_motor[i], temp_ar_ent[i], vel_sim)
    temp_ar_out1.append(result[0])
    temp_ar_out2.append(result[1])
    temp_ag_out1.append(result[2])
    temp_ag_out2.append(result[3])
    rad1.append(result[4])
    rad2.append(result[5])
print(rad1, rad2)
# DESVIO MÉDIO E DIFERENÇA MÁXIMA
listar = np.subtract(temp_ar_real, temp_ar_out2)/temp_ar_real*100
listad = np.subtract(temp_rad, temp_ag_out1)/temp_rad*100
listag = np.subtract(temp_ent_motor, temp_ag_out2)/temp_ent_motor*100
t_ar_mean = np.mean(listar)
ar = np.max(abs(listar))
rad = np.max(abs(listad))
ag = np.max(abs(listag))
t_rad_mean = np.mean(listad)
t_ag_mean = np.mean(listag)
print('Desvio percentual médio:\nTemperatura do ar de saída: {}%'
      '\nTemperatura água pós radiador 1: {}%\nTemperatura água pós radiador 2: {}%'
      '\nErro máximo:\nAr: {}%\nÁgua 1: {}%'
      '\nÁgua 2: {}%'.format(round(t_ar_mean,2), round(t_rad_mean,2), round(t_ag_mean,2),
                             round(ar,2),round(rad,2),round(ag,2)))

# ERRO MÉDIO ABSOLUTO
a = mean_absolute_error(temp_ar_real, temp_ar_out1)
b = mean_absolute_error(temp_ar_real, temp_ar_out2)
c = mean_absolute_error(temp_rad, temp_ag_out1)
d = mean_absolute_error(temp_ent_motor, temp_ag_out2)
print("Erro médio absoluto:\nAr - Radiador 1: {}\nAr - Radiador 2: {}\nÁgua - Radiador 1: {}\nÁgua - Radiador 2: {}".format(round(a,1),round(b,1),round(c,1),round(d,1)))

# ERRO MÉDIO ABSOLUTO PERCENTUAL
a = mean_absolute_percentage_error(temp_ar_real, temp_ar_out1)*100
b = mean_absolute_percentage_error(temp_ar_real, temp_ar_out2)*100
c = mean_absolute_percentage_error(temp_rad, temp_ag_out1)*100
d = mean_absolute_percentage_error(temp_ent_motor, temp_ag_out2)*100
print("Erro médio absoluto percentual:\nAr - Radiador 1: {}%\nAr - Radiador 2: {}%\nÁgua - Radiador 1: {}%\nÁgua - Radiador 2: {}%".format(round(a,1),round(b,1),round(c,1),round(d,1)))

# GRAFICO
plt.plot(rot, temp_rad, 'bo--', label='pós rad. 1 agua teste')
plt.plot(rot, temp_ag_out1, 'bo-', label='pós rad. 1 agua modelo')
plt.plot(rot, temp_ent_motor,'g--', label='pós rad. 2 agua teste')
plt.plot(rot, temp_ag_out2, 'g-', label='pós rad. 2 agua modelo')
plt.plot(rot, temp_ar_real, 'rx--', label='saída ar teste')
#plt.plot(rot, temp_ar_out2, 'rx-', label='saída ar modelo')
plt.plot(rot, temp_ar_out1, 'rx-', label='saída ar modelo')
plt.yticks(np.arange(40, 100, 5))
#plt.xticks(np.arange(4000, 12000, 1))
plt.grid()
plt.legend(loc=0, prop={'size': 12})
#plt.title('Correlação das Temperaturas de Saída')
plt.ylabel('Temperatura de Saída [°C]')
plt.xlabel('Rotação do Motor [rpm]')
#plt.ylim(50, 85)
plt.show()