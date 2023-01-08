def e_nut2(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    import numpy as np  # biblioteca para operações numéricas

    # VARIÁVEIS DE ENTRADA

    #  Geometria do Radiador
    l2 = prof  # profundidade do radiador [mm] 45
    l1 = compri  # (qf * n_p * (b + h_t) + b) <-- altura adaptada [mm]
    l3 = altura  # comprimento do radiador [mm]
    phi = 57  # angulação em relação à vertical [°]
    qf = 15  # quantidade de tubos por fluxo [-]
    n_p = 2  # quantidade de fluxo [-]
    b = 4.5  # espaçamento entre tubos [-]
    h_t = 1.3  # largura do tubo [mm]
    delta_w = 0.1  # espessura do tubo [mm]
    delta = 0.1  # espessura da aleta [mm]
    p_a = 10  # número de aletas por centímetro [-]
    l_p = 1  # espessura do agitador [mm]
    L_p = 0.5 # distância entre os agitadores [mm] 0.5
    theta = 20  # angulação do agitador [°]

    #  Temperaturas
    t_ar = temp_ar  # temperatura de entrada do ar [°C]
    tmed_ar = (((t_ar + 10)- t_ar)/2)
    t_ag = temp_ag  # temperatura de entrada da água [°C]
    tmed_ag = (t_ag - (t_ag - 10))/2

    # Propriedades do Fluxo - Ar e Água

    #  Propriedades da Água
    rho_ag = 965.3  # densidade [SI]
    mu_ag = 0.0003  # viscosidade [SI]
    k_ag = 0.661  # condutividade térmica [SI]
    c_p_ag = 4204  # calor específico [SI]

    #  Propriedades do Ar
    vel_ar = vel  # velocidade do fluxo de ar no radiador [m/s]
    mu_ar = 18.72E-6  # viscosidade
    c_p_ar = 1005  # calor específico
    k_ar = 0.026  # condutividade térmica
    rho_ar = 1.161  # densidade kg/m³

    #  Dados do Material do Radiador
    k_al = 240  # condutividade térmica, Alumínio comum

    # CARACTERIZAÇÃO GEOMÉTRICA

    #  Geometria do Radiador - Lado Ar
    p_f = (1 / p_a) * 10 # angulação das aletas [°]
    p_t = b + h_t #
    n_pg = (l1 - h_t) / p_t  # número de passagens, diferente do número de passes [-]
    n_f = l3*n_pg/p_f  # número total de aletas [-]
    n_louv = (((l2 / l_p) - 1) * n_f)  # número total de agitadores [-]
    s_f = np.sqrt(b**2 + p_f**2)  # largura do agitador
    l_louv = 0.85 * s_f  # comprimento do agitador
    l_h = l_p * np.sin(np.deg2rad(theta))  # altura do agitador
    a_f = (2 * (s_f * l2 + s_f * delta) * n_f + 2 * l_louv * delta * n_louv) / 1000000  # área total das aletas
    a_p = ((2 * (l2 - h_t) + np.pi * h_t) * (l3 * (n_pg + 1)) - 2 * delta * l2 * n_f) / 1000000  # área primária
    a_t2 = (a_p + a_f)  # área total de transferência de calor - lado ar
    a_c2 = ((b*l3*n_pg - (delta*(s_f - l_louv)+l_louv*l_h)*n_f)/1000000)  # área mínima de fluxo livre - lado ar
    a_fr = (l1 / 1000) * (l3 / 1000)  # área frontal do radiador

    #  Geometria do Radiador - Lado Água
    n_t = n_pg + 1  # número total de tubos [-]
    a_t1 = ((2 * (l2 - h_t) + np.pi * (h_t - 2*delta_w)) * l3 * n_t)/1000000  # área total de transferência de calor - lado ar
    m_dot_1 = (6 * 10**(-5) * rpm - 0.1194)  # fluxo de massa de água na bomba - valores para equação linearizada [kg/s]
    a_c1 = (((l2-h_t)*(h_t-2*delta_w)+(np.pi/4)*((h_t-2*delta_w)**2))*n_t/n_p)/1000000  # área mínima de fluxo livre - lado ar

    # CARACTERICAÇÃO DO ESCOAMENTO

    #  Ánalise do Escoamento - Lado Água
    p_r_ag = c_p_ag*mu_ag/k_ag  # Prandt
    g_1 = m_dot_1/a_c1  # velocidade de massa
    v_1 = g_1/rho_ag  # velocidade do escoamento
    d_h1 = (4*a_c1*l3/a_t1)  # diâmetro hidráulico
    re_ag = (g_1*d_h1/mu_ag)/1000  # Reynolds
    q_1 = 1.65E-3  # vazão volumétrica

    if re_ag >= 2100:
        f_1 = 1.58*np.log((re_ag)-3.28)**(-2)  # fator de fricção
        nus_ag = (f_1/2)*((re_ag-1000)*p_r_ag)/(1+12.7*((f_1/2)**(0.5)*(p_r_ag**(2/3)-1)))  # Nusselt
    else:
        f_1 = 1.6 / re_ag  # fator de fricção
        nus_ag = 7.541  # Nusselt

    h_1 = nus_ag*k_ag/(d_h1/1000)  # coeficiente de convectivo
    a_w = (2*l2*l3*n_t)/1000000  # área total de condução do tubo
    r_w = (delta_w/1000)/(k_al*a_w)  # resistência térmica da parede do tubo
    r_ag = 1/(h_1*a_t1)  # resistência térmica do lado água

    #  Ánalise do Escoamento - Lado Ar
    q_2 = a_c2*vel_ar  # vazão volumétrica
    m_dot_2 = q_2*rho_ar  # fluxo de massa
    g_2 = m_dot_2/a_c2  # velocidade de massa
    d_h2 = (4*a_c2*l2/a_t2)  # diâmetro hidráulico
    re_ar = (g_2*d_h2/mu_ar)/1000  # Reynolds
    p_r_ar = c_p_ar*mu_ar/k_ar  # Prandt
    sigma_2 = a_c2/a_fr
    a_hx2 = (l2*l3*b*n_pg)/1000000000
    beta_2 = a_t2/a_hx2
    l_s = (s_f/2)/1000
    p_t2 = 16.255
    jre_ar = 0.1*((re_ar**(-0.49))*((theta/90)**(0.27))*((p_f/l_p)**(-0.14))*((b/l_p)**(-0.29))*((l2/l_p)**(-0.23))*((l_louv/l_p)**(0.68))*((p_t2/l_p)**(-0.28))*((delta/l_p)**(-0.05)))  # fator j de Colburn
    h_2 = jre_ar*g_2*c_p_ar/((p_r_ar**(2/3)))  # coeficiente convectivo
    m_f = (2*h_2/(k_al*(delta/1000)))**(0.5)  # constante para cálculo da eficiência de aleta
    eta_f = np.tanh(m_f*(l_s))/(m_f*l_s)  # eficiência de uma aleta
    eta_o = 1-((a_f/a_t2)*(1-eta_f))  # eficência total das aletas
    r_ar = 1/(eta_o*h_2*a_t2)  # resistência térmica do lado ar

    #  Método e-NUT
    C_1 = m_dot_1*c_p_ag  # taxa de capacidade calorífica
    C_2 = m_dot_2*c_p_ar  # taxa de capacidade calorífica

    if C_1 > C_2:
        C_max = C_1  # taxa de capacidade calorífica máxima
        C_min = C_2  # taxa de capacidade calorífica mínima
    else:
        C_max = C_2
        C_min = C_1

    UA = 1/(r_ag+r_w+r_ar)  # coeficiente global de transferência de calor
    C_r = C_min/C_max  # taxa de capacidade calorífica crítica
    NTU = UA/C_min  # número de unidades de transferência
    epsilon_hx = (1-np.exp(((1/C_r)*NTU**(0.22))*(np.exp((-C_r)*NTU**(0.78))-1)))  # efetividade
    q = epsilon_hx*C_min*(t_ag-t_ar)  # troca de calor
    t_ag_out = (t_ag - (epsilon_hx*(C_min/C_1)*(t_ag-t_ar)))  # temperatura de saída da água
    t_ar_out = (t_ar + (epsilon_hx*(C_min/C_2)*(t_ag-t_ar)))*0.85  # temperatura de saída do ar
    return t_ar_out, t_ag_out, epsilon_hx, m_dot_1, m_dot_2, q

def temp_comp2(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    temp_ar_out1 = []
    temp_ar_out2 = []
    temp_ag_out1 = []
    temp_ag_out2 = []
    rad1 = e_nut2(prof, compri, altura, rpm, temp_ag, temp_ar, vel)
    rad2 = e_nut2(prof, compri, altura, rpm, rad1[1], temp_ar, vel)
    temp_ar_out1.append(rad1[0])
    temp_ar_out2.append(rad2[0])
    temp_ag_out1.append(rad1[1])
    temp_ag_out2.append(rad2[1])
    return temp_ar_out1, temp_ar_out2, temp_ag_out1, temp_ag_out2, rad1[5], rad2[5]

def hx12(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    hx = []
    cool_flow = []
    for i in rpm:
        rad1 = e_nut2(prof, compri, altura, i, temp_ag, temp_ar, vel)
        rad2 = e_nut2(prof, compri, altura, i, rad1[1], temp_ar, vel)
        q_tot = (rad1[5] + rad2[5])/1000
        flow = rad1[3]
        hx.append(q_tot)
        cool_flow.append(flow)
    return hx, cool_flow

def hx22(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    hx = []
    air_flow = []
    cool_flow_plot = []
    for i in vel:
        rad1 = e_nut2(prof, compri, altura, rpm, temp_ag, temp_ar, i)
        rad2 = e_nut2(prof, compri, altura, rpm, rad1[1], temp_ar, i)
        q_tot = (rad1[5] + rad2[5])/1000
        flow1 = rad1[4]
        flow2 = rad1[3]
        hx.append(q_tot)
        air_flow.append(flow1)
        cool_flow_plot.append(flow2)
    return hx, air_flow, cool_flow_plot
