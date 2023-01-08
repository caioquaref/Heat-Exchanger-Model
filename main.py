# BIBLIOTECAS PARA O PROGRAMA
def e_nut(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    import numpy as np  # biblioteca para operações numéricas
    import pyromat as pm
    pm.config['unit_temperature'] = 'C'
    pm.config['unit_length'] = 'mm'
    pm.config['unit_pressure'] = 'kPa'
    pm.config['unit_energy'] = 'J'

    # from CoolProp.CoolProp import Props

    # VARIÁVEIS DE ENTRADA

    #   Geometria do Radiador
    l2 = prof # profundidade do radiador [mm] 45
    phi = 57  # angulação em relação à vertical [°]
    qf = 15  # quantidade de tubos por fluxo [-]
    n_p = 2  # quantidade de fluxo [-]
    b = 4.5  # espaçamento entre tubos [-]
    h_t = 1.3  # largura do tubo [mm]
    delta_w = 0.1  # espessura do tubo [mm]
    delta = 0.1  # espessura da aleta [mm]
    p_a = 10  # número de aletas por centímetro [-]
    le_p = 1  # espessura do agitador [mm]
    l_p = 0.5  # distância entre os agitadores [mm]
    theta = 20  # angulação do agitador [°]

    #   Propriedades do Fluxo - Ar e Água
    vel_ar = vel  # velocidade do fluxo de ar no radiador [m/s]
    t_ar = temp_ar  # temperatura de entrada do ar [°C]
    tmed_ar = (((t_ar + 40)- t_ar)/2)
    t_ag = temp_ag  # temperatura de entrada da água [°C]
    tmed_ag = (t_ag - (t_ag - 10))/2

    # CARACTERIZAÇÃO DO RADIADOR

    #   Geometria do Radiador - Lado Ar
    #dummy = altura
    l1 = altura  #(qf * n_p * (b + h_t) + b) # altura adaptada [mm]
    l3 = compri  #* np.sin(np.deg2rad(phi))  # comprimento do radiador [mm]

    p_f = (1 / p_a) * 10  # angulação das aletas [°]
    p_t = b + h_t
    n_pg = (l1 - h_t) / p_t  # número de passagens, diferente do número de passes [-]
    n_f = l3 * n_pg / p_f

    n_louv = (((l2 / le_p) - 1) * n_f)
    s_f = np.sqrt(b**2 + p_f**2)
    l_louv = 0.85 * s_f
    l_h = le_p * np.sin(np.deg2rad(theta))

    a_f = (2 * (s_f * l2 + s_f * delta) * n_f + 2 * l_louv * delta * n_louv) / 1000000
    a_p = ((2 * (l2 - h_t) + np.pi * h_t) * (l3 * (n_pg + 1)) - 2 * delta * l2 * n_f) / 1000000
    a_t2 = (a_p + a_f)*1.2
    a_c2 = ((b*l3*n_pg - (delta*(s_f - l_louv)+l_louv*l_h)*n_f)/1000000)
    a_fr = (l1 / 1000) * (l3 / 1000)

    #   Geometria do Radiador - Lado Água
    n_t = n_pg + 1
    a_t1 = ((2 * (l2 - h_t) + np.pi * (h_t - delta_w)) * l3 * n_t)/1000
    m_dot_1 = 6 * 10**(-5) * rpm - 0.1194  # fluxo de massa de água na bomba - valores para equação linearizada [kg/s]
    a_c1 = (((l2-h_t)*(h_t-2*delta_w)+(np.pi/4)*((h_t-2*delta_w)**2))*n_t/n_p)/1000000

    # dados do material
    k_al = 240

    # propriedades da agua
    H2O = pm.get('mp.H2O')
    r = H2O.d(T=tmed_ag)
    c = H2O.cp(T=tmed_ag)
    rho_ag = r[0]
    mu_ag = 0.0003145
    k_ag = 0.6613
    c_p_ag = c[0]

    # analise do escoamento
    p_r_ag = c_p_ag*mu_ag/k_ag
    g_1 = m_dot_1/a_c1
    v_1 = g_1/rho_ag
    d_h1 = (4*a_c1*l3*n_p/a_t1)
    re_ag = (g_1*d_h1/mu_ag)/1000
    q_1 = (m_dot_1/rho_ag)
    f_1 = (1.58*np.log(re_ag)-3.28)**(-2)
    nus_ag = (f_1/2)*((re_ag-1000)*p_r_ag)/(1+12.7*((f_1/2)**(0.5)*(p_r_ag**(2/3)-1)))
    h_1 = nus_ag*k_ag/(d_h1/1000)
    a_w = (2*l2*l3*n_t)/1000000
    r_w = (delta_w/1000)/(k_al*a_w)
    r_ag = 1/(h_1*a_t1)

    # propriedades do ar
    air = pm.get('ig.air')
    c = air.cp(T=tmed_ar,p=100)
    r = air.d(T=tmed_ar,p=100)
    mu_ar = 0.00001872
    c_p_ar = c[0] #1005
    k_ar = 0.02588
    rho_ar =  r[0] #1.161

    # analise do escoamento
    q_2 = a_c2*vel_ar
    m_dot_2 = q_2*rho_ar
    g_2 = m_dot_2/a_c2
    d_h2 = (4*a_c2*l2/a_t2)
    re_ar = (g_2*d_h2/mu_ar)/1000
    p_r_ar = c_p_ar*mu_ar/k_ar
    sigma_2 = a_c2/a_fr
    a_hx2 = (l2*l3*b*n_pg)/1000000000
    beta_2 = a_t2/a_hx2
    l_s = (s_f/2)/1000
    jre_ar = (re_ar**(-0.49))*((theta/90)**(0.27))*((p_f/le_p)**(0.14))*((b/le_p)**(-0.29))*((l2/le_p)**(-0.23))*((l_louv/le_p)**(0.68))*((p_t/le_p)**(-0.28))*((delta/le_p)**(-0.05))
    h_2 = jre_ar*g_2*c_p_ar/((p_r_ar**(2/3)))
    m_f = (2*h_2/(k_al*(delta/1000)))**(0.5)
    eta_f = np.tanh(m_f*(l_s))/(m_f*l_s)
    eta_o = 1-((a_f/a_t2)*(1-eta_f))
    r_ar = 1/(eta_o*h_2*a_t2)  #*0.92

    # metodo e-nut
    C_1 = m_dot_1*c_p_ag
    C_max = C_1
    UA = 1/(r_ag*r_w+r_ar)
    C_2 = m_dot_2*c_p_ar
    C_min = C_2
    C_r = C_min/C_max
    NTU = UA/C_min
    epsilon_hx = (1-np.exp(((1/C_r)*NTU**(0.22))*(np.exp((-C_r)*NTU**(0.78))-1)))
    q = epsilon_hx*C_min*(t_ag-t_ar)
    t_ag_out = (t_ag - (epsilon_hx*(C_min/C_1)*(t_ag-t_ar)))  #*1.28
    t_ar_out = (t_ar + (epsilon_hx*(C_min/C_2)*(t_ag-t_ar)))*0.45
    return q, t_ar_out, t_ag_out, epsilon_hx, m_dot_1, m_dot_2

def temp_comp(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    temp_ar_out1 = []
    temp_ar_out2 = []
    temp_ag_out1 = []
    temp_ag_out2 = []
    rad1 = e_nut(prof, compri, altura, rpm, temp_ag, temp_ar, vel)
    rad2 = e_nut(prof, compri, altura, rpm, rad1[2], temp_ar, vel)
    temp_ar_out1.append(rad1[1])
    temp_ar_out2.append(rad2[1])
    temp_ag_out1.append(rad1[2])
    temp_ag_out2.append(rad2[2])
    return temp_ar_out1, temp_ar_out2, temp_ag_out1, temp_ag_out2

def hx1(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    hx = []
    cool_flow = []
    for i in rpm:
        rad1 = e_nut(prof, compri, altura, i, temp_ag, temp_ar, vel)
        rad2 = e_nut(prof, compri, altura, i, rad1[2], temp_ar, vel)
        q_tot = (rad1[0] + rad2[0]) / 1000
        flow = rad1[4]
        hx.append(q_tot)
        cool_flow.append(flow)
    return hx, cool_flow

def hx2(prof, compri, altura, rpm, temp_ag, temp_ar, vel):
    hx = []
    air_flow = []
    cool_flow_plot = []
    for i in vel:
        rad1 = e_nut(prof, compri, altura, rpm, temp_ag, temp_ar, i)
        rad2 = e_nut(prof, compri, altura, rpm, rad1[2], temp_ar, i)
        q_tot = (rad1[0] + rad2[0]) / 1000
        flow1 = rad1[5]
        flow2 = rad1[4]
        hx.append(q_tot)
        air_flow.append(flow1)
        cool_flow_plot.append(flow2)
    return hx, air_flow, cool_flow_plot

