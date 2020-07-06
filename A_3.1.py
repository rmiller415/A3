# -*- coding: utf-8 -*-
"""
Created on Sun May 31 09:57:25 2020

@author: Richard Miller
"""

######################IMPORTS###############################
import numpy as np
import scipy.optimize as sciop
import os
import winsound as ws
import random
import time
start = time.time()

######################DEFINE VARIABLES######################

w_3 = 0.1
w_2 = -0.4
l_b = 7.2
l_0 = 3.8
i = 50
j = 1


#####################INPUT SEQUENCE########################
raw_input = "MSETITVNCPTCGKTVVWGEISPFRPFCSKRCQLIDLGEWAAEEKRIPSSGDLSESDDWSEEPKQ"
N = len(raw_input)


#####################PARSE INPUT SEQUENCE##################
def Sequence_Parser(raw_input):
    sequence_raw = []
    sequence_numeric = []
    
    
    p = 0
    while p in range(N):
        sequence_raw.append(str(raw_input[p]))
        p += 1
        
        
        
    p = 0   
    while p in range(N):
        
        
        if sequence_raw[p] == 'E':
            sequence_numeric.append(-1)
            
        elif sequence_raw[p] == 'D':
            sequence_numeric.append(-1)
            
        elif sequence_raw[p] == 'R':
            sequence_numeric.append(1)

        elif sequence_raw[p] == 'K':
            sequence_numeric.append(1)

        else:
            sequence_numeric.append(0)
            
            
        p += 1
    return sequence_numeric


#####################DEFINE OMEGA FUNCTIONS##############
def Omega(i,j):
    m_1 = j
    omega_1 = 0
    
    m_2 = (j + 1)
    omega_2 = 0
    
    m_3 = (i + 1)
    omega_3_constant = (i - j)**2
    omega_3 = 0
    
    m_4 = (i + 1)
    omega_4 = 0
    
    while m_1 <= i:
        n_1 = 1
        while n_1 <= (j - 1):
            omega_1_temp = ((m_1 - j)**2)/((m_1 - n_1)**(5/2))
            omega_1 = omega_1 + omega_1_temp
            n_1 += 1
        m_1 += 1
        
      
    while m_2 <= i:
        n_2 = j
        while n_2 <= (m_2 - 1):
            omega_2_temp = (m_2 - n_2)**(-1/2)
            omega_2 = omega_2 + omega_2_temp
            n_2 += 1
        m_2 += 1
        
    
    while m_3 <= N:
        n_3 = 1
        while n_3 <= (j - 1):
            omega_3_temp = omega_3_constant/((m_3 - n_3)**(5/2))
            omega_3 = omega_3 + omega_3_temp
            n_3 += 1
        m_3 += 1
        
    while m_4 <= N:
        n_4 = j
        while n_4 <= i:
            omega_4_temp = ((i - n_4)**2)/((m_4 - n_4)**(5/2))
            omega_4 = omega_4 + omega_4_temp
            n_4 += 1
        m_4 += 1
    result = w_2*(omega_1 + omega_2 + omega_3 + omega_4)
    
    return result



#####################DEFINE T FUNCTIONS#####################
def T_ij(i,j):
    l_1 = (j + 1)
    T_1 = 0
    
    l_2 = (j + 2)
    T_2 = 0
    
    l_3 = (i + 1)
    T_3 = 0
    
    l_4 = (i + 1)
    T_4 = 0
    
    l_5 = (i + 1)
    T_5 = 0
    
    l_6 = (j + 3)
    T_6 = 0
    
    l_7 = (i + 1)
    T_7 = 0
    
    l_8 = (i + 2)
    T_8 = 0
    
    T_const = (i - j)**2
    
    while l_1 <= i:
        m_1 = 2
        while m_1 <= j:
            n_1 = 1
            while n_1 <= (m_1 - 1):
                T_1_temp = (l_1 - j)**2/(((l_1 - m_1)**(5/2))*(m_1 - n_1)**(3/2))
                T_1 = T_1 + T_1_temp
                n_1 += 1
            m_1 += 1
        l_1 += 1
    
    while l_2 <= i:
        m_2 = (j + 1)
        while m_2 <= (l_2 - 1):
            n_2 = 1
            while n_2 <= j:
                T_2_temp = ((l_2 - m_2)**2)/(((l_2 - m_2)**(5/2))*(m_2 - n_2)**(3/2)) + ((m_2 - j)**2)/(((l_2 - m_2)**(3/2))*(m_2 - n_2)**(5/2))
                T_2 = T_2 + T_2_temp
                n_2 += 1
            m_2 += 1
        l_2 += 1
    
    while l_3 <= N:
        m_3 = 2
        while m_3 <= j:
            n_3 = 1
            while n_3 <= (m_3 - 1):
                T_3_temp = T_const/(((l_3 - m_3)**(5/2))*(m_3 - n_3)**(3/2))
                T_3 = T_3 + T_3_temp
                n_3 += 1
            m_3 += 1
        l_3 += 1
    
    while l_4 <= N:
        m_4 = i
        while m_4 <= (l_4 - 1):
            n_4 = 1
            while n_4 <= j:
                T_4_temp = T_const/(((l_4 - m_4)**(3/2))*(m_4 - n_4)**(5/2))
                T_4 = T_4 + T_4_temp
                n_4 += 1
            m_4 += 1
        l_4 += 1
        
    while l_5 <= N:
        m_5 = (j + 1)
        while m_5 <= (i - 1):
            n_5 = 1
            while n_5 <= j:
                T_5_temp = ((i - m_5)**2)/(((l_5 - m_5)**(5/2))*(m_5 - n_5)**(3/2)) + ((m_5 - j)**2)/(((l_5 - m_5)**(3/2))*(m_5 - n_5)**(5/2))
                T_5 = T_5 + T_5_temp
                n_5 += 1
            m_5 += 1
        l_5 += 1
    
    while l_6 <= i:
        m_6 = (j + 2)
        while m_6 <= (l_6 - 1):
            n_6 = (j + 1)
            while n_6 <= (m_6 - 1):
                T_6_temp = ((l_6 - m_6)**2)/(((l_6 - m_6)**(5/2))*(m_6 - n_6)**(3/2)) + ((m_6 - n_6)**2)/(((l_6 - m_6)**(3/2))*(m_6 - n_6)**(5/2))
                T_6 = T_6 + T_6_temp
                n_6 += 1
            m_6 += 1
        l_6 += 1
        
    while l_7 <= N:
        m_7 = (j + 2)
        while m_7 <= i:
            n_7 = (j + 1)
            while n_7 <= (m_7 - 1):
                T_7_temp = ((i - m_7)**2)/(((l_7 - m_7)**(5/2))*(m_7 - n_7)**(3/2)) + ((m_7 - n_7)**2)/(((l_7 - m_7)**(3/2))*(m_7 - n_7)**(5/2))
                T_7 = T_7 + T_7_temp
                n_7 += 1
            m_7 += 1
        l_7 += 1
    
    while l_8 <= N:
        m_8 = (i + 1)
        while m_8 <= (l_8 - 1):
            n_8 = (j + 1)
            while n_8 <= i:
                T_8_temp = ((i - n_8)**2)/(((l_8 - m_8)**(3/2))*(m_8 - n_8)**(5/2))
                T_8 = T_8 + T_8_temp
                n_8 += 1
            m_8 += 1
        l_8 += 1
        
    result = T_1 + T_2 + T_3 + T_4 + T_5 + T_6 + T_7 + T_8
    return result


#####################DEFINE Q FUNCTIONS#####################
def Q_ij(i,j,raw_input):
    q = Sequence_Parser(raw_input)
    A_coef = 0.5*((6*np.pi)**(1/2))
    Q_coef = (i - j)**2
    
    m_1 = j
    Q_1 = 0
    
    m_2 = (j + 1)
    Q_2 = 0
    
    m_3 = (i + 1)
    Q_3 = 0
    
    m_4 = (i + 1)
    Q_4 = 0
    
    while m_1 <= i:
        q_m = q[(m_1 - 1)]
        n_1 = 1
        while n_1 <= (j - 1):
            q_n = q[(n_1 - 1)]
            Q_1_temp = q_m*q_n*A_coef*((m_1 - j)**2)/((m_1 - n_1)**(3/2))
            Q_1 = Q_1 + Q_1_temp
            n_1 += 1
        m_1 += 1
    
    while m_2 <= i:
        n_2 = j
        q_m = q[(m_2 - 1)]
        while n_2 <= (m_2 - 1):
            q_n = q[(n_2 - 1)]
            Q_2_temp = q_m*q_n*A_coef*((m_2 - n_2)**2)/((m_2 - n_2)**(3/2))
            Q_2 = Q_2 + Q_2_temp
            n_2 += 1
        m_2 += 1
        
    while m_3 <= N:
        n_3 = 1
        q_m = q[(m_3 - 1)]
        while n_3 <= (j - 1):
            q_n = q[(n_3 - 1)]
            Q_3_temp = q_m*q_n*A_coef*Q_coef/((m_3 - n_3)**(3/2))
            Q_3 = Q_3 + Q_3_temp
            n_3 += 1
        m_3 += 1
    
    while m_4 <= N:
        n_4 = j
        q_m = q[(m_4 - 1)]
        while n_4 <= i:
            q_n = q[(n_4 - 1)]
            Q_4_temp = q_m*q_n*A_coef*((i - n_4)**2)/((m_4 - n_4)**(3/2))
            Q_4  = Q_4 + Q_4_temp
            n_4 += 1
        m_4 += 1
    
    result = Q_1 + Q_2 + Q_3 + Q_4
    return result


#####################INITIAL GUESS FUNCTION################
def User_Input():
    while_cond = True
    
    
    while while_cond == True:
        user_input = input("Input an initial guess:")
        
        
        try:
            float(user_input)
            while_cond = False
        except ValueError:
            print("This is not a valid input, try again.")
            
    
    return user_input


#####################RUN MINIMIZER##########################
def Function(x_ij,i,j,raw_input):
    Q_ij_ = Q_ij(i,j,raw_input)
    T_ij_ = T_ij(i,j)
    O_ij_ = Omega(i,j)
    FiTC = 3/2
    STC = O_ij_*((3/(2*np.pi))**(3/2))/(i - j)
    TTC = T_ij_*((3/(2*np.pi))**3)*(w_3/(2*(i - j)))   
    FoTC = Q_ij_*(2*l_b)/(np.pi*l_0*(i - j))
    function = FiTC*(x_ij - np.log(x_ij)) + STC*x_ij**(-3/2) + TTC*x_ij**(-3) + FoTC*x_ij**(-1/2)
    return function


def Minimize(i,j,raw_input,x_ij):
    result = sciop.minimize(Function, x_ij, args = (i,j,raw_input), method='SLSQP')
    return result


#########################ITERATE OVER i,j###################
def Crawler(i,j,raw_input):
    arr = np.zeros([65,65], dtype=float)
    x_ij = User_Input()
    i = 2
    while i <= N:
        j = 1
        while j <= i:
            if j == (i - 1):
                arr[(i-1),(j-1)] = 1
            elif j == i:
                arr[(i-1),(j-1)] = 0
            else:
                temp = Minimize(i,j,raw_input, x_ij)
                arr[(i-1),(j-1)] = temp.x[0]
            j += 1
        i+= 1
    np.save(os.path.join('C:\Folder', 'a3results1.npy'), arr)
    return arr

RESULT = Crawler(i,j,raw_input)
finish = time.time()
k = 0
first_pass = True
while k <= 10:
    freq = 414
    dur = 500
    freq_pn_test_result = 0
    freq_pn_test = random.randint(1,51)
    dur_pn_test_result = 0
    dur_pn_test = random.randint(1,51)
    
    if first_pass == False:
        if freq_pn_test%2 == 0:
            freq_pn_test_result = (-1)
        else:
            freq_pn_test_result = 1 
        freq = freq + freq_pn_test_result*random.randint(0,400)
        
        
        if dur_pn_test%2 == 0:
            dur_pn_test_result = (-1)
        else:
            freq_pn_test_result = 1
            
        dur = dur + freq_pn_test_result*random.randint(0,250)
        
        ws.Beep(freq,dur)
    else:
        ws.Beep(freq,dur)
    k += 1
    first_pass = False

print("------------------%s seconds----------------" %(finish-start))