import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Parâmetros
mi = 2  # kg/metro
l1 = 0.102040  # metros
g = 9.8  # m/s^2
wp = np.sqrt(g / l1)
lamb = 1 + 1 / 20 + 7 / 20 + 9 / 20 + 8 / 20
l2 = l1/lamb
massa1 = mi * l1
massa2 = mi * l2

# Parâmetros de simulação
tf = 5
t = np.linspace(0, tf, 10000)

# Condições iniciais 1
teta1_01 = 0.2  # rad
teta2_01 = 0.1  # rad

# Condições iniciais do sistema 1
x01 = [teta1_01, teta2_01, 0, 0]

# Condições iniciais 2
teta1_02 = 0.5  # rad
teta2_02 = 0.3  # rad

# Condições iniciais do sistema 2
x02 = [teta1_02, teta2_02, 0.1, 0.1]

def sistema_naoLinear(t, x):
    teta1, teta2, omega1, omega2 = x

    # Equações do sistema
    sen_diff = np.sin(teta1 - teta2)
    cos_diff = np.cos(teta1 - teta2)
    sen_diff2 = np.sin(2*(teta1-teta2))
    cos_diff2 = np.cos(2*(teta1-teta2))
    omega1_sq = omega1**2

    eq1 = omega1
    eq2 = omega2
    eq3 = (12*lamb/(4*lamb+12+9*cos_diff**2))*((-3/(8*lamb)*sen_diff2*(omega1-omega2))-(omega2/(2*lamb**2)*sen_diff*(omega1-omega2))-(3/(4*lamb)*(omega1**2)*omega2*sen_diff)-(omega1*(omega2**2)/(2*lamb**2)*sen_diff)+((wp**2)*(3/(2*lamb)*np.sin(teta2)-(1/lamb+1/2)*np.sin(teta1))))
    eq4 = (3*lamb*omega1/2*sen_diff*(omega1-omega2))+(3*lamb*(omega1**2)*omega2/2*sen_diff)-(3*lamb*(wp**2)*np.sin(teta2))+(eq3)*3*(lamb/2)*cos_diff

    return [eq1, eq2, eq3, eq4]

# Resolvendo as equações
s_n = solve_ivp(sistema_naoLinear, [0, tf], x01, method='DOP853', dense_output=True)
x_n = s_n.sol(t)
# Extrai a primeira equação (eq1) da solução
teta1 = x_n[0]
teta2 = x_n[1]
omega1 = x_n[2]
omega2 = x_n[3]

# Plota a primeira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('Theta1')
plt.title('Gráfico de Theta1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a terceira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('omega1')
plt.title('Gráfico de omega1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a segunda equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Theta 2')
plt.title('Gráfico de Theta 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a quarta equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Omega 2')
plt.title('Gráfico de Omega 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()




# Resolvendo as equações
s_n1 = solve_ivp(sistema_naoLinear, [0, tf], x01, method='RK45', dense_output=True)
x_n1 = s_n1.sol(t)

# Extrai a primeira equação (eq1) da solução
teta1 = x_n1[0]
teta2 = x_n1[1]
omega1 = x_n1[2]
omega2 = x_n1[3]

# Plota a primeira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta1, label='Theta 1')
plt.xlabel('Tempo')
plt.ylabel('Theta 1')
plt.title('Gráfico de Theta1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a terceira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega1, label='Theta 1')
plt.xlabel('Tempo')
plt.ylabel('Omega 1')
plt.title('Gráfico de omega1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a segunda equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Theta 2')
plt.title('Gráfico de Theta 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a quarta equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Omega 2')
plt.title('Gráfico de Omega 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()


#Resolvendo as equações
s_n = solve_ivp(sistema_naoLinear, [0, tf], x02, method='DOP853', dense_output=True)
x_n = s_n.sol(t)
# Extrai a primeira equação (eq1) da solução
teta1 = x_n[0]
teta2 = x_n[1]
omega1 = x_n[2]
omega2 = x_n[3]
# Plota a primeira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('Theta1')
plt.title('Gráfico de Theta1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a terceira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('Omega 1')
plt.title('Gráfico de omega1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a segunda equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Theta 2')
plt.title('Gráfico de Theta 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a quarta equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Omega 2')
plt.title('Gráfico de Omega 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()

# Resolvendo as equações
s_n1 = solve_ivp(sistema_naoLinear, [0, tf], x02, method='DOP853', dense_output=True)
x_n1 = s_n1.sol(t)
# Extrai a primeira equação (eq1) da solução
teta1 = x_n1[0]
teta2 = x_n1[1]
omega1 = x_n1[2]
omega2 = x_n1[3]
# Plota a primeira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('Theta 1')
plt.title('Gráfico de Theta1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a terceira equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega1, label='Theta1')
plt.xlabel('Tempo')
plt.ylabel('Omega 1')
plt.title('Gráfico de omega1 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a segunda equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, teta2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Theta 2')
plt.title('Gráfico de Theta 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()
# Plota a quarta equação em função do tempo
plt.figure(figsize=(10, 6))
plt.plot(t, omega2, label='Theta 2')
plt.xlabel('Tempo')
plt.ylabel('Omega 2')
plt.title('Gráfico de Omega 2 em função do tempo')
plt.legend()
plt.grid(True)
plt.show()



