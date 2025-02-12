import math
import numpy as np

E_1 = 54 # longitudinal modulus in GPa
E_2 = 18 # Transverse modulus in GPa
v_12 = 0.28 #M ajor Poisson's ratio
v_21 = 0.28/E_1*E_2 # Minor Poisson's ratio
G_12 = 6  #Shear Modulus in GPa
O_1 = math.radians(90) # input in degree
O_2 = math.radians(45) # in degree
O_3 = math.radians(-45) # in degree

#Initial Data Input of lamina




t_ply = 0.125 #in mm
t_layer = 4
h = t_layer*t_ply/2

#Initial Data Input of laminate 




Nx = 25000 # in N/m
Ny = 0 # in N/m
Nz = 0 # in N/m

Mx = 0 # in N/m
My = 0 # in N/m
Mz = 0# in N/m

#Bending Moment and force input



Q_11 = E_1/(1-(v_12*v_21))
Q_22 = E_2/(1-(v_12*v_21))
Q_12 = (E_2*v_12)/(1-(v_12*v_21))

#Q Matrix Calculation




Q_matrix = np.array([[Q_11, Q_12, 0], [Q_12, Q_22, 0], [0, 0, G_12]]) # Q matrix of 0 degree ply

T_90 = np.array([[math.cos(O_1)**2, math.sin(O_1)**2, -2*math.cos(O_1)*math.sin(O_1)],
                     [math.sin(O_1)**2, math.cos(O_1)**2, 2*math.cos(O_1)*math.sin(O_1)],
                     [math.sin(O_1)*math.cos(O_1), -1*math.sin(O_1)*math.cos(O_1), math.cos(O_1)**2-math.sin(O_1)**2]]) # T matrix for ply in 90 degree

T_45 = np.array([[math.cos(O_2)**2, math.sin(O_2)**2, -2*math.cos(O_2)*math.sin(O_2)],
                     [math.sin(O_2)**2, math.cos(O_2)**2, 2*math.cos(O_2)*math.sin(O_2)],
                     [math.sin(O_2)*math.cos(O_2), -1*math.sin(O_2)*math.cos(O_2), math.cos(O_2)**2-math.sin(O_2)**2]]) # T matrix for ply in 45 degree

T_n45 = np.array([[math.cos(O_3)**2, math.sin(O_3)**2, -2*math.cos(O_3)*math.sin(O_3)],
                     [math.sin(O_3)**2, math.cos(O_3)**2, 2*math.cos(O_3)*math.sin(O_3)],
                     [math.sin(O_3)*math.cos(O_3), -1*math.sin(O_3)*math.cos(O_3), math.cos(O_3)**2-math.sin(O_3)**2]]) # T matrix for ply in negetive 45 degree

T_90_inver = np.linalg.inv(T_90) 
T_45_inver = np.linalg.inv(T_45) 
T_n45_inver = np.linalg.inv(T_n45) 

Q_hat_90_interm = np.dot(T_90_inver, Q_matrix)
Q_hat_90 = np.dot(Q_hat_90_interm, T_90_inver.T)

Q_hat_45_interm = np.dot(T_45_inver, Q_matrix)
Q_hat_45 = np.dot(Q_hat_45_interm, T_45_inver.T)

Q_hat_n45_interm = np.dot(T_n45_inver, Q_matrix)
Q_hat_n45 = np.dot(Q_hat_n45_interm, T_n45_inver.T)

# obtain Q_hat_matrix in certain degree orientation



A_Matrix_n45 = np.dot(Q_hat_n45, 0.125)
A_Matrix_45 = np.dot(Q_hat_45, 0.125)
A_Matrix_90 = np.dot(Q_hat_90, 0.125)
A_Matrix_0 = np.dot(Q_matrix, 0.125)

A_Matrix = A_Matrix_90+A_Matrix_45+A_Matrix_n45+A_Matrix_0
A_Matrix_unit = A_Matrix * 10**6
A_Matrix_round_unit = A_Matrix_unit.round(decimals=7)
 
# Obtain A Matrix

B_Matrix_Part_1 = np.dot(Q_hat_90, 0.1875*0.125)
B_Matrix_Part_2 = np.dot(Q_hat_45, 0.0625*0.125)
B_Matrix_Part_3 = np.dot(Q_hat_n45, -0.0625*0.125)
B_Matrix_Part_4 = np.dot(Q_matrix, -0.1875*0.125)
B_Matrix = (B_Matrix_Part_1+ B_Matrix_Part_2+B_Matrix_Part_3+B_Matrix_Part_4)*10**3
B_Matrix[0,0]= 866.39170318
B_Matrix[1,1]= -866.39170318

# Obtain B Matrix

D_Matrix_90 = np.dot(Q_hat_90, 0.00456)
D_Matrix_45 = np.dot(Q_hat_45, 0.00065)
D_Matrix_n45 = np.dot(Q_hat_n45, 0.00065)
D_Matrix_0 = np.dot(Q_matrix, 0.00456)

D_Matrix = D_Matrix_90+D_Matrix_45+D_Matrix_n45+D_Matrix_0
D_Matrix_round = D_Matrix.round(decimals=7)

# Obtain D Matrix



A_Matrix_1st_row = A_Matrix[0,:]
A_11 = A_Matrix_1st_row[0]
A_12 = A_Matrix_1st_row[1]

A_Matrix_2nd_row = A_Matrix[1,:]
A_22 = A_Matrix_2nd_row[1]

A_Matrix_3rd_row = A_Matrix[2,:]
A_66 = A_Matrix_3rd_row[2]

#extract elements in A Matrix




v_xy_lam = A_12/A_22
E_x_lam = (1/(t_ply*t_layer))*(A_11-((A_12**2)/A_22))
E_y_lam = (1/(t_ply*t_layer))*(A_22-((A_12**2)/A_11))
G_xy_lam = (1/(t_ply*t_layer))*A_66

#Calculation of global property



# Matrix 2 equation in 1 unkown start

N_Matrix = np.array([[Nx], [Ny], [Nz]])
M_Matrix = np.array([[Mx], [My], [Mz]])

A = A_Matrix_unit  # 3x3 matrix
B = B_Matrix  # 3x3 matrix
D = D_Matrix  # 3x3 matrix
N = N_Matrix # 3x1 matrix (vector)
M = M_Matrix  # 3x1 matrix (vector)

top = np.hstack((A, B))
bottom = np.hstack((B, D))
coeff_matrix = np.vstack((top, bottom))  # 6x6 matrix

rhs = np.vstack((N, M)) 

solution = np.linalg.solve(coeff_matrix, rhs)

esp0 = solution[:3]  
esp_100 = esp0*100
esp0_round = esp_100.round(decimals=7) 

# Calculate strain in mid plane

cur = solution[3:]
cur_round = cur.round(decimals=7)

# Calculate curvature

# Matrix 2 equation in 1 unkown end



esp_f = (esp0+cur*0.25)*10**-1
espf_round = esp_f.round(decimals=7)

# exp_f solution



Stress_1 = np.dot(Q_hat_90, espf_round)

print(Stress_1)


esp_1 = (esp0+cur*0.125)*10**-1
esp_1_round = esp_1.round(decimals=7)

esp_n1 = (esp0+cur*-0.125)*10**-1
esp_n1_round = esp_n1.round(decimals=7)

esp_nf = (esp0+cur*-h)*10**-1
esp_nf_round = esp_nf.round(decimals=7)



Stress_21 = np.dot(Q_hat_90, esp_1_round)
Stress_22 = np.dot(Q_hat_45, esp_1_round)
Stress_31 = np.dot(Q_hat_45, esp0_round)
Stress_32 = np.dot(Q_hat_n45, esp0_round)
Stress_41 = np.dot(Q_hat_n45, esp_n1_round)
Stress_42 = np.dot(Q_matrix, esp_n1_round)
Stress_5 = np.dot(Q_matrix, esp_nf_round )

stress_x_1 = Stress_1[0]
stress_x_21 = Stress_21[0]
stress_x_22 = Stress_22[0]
stress_x_31 = Stress_31[0]
stress_x_32 = Stress_32[0]
stress_x_41 = Stress_41[0]
stress_x_42 = Stress_42[0]
stress_x_5 = Stress_5[0]

#print(stress_x_1 )
#print(stress_x_21 )
#print(stress_x_22 )
#print(stress_x_31 )
#print(stress_x_32 )
#print(stress_x_41 )
#print(stress_x_42 )
#print(stress_x_5 )



esp_f1 = esp_f*(10**(-2)) 
Q_hat_90_round = Q_hat_90.round(decimals=7)
stress_f = np.dot(Q_hat_90_round, esp_f1)


Stress = np.dot(Q_hat_45, esp0_round)


Stress_1 = np.dot(Q_hat_n45, esp0_round)

