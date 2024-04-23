#to simulate the hamiltonian of graphene nano-chunk
import numpy as np
import matplotlib.pyplot as plt
import LIBRARY as lb

U = 2
t = 3
nr = 3# no. of rows of centres
N = int(2*(((nr+3)/2)**2 - 1))#no. of lattice sites



H_up = lb.H_con_dia(N, U, t, nr)
H_down = np.array(H_up)

H_up[1, 1] = U*0.4
H_down[1, 1] = U*0.6

H_up[5, 5] = U*0.4
H_down[5, 5] = U*0.6

H_up[12, 12] = U*0.4
H_down[12, 12] = U*0.6

H_up[11, 11] = U*0.4
H_down[11, 11] = U*0.6

# H_up[5, 5] = U*0.2
# H_down[5, 5] = U*0.8

# H_up[8, 8] = U*0.3
# H_down[8, 8] = U*0.7


H_f = lb.H_block(H_up, H_down)
# print(H_f)
np.set_printoptions(precision=2)
np.savetxt('h_f_dia_init_1.txt',H_f, delimiter='\t',fmt='%1.2f')
# print(np.diag(H_f))

tol = 0.01



# print(den, count, diff= scf(H_f, tol))
den, lattice_index, spinup_list, spindown_list, netspin_list, count, diff, countlist, difflist= lb.scf(H_f, tol)

den_lat = []
for i in range (0, N):
    den_lat.append([den[i] + den[i+N]])
den_lat = np.array(den_lat)



out = np.append(lattice_index, netspin_list, axis = 1)
out = np.append(out, den_lat, axis = 1)
out = np.append(out, spinup_list, axis = 1)
out = np.append(out, spindown_list, axis = 1)

plt.plot(countlist, difflist)
plt.xlabel("number of iteration")
plt.ylabel("norm of difference")
plt.show()

np.set_printoptions(precision=2)
np.savetxt('output_dia_init_1.txt',out, delimiter='\t',fmt='%1.2f', header='lattice site\tnet spin density\telectron density\tspin-up electron density\tspin-down electron density', footer = "Number of iterations = "+str(count))







    

