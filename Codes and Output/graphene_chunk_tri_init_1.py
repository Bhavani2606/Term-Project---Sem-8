import numpy as np
import matplotlib.pyplot as plt
import LIBRARY as lb
import math
from operator import itemgetter

U = 0.9
t = 1
nr = 3 #no. of rows
N = (nr+3)*(nr+1) - 2 #no. of lattice sites




H_up = lb.H_con_tri(N, U, t)
H_down = np.array(H_up)


H_up[2, 2] = 0.4*U
H_down[2,2] = 0.6*U

H_up[12, 12] = 0.4*U
H_down[12, 12] = 0.6*U

H_up[17, 17] = 0.4*U
H_down[17, 17] = 0.6*U

H_up[1, 1] = 0.9*U
H_down[1, 1] = 0.1*U

H_up[9, 9] = 0.8*U
H_down[9, 9] = 0.2*U

H_up[21, 21] = 0.6*U
H_down[21, 21] = 0.4*U

H_up[14, 14] = 0.9*U
H_down[14, 14] = 0.1*U

# H_up[10, 10] = 0.8*U
# H_down[10, 10] = 0.2*U

H_f = lb.H_block(H_up, H_down)
# print(H_f)
np.set_printoptions(precision=2)
np.savetxt('h_f_tri_init_1.txt',H_f, delimiter='\t',fmt='%1.2f')
# print(np.diag(H_f))
# for i in range(0, int(len(H_f)/2)):
#     for j in range (0, i):
#         if H_f[i, j]!=0:
#             print(i+1, j+1)

tol = 0.01

den, lattice_index, spinup_list, spindown_list, netspin_list, count, diff, countlist, difflist = lb.scf(H_f, tol)

den_lat = []
for i in range (0, N):
    den_lat.append([den[i] + den[i+N]])
den_lat = np.array(den_lat)
mag = np.sum(abs(netspin_list))

out = np.append(lattice_index, netspin_list, axis = 1)
out = np.append(out, den_lat, axis = 1)
out = np.append(out, spinup_list, axis = 1)
out = np.append(out, spindown_list, axis = 1)

plt.plot(countlist, difflist)
plt.xlabel("number of iteration")
plt.ylabel("norm of difference")
plt.show()
# out = np.vstack((title, out))
# print(out)

np.set_printoptions(precision=2)
np.savetxt('output_tri_init_1.txt',out, delimiter='\t',fmt='%1.2f', header='lattice site\tnet spin density\telectron density\tspin-up electron density\tspin-down electron density', footer = "Number of iterations = "+str(count)+"\n Net magnetization = "+str(mag))

#to generate a plot showing the honeycomb lattice used in this simulation and the various spin states in the lattice sites.

fig, ax = plt.subplots()
def hex(ax, nr, radius):
    points = []
    xlist = []
    for i in range (0, nr):#for ith row
        xlist1 = []
        for j in range (0, nr-i):#for the jth hexagon in the ith row
            if i == 0:
                centre_x = 1.5*radius + math.sqrt(3)*j #x coordinate centre
                centre_y = 1.5*radius +i*1.5*radius #y coordinate of centre
                xlist1.append(centre_x)
                pointsin = lb.generate_hexagon([centre_x, centre_y], radius)
                for l in pointsin:
                    points.append(l)
                lb.plot_hexagon(ax, pointsin)
            else:
                centre_x = (xlist[i-1][j] + xlist[i-1][1+j])/2
                centre_y = 1.5*radius +i*1.5*radius
                xlist1.append(centre_x)
                pointsin = lb.generate_hexagon([centre_x, centre_y], radius)
                for l in pointsin:
                    points.append(l)
                lb.plot_hexagon(ax, pointsin)
        xlist.append(xlist1)
    points = np.round_(np.array(points), decimals = 2)
    coord = np.unique(points, axis = 0)
    return coord


coords = np.array(hex(ax, 3, 1))


def lattice_index(coord):
    coord = np.array(coord)
    ind_coord = []
    m = 1
    while len(coord) != 0:
        x_max = np.max(coord[:, 0])
        # print("xmax", x_max)
        temp = []
        n_coord = len(coord)
        i_list = []
        for i in range (0, n_coord):
            if coord[i, 0] == x_max:
                # print("valid", len(coord))
                temp.append(coord[i])
                i_list.append(i)
        coord = np.delete(coord, i_list, axis = 0)
        temp1 = np.array(sorted(temp, key=itemgetter(1)))
        for k in range(len(temp1), 0, -1):
            ind_coord.append(temp1[k-1])
            m+=1
        # print(len(coord))
    return ind_coord


A = np.array(lattice_index(coords))

check = []
for p in range(len(netspin_list)):
    if netspin_list[p] < 0:
        circle = plt.Circle(A[p], abs(netspin_list[p])*2, fill=True, color='r')
        # print(p, "down")
    elif netspin_list[p] > 0:
        circle = plt.Circle(A[p], abs(netspin_list[p])*2, fill=True, color='b')
        # print(p, "up")
    else:
        continue
    ax.add_artist(circle)
    check.append([A[p], netspin_list[p]])

ax.set_aspect('equal', adjustable='box')
    
#Show the plot
plt.xlabel("X-Axis")
plt.ylabel("Y-Axis")
plt.show()
# print(check)