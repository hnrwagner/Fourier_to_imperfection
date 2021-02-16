# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 09:51:16 2019

###############################################################################
###############################################################################

@author: H.N.Ronald Wagner 

###############################################################################
###############################################################################

@description

###############################################################################

This python script calculates the mean and variance imperfection pattern from the input imperfection pattern 
AND the corresponding nodes.files for ABAQUS for an First-Order Second-Moment Analysis (FOSM).

Make sure that this python script is in the same dictonary as all the input data (Fourier Coefficient - "data.FC_A" & "data.FC_B")
If all data are in the correct dictonary just "Run" the script in a:

Python Development Enviroment for example "Spyder" which is part of the Anaconda Distribution as Python Science Platform from the website:

https://www.anaconda.com/distribution/

###############################################################################
###############################################################################

@input

###############################################################################

Name of the shell identifier - myShellName

Cylinder geometry:
Length - myLength
Radius - myRadius 
Wall Thickness - myThickness

The number of nodes in axial and circumferential direction

na
nc

The name of the input data (Fourier Coefficient - "data.FC_A" & "data.FC_B")

Shells

###############################################################################
###############################################################################

@output

###############################################################################

Text files which contain the nodes for ABAQUS for the given elements file ("elements.txt")

Files for:
    
mean imperfection pattern

variance imperfection pattern

figures of the imperfection pattern

###############################################################################
###############################################################################


"""

import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import matrix_rank

###############################################################################
    
# Cylinder Geometry Identifier

myShellName = "N"

# File Names for the analysis 

myFileName = str(myShellName)+"_mean_Nodes.txt"
myImageName = str(myShellName)+"_mean_Shell.png"

# Length of the cylinder

myLength = 196.85

# Radius of the cylinder

myRadius = 101.6

# Wall thickness of the cylinder

myThickness = 0.1

# Number of Nodes in axial and circumferential direction

na = 61

nc = 240


# Displacement field

w = np.zeros((na,nc))

xyz = np.zeros((na*nc,3))

# Fourier approach (1 - phase shift, 2 - cos , 3 - sin)

case = 2
    
###############################################################################

myFC_A_v = []
myFC_B_v = []


#Shells = ["Z07","Z08","Z10","Z11","Z12"]

Shells = ["N_6","N_9","N_11"]

#Shells = ["ST_1","ST_2","ST_3","ST_4","ST_5","ST_6"]

for i in Shells:
    A = np.loadtxt(str(i)+"_measured.FC_A")
    myFC_A_v.append(A)
    B = np.loadtxt(str(i)+"_measured.FC_B")
    myFC_B_v.append(B)


n1 = len(A)
n2 = len(A[0])

n_imp = len(myFC_A_v)

phim = np.zeros((n_imp,n1,n2))

eps = np.zeros((n_imp,n1,n2))

phim_list = np.zeros(n_imp)
X_merge = []
"""
###############################################################################

@equation

The following equations are described in

"Probabilistic and deterministic lower-bound design benchmarks for cylindrical shells under axial compression"

Thin-Walled Structures
Volume 146, January 2020, 106451
###############################################################################
###############################################################################
"""
if case == 1:
    for j in range(0,n_imp,1):
        for x in range(0,n1,1):
            for y in range(0,n2,1):
                #print(B_t[x][y])
                if (myFC_A_v[j][x][y] >= 0):
                    phim[j][x][y] = np.arctan(myFC_B_v[j][x][y]/myFC_A_v[j][x][y])
                elif (myFC_A_v[j][x][y] <= 0):
                    phim[j][x][y] = np.arctan(myFC_A_v[j][x][y]/myFC_B_v[j][x][y]) + np.pi
                elif (myFC_A_v[j][x][y] == 0):
                    phim[j][x][y] = np.sign(myFC_B_v[j][x][y])*np.pi/2.0
                        
        #                
        for x in range(0,n1,1):
            for y in range(0,n2,1):
                eps[j][x][y] = np.sqrt(myFC_A_v[j][x][y] * myFC_A_v[j][x][y] + myFC_B_v[j][x][y] * myFC_B_v[j][x][y]) 
        
        ###############################################################################
        
        phim_list = phim[j]
        phim_list = phim_list.ravel()
        phim_list = phim_list.tolist()
        eps_list = eps[j]
        eps_list = eps_list.ravel()
        eps_list = eps_list.tolist()
        X = phim_list + eps_list
        X = np.nan_to_num(X) 
        X_merge.append(X)
elif case == 1 or case == 2:
    for j in range(0,n_imp,1):
        for x in range(0,n1,1):
            for y in range(0,n2,1):
                phim[j][x][y] = myFC_A_v[j][x][y]
                eps[j][x][y] = myFC_B_v[j][x][y] 
    
    ###############################################################################
    
        phim_list = phim[j]
        phim_list = phim_list.ravel()
        phim_list = phim_list.tolist()
        eps_list = eps[j]
        eps_list = eps_list.ravel()
        eps_list = eps_list.tolist()
        X = phim_list + eps_list
        X = np.nan_to_num(X) 
        X_merge.append(X)
            
    
    
    ###############################################################################
    
    # Calculate mean vector
    
    ###############################################################################
    
r = len(X_merge[0])
m = len(X_merge)
mean_v =   np.zeros((r))  

k = 0

for i in range(0,r,1):
    for j in range(0,m,1):
        k = X_merge[j][i] + k
    mean_v[i] = k
    mean_v[i] = 1/m*mean_v[i]
    k = 0
    
    #mean_v2 = np.nan_to_num(mean_v)    
    
    ###############################################################################
    
    # Calculate displacement field
    
    ###############################################################################
    
phim_mean = np.zeros((n1,n2))

eps_mean = np.zeros((n1,n2))


phim_mean = mean_v[:len(mean_v)//2]
phim_mean = np.reshape(phim_mean, (n1,n2))
#phim_mean = np.transpose(phim_mean)
eps_mean = mean_v[len(mean_v)//2:]
eps_mean = np.reshape(eps_mean, (n1,n2))
#eps_mean = np.transpose(eps_mean)



for i in range(0,na,1):
    for j in range(0,nc,1):
        x = myLength*(i)/(na)
        y = 2.0*np.pi*myRadius*(j)/(nc)

        if (case == 1):
            for k in range(0,n1,1):
                for l in range(0,n2,1):
                    w[i][j] = w[i][j] + eps_mean[k][l]*np.cos((k)*np.pi*x/myLength)*np.cos((l)*y/myRadius-phim_mean[k][l])
                    #w[i][j] = w[i][j] + eps_mean[k][l]
                    
        elif (case == 2):
            for k in range(0,n1,1):
                for l in range(0,n2,1):
                    w[i][j] = w[i][j] + np.cos((k)*np.pi*x/myLength)*phim_mean[k][l]*np.cos((l)*y/myRadius) +eps_mean[k][l]*np.sin((l)*y/myRadius)
    
        elif (case == 3):
            for k in range(0,n1,1):
                for l in range(0,n2,1):
                    w[i][j] = w[i][j] + np.sin((k)*np.pi*x/myLength)*phim_mean[k][l]*np.cos((l)*y/myRadius) +eps_mean[k][l]*np.sin((l)*y/myRadius)
 
       
        w[i][j] = w[i][j] * myThickness
        xyz[(i)*nc+j][0] = (myRadius-w[i][j])*np.cos(y/myRadius)
        xyz[(i)*nc+j][1] = (myRadius-w[i][j])*np.sin(y/myRadius)
        xyz[(i)*nc+j][2] = x*1.016666671
##    
##################################################################################
##   
        
myFileName 
myImageName       
nm = len(xyz)

ABAQUS_NODES = np.zeros((nm,4))
#
for i in range(0,nm,1):

    ABAQUS_NODES[i][0] = i+1
    ABAQUS_NODES[i][1] = xyz[i][0]
    ABAQUS_NODES[i][2] = xyz[i][1]
    ABAQUS_NODES[i][3] = xyz[i][2]
    

np.savetxt(myFileName,ABAQUS_NODES,fmt ='%i, %10.5f, %10.5f, %10.5f')

print("\n")
print("###############################################################################")
print("End of Calculation Number 1")
print("###############################################################################")

cmap = cm.get_cmap('jet')
p = plt.pcolor(w, cmap = cmap)
cb = plt.colorbar(p, orientation = 'vertical')
plt.xlabel("Circumferential Nodes",fontsize=14)
plt.ylabel("Axial Nodes",fontsize=14)
#plt.xticks(np.arange(0,405,45))
cb.set_label('Radial Displacement [mm]',fontweight='bold',fontsize=14)
plt.savefig(myImageName, dpi = 1000)
cb.remove()
p.remove()



#############################################################################

# Calculate Covariance Matrix COV

#############################################################################


COV = np.zeros((r,r))
summe = 0

for i in range(0,r,1):
    for j in range(0,r,1):
        for k in range(0,m,1):
            summe = ((X_merge[k][i]-mean_v[i])*(X_merge[k][j]-mean_v[j]))+ summe
        COV[i][j] = summe
        COV[i][j] = 1/(m-1)*COV[i][j]
        summe = 0
        

##############################################################################

# Calculate eigen values e, eigen vectors v and rank of Covariance Matrix COV

##############################################################################

e, v = LA.eig(COV)

rank = matrix_rank(COV)
    
##############################################################################

# Define spectral matrix D (eigen values)

##############################################################################

ra = np.zeros(rank)  

for i in range(0,rank,1):
    ra[i] = e[rank-1-i]
      
D = np.zeros((rank,rank))    

D = np.diag((np.sqrt(ra)))
                
##############################################################################

# Define eigen vector matrix Q

##############################################################################

Q = np.zeros((r,rank))      
    
for i in range(0,rank,1):
    for j in range(0,r,1):
        Q[j][i] = v[j][rank-1-i]

##############################################################################

# Calculate transformation matrix B_m

##############################################################################

B_m = np.matmul(Q,D)

##############################################################################

# Define standard deviation parameter Sigma for transformation

##############################################################################

Sigma = 1.5

##############################################################################

# Perform Mahalanbis Transformation

##############################################################################

x_v = np.zeros((r,rank*2))
z_v = np.zeros((rank,rank*2))

k = 0
for i in range(0,rank,1):
    z_v[i,k] = Sigma
    z_v[i,k+rank] = -Sigma
    k=k+1

for i in range(0,2*rank,1):
    x_v[:,i] = B_m.dot(z_v[:,i])+mean_v

##############################################################################

# Calculate Node files for Variance Approximation

##############################################################################


kk = 0
for k in range(0,2*rank,1):
    myFileName2 = str(myShellName)+"_var_"+str(kk+1)+"_Nodes_.txt"
    myImageName2 = str(myShellName)+"_var_"+str(kk+1)+"_Nodes_.png"
    
    phim_mean = np.zeros((n1,n2))
    
    eps_mean = np.zeros((n1,n2))
    
    
    phim_mean = x_v[:,k][:len(x_v[:,k])//2]
    phim_mean = np.reshape(phim_mean, (n1,n2))
    #phim_mean = np.transpose(phim_mean)
    eps_mean = x_v[:,k][len(x_v[:,k])//2:]
    eps_mean = np.reshape(eps_mean, (n1,n2))
    #eps_mean = np.transpose(eps_mean)
    

    for i in range(0,na,1):
        for j in range(0,nc,1):
            x = myLength*(i)/(na)
            y = 2.0*np.pi*myRadius*(j)/(nc)
    
            if (case == 1):
                for k in range(0,n1,1):
                    for l in range(0,n2,1):
                        w[i][j] = w[i][j] + eps_mean[k][l]*np.cos((k)*np.pi*x/myLength)*np.cos((l)*y/myRadius-phim_mean[k][l])
                        #w[i][j] = w[i][j] + eps_mean[k][l]
                        
            elif (case == 2):
                for k in range(0,n1,1):
                    for l in range(0,n2,1):
                        w[i][j] = w[i][j] + np.cos((k)*np.pi*x/myLength)*phim_mean[k][l]*np.cos((l)*y/myRadius) +eps_mean[k][l]*np.sin((l)*y/myRadius)
        
            elif (case == 3):
                for k in range(0,n1,1):
                    for l in range(0,n2,1):
                        w[i][j] = w[i][j] + np.sin((k)*np.pi*x/myLength)*phim_mean[k][l]*np.cos((l)*y/myRadius) +eps_mean[k][l]*np.sin((l)*y/myRadius)
     
           
            w[i][j] = w[i][j] * myThickness
            xyz[(i)*nc+j][0] = (myRadius-w[i][j])*np.cos(y/myRadius)
            xyz[(i)*nc+j][1] = (myRadius-w[i][j])*np.sin(y/myRadius)
            xyz[(i)*nc+j][2] = x*1.016666671

    ##    
    ##################################################################################
    ##    
    nm = len(xyz)
    
    ABAQUS_NODES = np.zeros((nm,4))
    #
    for i in range(0,nm,1):
    
        ABAQUS_NODES[i][0] = i+1
        ABAQUS_NODES[i][1] = xyz[i][0]
        ABAQUS_NODES[i][2] = xyz[i][1]
        ABAQUS_NODES[i][3] = xyz[i][2]
        
    
    np.savetxt(myFileName2,ABAQUS_NODES,fmt ='%i, %10.5f, %10.5f, %10.5f')
    
    print("\n")  
    cmap = cm.get_cmap('jet')
    p = plt.pcolor(w, cmap = cmap)
    cb = plt.colorbar(p, orientation = 'vertical')
    plt.xlabel("Circumferential Nodes",fontsize=14)
    plt.ylabel("Axial Nodes",fontsize=14)
    plt.xticks(np.arange(0,405,45))
    cb.set_label('Radial Displacement [mm]',fontweight='bold',fontsize=14)
    plt.savefig(myImageName2, dpi = 1000)
    cb.remove()
    p.remove()
    print("###############################################################################")
    print("End of Calculation Number "+str(kk+2)+" out of "+str(2*rank+1))
    print("###############################################################################")
    kk = kk+1
    

