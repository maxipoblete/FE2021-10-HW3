import numpy as np
from quad4 import quad4
from quad4 import quad4_post
from numpy import array, pi, zeros, ix_,around
import matplotlib.pylab as plt
from scipy.linalg import solve

fid = open("2D_2.msh","r")


LINE_ELEMENT = 1
TRI_ELEMENT  = 2
QUAD_ELEMENT = 3 

Empotrado = 1
BordeNatural = 2
Placa = 3 
Extremos = 4 

#---------------------------------------------
#-------------------NODOS---------------------
#---------------------------------------------

while True:
    line = fid.readline()
    
    if line.find("$Nodes")>=0:
        break


Nnodes = int(fid.readline())   # = 8
print (f'Nnodes ={Nnodes}')


xy = np.zeros([Nnodes,2])
for i in range(Nnodes):
    line = fid.readline() 
    sl = line.split()
    xy[i,0] = float(sl[1])
    xy[i,1] = float(sl[2])


print (f' xy = {xy}')


#---------------------------------------------
#----------------CONEXIONES-------------------
#---------------------------------------------


while True:
    line = fid.readline()
    if line.find("$Elements")>=0:
        break

Nelements = int(fid.readline())  
print (f'Nelements = {Nelements}')



conec = np.zeros((Nelements,4), dtype=np.int32)
fixed_nodes = []
natural_nodes = []
Nquads= 0
Quadrangles = []
Extremos_Element = []
Placa_Element = [] 


for i in range(Nelements):
    line = fid.readline() 
    sl = line.split()
    
    
    element_number  =   np.int32(sl[0]) -1 
    element_type    =   np.int32(sl[1])
    physical_grp    =   np.int32(sl[3])
    entity_number   =   np.int32(sl[4])


    if element_type  == LINE_ELEMENT and \
        physical_grp == Empotrado: 

        n1 = np.int32(sl[5]) -1 
        n2 = np.int32(sl[6]) -1 
        fixed_nodes += [n1, n2]
        
    if element_type  == LINE_ELEMENT and \
        physical_grp == BordeNatural: 

        n1 = np.int32(sl[5]) -1 
        n2 = np.int32(sl[6]) -1 
        natural_nodes += [n1, n2]


    if element_type == QUAD_ELEMENT and \
        physical_grp == Placa or physical_grp == Extremos: 
        
        n0 = np.int32(sl[5]) -1
        n1 = np.int32(sl[6]) -1
        n2 = np.int32(sl[7]) -1
        n3 = np.int32(sl[8]) -1
        
        Quadrangles.append(element_number)
        Nquads += 1 

        conec[element_number, :] = [n0, n1, n2, n3]
        
        if physical_grp == Extremos:
            Extremos_Element.append(element_number)
            
        if physical_grp == Placa:
            Placa_Element.append(element_number)
            
            
            
        
        

print (conec)
print ("Fin del Archivo")

fixed_nodes = np.unique(fixed_nodes)
natural_nodes = np.unique(natural_nodes)
# print (f'conec = {conec}')
# print (f'fixed_nodes = {fixed_nodes}')


#---------------------------------------------
#---------------------------------------------
#----------------Ensamblaje-------------------
#---------------------------------------------
#---------------------------------------------

NDOFs_per_node = 2
NDOFs          = 2 *Nnodes




#---------------------------------------------
#----------------Restricciones----------------
#---------------------------------------------
constrained_DOFs = []

for n in fixed_nodes:
    constrained_DOFs += [2*n, 2*n +1]
    
free_DOFs = np.arange(NDOFs)
free_DOFs = np.setdiff1d(free_DOFs,constrained_DOFs)  #Resta de subconjuntos 
print (f'free_DOFs = {free_DOFs}')
print (f'constrained_DOFs = {constrained_DOFs}')


#---------------------------------------------
#----------------Propiedades------------------
#---------------------------------------------

properties = {}

ρ = 2500
g = 9.81

properties["E"]     = 20e9
properties["nu"]    = 0.25
properties["bx"]    = 0
properties["by"]    = 0
properties["t"]     = 4e-3


properties1 = {}
properties1["E"]     = 20e9
properties1["nu"]    = 0.25
properties1["bx"]    = 0
properties1["by"]    = 0
properties1["t"]     = 5e-3


print (NDOFs)


K = zeros((NDOFs, NDOFs))
f = zeros((NDOFs, 1))
u = zeros((NDOFs,1)) 


for n in natural_nodes:
    f[2*n]=1000.0





#---------------------------------------------
#----------------KKKKKKKKKKK------------------
#---------------------------------------------
for e in Quadrangles:      #range(1, Nelements):
    
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    print(f"e={e}   ni={ni}  nj={nj}   nk={nk}")

    xy_e = xy[[ni,nj,nk, nl],:]
    # print (f"xy_e = {xy_e}")
    
    
    if e in Extremos_Element:
        ke,fe = quad4(xy_e,properties1)
        # print (f"ke = {ke}" )
    
    if e in Placa_Element:
        ke,fe = quad4(xy_e,properties)
        # print (f"ke = {ke}" )
    
    
    d = [2*ni, 2*ni+1 ,2*nj, 2*nj+1 ,2*nk , 2*nk+1,2*nl , 2*nl+1]
    
    for i in range(4*NDOFs_per_node):
        p = d[i]
        for j in range(4*NDOFs_per_node):
            q = d[j]
            K[p,q] += ke[i,j]
        f[p] += fe[i]


#---------------------------------------------
#----------------Solucion---------------------
#---------------------------------------------

Kff = K[ix_(free_DOFs, free_DOFs)]
Kfc = K[ix_(free_DOFs, constrained_DOFs)]
Kcf = K[ix_(constrained_DOFs, free_DOFs)]
Kcc = K[ix_(constrained_DOFs, constrained_DOFs)] 
ff = f[free_DOFs]
fc = f[constrained_DOFs]

u[free_DOFs]=solve(Kff,ff)
R = Kcf @ u[free_DOFs] + Kcc @ u[constrained_DOFs] - fc

print (f'u = {u}')


#---------------------------------------------
#------------Grafico Viga inicial-------------
#---------------------------------------------

for e in Quadrangles:      #range(1, Nelements):
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    
    xy_e = xy[[ni, nj, nk, nl, ni],:]
    plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.plot(xy[:, 0], xy[:, 1] ,'.')  #Punto nodos


#---------------------------------------------
#------------Grafico Viga Deforma-------------
#---------------------------------------------


factor = 1e1
uv = u.reshape([-1,2])  #Reordena

for e in Quadrangles:
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    
    xy_e = xy[[ni, nj, nk,nl, ni],:] + factor*uv[[ni, nj, nk, nl, ni],:]
    plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.plot(xy[:, 0]  + factor*uv[:, 0], xy[:, 1]+ factor*uv[:, 1] ,'.')  #Punto nodos

plt.axis('equal')






from gmsh_post import write_node_data_2,write_element_data

nodes = np.arange(1,Nnodes+1)
write_node_data_2("desplazamientos.msh" , nodes, uv[:,0],uv[:,1] , "Despl")


σxx = np.zeros(Nquads+1)
σyy = np.zeros(Nquads+1)
σxy = np.zeros(Nquads+1)

i=0

for e in Quadrangles:
    
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    
    xy_e= xy[[ni, nj, nk, nl, ni],:]
    uv_e= uv[[ni, nj, nk, nl],:] 
    u_e= uv_e.reshape((-1))

    if e in Extremos_Element:
        εe, σe = quad4_post(xy_e, u_e, properties1)
        # print (f"ke = {ke}" )
    
    if e in Placa_Element:
        εe, σe = quad4_post(xy_e, u_e, properties)
        # print (f"ke = {ke}" )
        
        
    print (f"i = {i} e = {e}  σe = {σe}")
    
    σxx[i] = σe[0]
    σyy[i] = σe[1]
    σxy[i] = σe[2]
    
    i+=1 

elementos = np.array(Quadrangles) +1 

write_element_data("sx.msh", elementos, σxx , "Sigma_x" )
    
    
    
    
    
    

    












