import numpy as np
from numpy import arange
from quad9 import quad9
from quad9 import quad9_line_load,quad9_post
from numpy import array, pi, zeros, ix_,around
import matplotlib.pylab as plt
from scipy.linalg import solve
from scipy.linalg import norm
from gmsh_post import write_node_data




fid = open("Q4/MG_Q4.msh","r")



TRI_ELEMENT  = 2
LINE_ELEMENT = 8
QUAD_ELEMENT = 10 

Empotrado    = 1
BordeNatural = 2
Placa        = 3 
Extremos     = 4 

#---------------------------------------------
#-------------------NODOS---------------------
#---------------------------------------------

while True:
    line = fid.readline()
    
    if line.find("$Nodes")>=0:
        break


Nnodes = int(fid.readline())   # = 8
# print (f'Nnodes ={Nnodes}')


xy = np.zeros([Nnodes,2])
for i in range(Nnodes):
    line = fid.readline() 
    sl = line.split()
    xy[i,0] = float(sl[1])
    xy[i,1] = float(sl[2])


# print (f' xy = {xy}')


#---------------------------------------------
#----------------CONEXIONES-------------------
#---------------------------------------------


while True:
    line = fid.readline()
    if line.find("$Elements")>=0:
        break

Nelements = int(fid.readline())  
# print (f'Nelements = {Nelements}')



conec            = np.zeros((Nelements,9), dtype=np.int32)
fixed_nodes      = []
natural_nodes    = []
Nquads           = 0
Quadrangles      = []
Extremos_Element = []
Placa_Element    = [] 


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
        natural_nodes.append([n1, n2])


    if element_type == QUAD_ELEMENT and \
        physical_grp == Placa or physical_grp == Extremos: 
        
        n0 = np.int32(sl[5])  - 1
        n1 = np.int32(sl[6])  - 1
        n2 = np.int32(sl[7])  - 1
        n3 = np.int32(sl[8])  - 1
        n4 = np.int32(sl[9])  - 1
        n5 = np.int32(sl[10]) - 1
        n6 = np.int32(sl[11]) - 1
        n7 = np.int32(sl[12]) - 1
        n8 = np.int32(sl[13]) - 1
        
        Quadrangles.append(element_number)
        Nquads += 1 

        conec[element_number, :] = [n0, n1, n2, n3, n4, n5, n6, n7, n8]
        
        if physical_grp == Extremos:
            Extremos_Element.append(element_number)
            
        if physical_grp == Placa:
            Placa_Element.append(element_number)
            
            
            
        
        

# print (conec)
print ("Fin del Archivo")

fixed_nodes = np.unique(fixed_nodes)
# natural_nodes = np.unique(natural_nodes)
# print (f'conec = {conec}')
# print (f'fixed_nodes = {fixed_nodes}')


#---------------------------------------------
#---------------------------------------------
#----------------Ensamblaje-------------------
#---------------------------------------------
#---------------------------------------------

NDOFs_per_node = 2
NDOFs          = 2 * Nnodes




#---------------------------------------------
#----------------Restricciones----------------
#---------------------------------------------
constrained_DOFs = []

for n in fixed_nodes:
    constrained_DOFs += [2*n, 2*n +1]
    
free_DOFs = np.arange(NDOFs)
free_DOFs = np.setdiff1d(free_DOFs,constrained_DOFs)  #Resta de subconjuntos 
# print (f'free_DOFs = {free_DOFs}')
# print (f'constrained_DOFs = {constrained_DOFs}')


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


# print (NDOFs)


K = zeros((NDOFs, NDOFs))
f = zeros((NDOFs, 1))
u = zeros((NDOFs,1)) 







#---------------------------------------------
#----------------KKKKKKKKKKK------------------
#---------------------------------------------
for e in Quadrangles:      #range(1, Nelements):
    
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    nm = conec[e,4]
    nn = conec[e,5]
    no = conec[e,6]
    np = conec[e,7]
    nq = conec[e,8]  #NODO CENTRAL
    
    # print(f"e={e}   ni={ni}  nj={nj}   nk={nk}")
    
    xy_e = xy[[ni,nj,nk,nl,nm,nn,no,np,nq],:]
    # xy_e = xy[[ ni, nm, nj, nn, nk , no , nl , np, nq  ],:]
    # print (f"xy_e = {xy_e}")
    
    
    if e in Extremos_Element:
        ke,fe = quad9(xy_e,properties1)
        # print (f"ke = {ke}" )
    
    if e in Placa_Element:
        ke,fe = quad9(xy_e,properties)
        # print (f"ke = {ke}" )
    
    
    d = [2*ni , 2*ni+1 ,
         2*nj , 2*nj+1 ,
         2*nk , 2*nk+1 ,
         2*nl , 2*nl+1 ,
         2*nm , 2*nm+1 ,
         2*nn , 2*nn+1 ,
         2*no , 2*no+1 ,
         2*np , 2*np+1 ,
         2*nq , 2*nq+1 ]
    
    
    for i in range( 9 * NDOFs_per_node ):
        p = d[i]
        for j in range( 9 * NDOFs_per_node ):
            q = d[j]
            K[p,q] += ke[i,j]
        # f[p] += fe[i]





properties_load       = {}
properties_load["t"]  = 5e-3
properties_load["tx"] = 1000 / (properties_load["t"]*4e-2)
properties_load["ty"] = 0


for NN in natural_nodes:
    
    ni = NN[0]
    nj = NN[1]

    xy_e = xy[[ni,nj],:]
    #print(f" ni = {ni}  nj = {nj}  xy_e={xy_e}")
    
    
    #print(f" xy = {xy}")
    fe = quad9_line_load(xy,properties_load, ni, nj)
    
    d = [2*ni,2*ni+1,2*nj,2*nj+1]
    
    for i in range(4):
        p=d[i]
        f[p]+=fe[i]
        
# print (f"f={f}")







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

# print (f'u = {u}')


#---------------------------------------------
#------------Grafico Viga inicial-------------
#---------------------------------------------

for e in Quadrangles:      #range(1, Nelements):
    
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    nm = conec[e,4]
    nn = conec[e,5]
    no = conec[e,6]
    np = conec[e,7]
    nq = conec[e,8]
    
    xy_e = xy[[ni,nm,nj,nn,nk,no,nl,np,ni],:]
    
    plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.plot(xy[:, 0], xy[:, 1] ,'.')  #Punto nodos


#---------------------------------------------
#------------Grafico Viga Deforma-------------
#---------------------------------------------


factor = 1e2
uv = u.reshape([-1,2])  #Reordena

for e in Quadrangles:
    
    ni = conec[e,0]
    nj = conec[e,1]
    nk = conec[e,2]
    nl = conec[e,3]
    nm = conec[e,4]
    nn = conec[e,5]
    no = conec[e,6]
    np = conec[e,7]
    nq = conec[e,8]
    
    # xy_e = xy[[ni, nj, nk,nl, ni],:] + factor * uv[[ni, nj, nk, nl, ni],:]
    xy_e = xy[[ni,nm,nj,nn,nk,no,nl,np,ni],:] + factor * uv[[ni,nm,nj,nn,nk,no,nl,np,ni],:]

    plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.plot(xy[:, 0]  + factor*uv[:, 0], xy[:, 1]+ factor*uv[:, 1] ,'.')  #Punto nodos

plt.axis('equal')




from gmsh_post import write_node_data_2,write_element_data
nodes = arange(1,Nnodes+1)
write_node_data_2("disp.msh" , nodes, uv[:,0],uv[:,1] , "Despl")





#---------------------------------------------
#------------TENSUONES Viga Deforma-----------
#---------------------------------------------


σxx = zeros(Nnodes)
σyy = zeros(Nnodes)
σxy = zeros(Nnodes)



i=0

Lista_con_Elementos_que_llegan_a_nodo_i =[]

for nodoi in range(Nnodes):
    Elementos_que_llegan_a_nodo_i = []
    
    for el in Quadrangles:       
        Nodos_de_elemento_e = conec[el]
    
        if nodoi in Nodos_de_elemento_e:
            Elementos_que_llegan_a_nodo_i.append(el)
            
    Lista_con_Elementos_que_llegan_a_nodo_i.append(Elementos_que_llegan_a_nodo_i)
            


for i,j in enumerate(Lista_con_Elementos_que_llegan_a_nodo_i):
    
    Cantidad = len(j)
    
    σxxi = 0
    σyyi = 0 
    σxyi = 0 
    
    for e in j:
               
        nodos_de_e = conec[e]
        list1 = nodos_de_e.tolist()
        Indice = list1.index(i) 
               
        ni = conec[e,0]
        nj = conec[e,1]
        nk = conec[e,2]
        nl = conec[e,3]
        nm = conec[e,4]
        nn = conec[e,5]
        no = conec[e,6]
        np = conec[e,7]
        nq = conec[e,8]      
        
        
        xy_e = xy[[ni,nj,nk,nl,nm,nn,no,np,nq],:]
        uv_e = uv[[ni,nj,nk,nl,nm,nn,no,np,nq],:] 


        u_e  = uv_e.reshape((-1))
    
        if e in Extremos_Element:
            εe, σe = quad9_post(xy_e, u_e, properties1)
        
        if e in Placa_Element:
            εe, σe = quad9_post(xy_e, u_e, properties)
    
    
        Array_σe = array(σe)
        
        Extraccion_σxx = Array_σe[:, 0] 
        Extraccion_de_interes_σxx = Extraccion_σxx[Indice]
        
        Extraccion_σyy = Array_σe[:, 1] 
        Extraccion_de_interes_σyy = Extraccion_σyy[Indice]
        
        Extraccion_σxy = Array_σe[:, 2] 
        Extraccion_de_interes_σxy = Extraccion_σxy[Indice]
        
        σxxi += Extraccion_de_interes_σxx
        σyyi += Extraccion_de_interes_σyy
        σxyi += Extraccion_de_interes_σxy
        
    Promedio_σxxi = σxxi/Cantidad
    Promedio_σyyi = σyyi/Cantidad 
    Promedio_σxyi = σxyi/Cantidad 
    σxx[i] = Promedio_σxxi
    σyy[i] = Promedio_σyyi
    σxy[i] = Promedio_σxyi
    
    

elementos = array(Quadrangles) + 1 
nodes = arange(1,Nnodes+1)

write_node_data("sxx.msh", nodes, σxx, "Sigma_x")
write_node_data("syy.msh", nodes, σyy, "Sigma_y")
write_node_data("sxy.msh", nodes, σxy, "Sigma_xy")
















# write_node_data("uy.msh", nodes, uv[:,1], "Despl. Y")
# write_node_data_2("desplazamientos_viga.msh", nodes, uv[:,0], uv[:,1], "Desplazamientos δFEM Viga")

#write_tension_Nodes_data("sx.msh", σxx , σxx , "Sigma_x")
    




    
    
       
    
    
    





# for i in Quadrangles:
    
    
 
    
#     ni = conec[e,0]
#     nj = conec[e,1]
#     nk = conec[e,2]
#     nl = conec[e,3]
    
    
    
#     xy_e = xy[[ni, nj, nk, nl],:]
#     uv_e = uv[[ni, nj, nk, nl],:] 
#     u_e  = uv_e.reshape((-1))

#     if e in Extremos_Element:
#         εe, σe = quad4_post(xy_e, u_e, properties1)
    
#     if e in Placa_Element:
#         εe, σe = quad4_post(xy_e, u_e, properties)



#     i+=1 

# for i in (σe):
#     print(i)



