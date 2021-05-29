# $NodeData
# number-of-string-tags
# < "string-tag" >
# …
# number-of-real-tags
# < real-tag >
# …
# number-of-integer-tags
# < integer-tag >
# …
# node-number value …
# …
# $EndNodeData


def write_node_data(fname,nodes,data,nombre_datos="\n"):
    
    Ndata = len(nodes)
    
    fid = open(fname,"w")
    
    
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")   
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")
    fid.write("3\n")
    fid.write("0\n")
    fid.write("1\n")            #DIMENSION DE LOS DATOS
    fid.write(f"{Ndata}\n")     #NUMERO DE DATOS
    
     
    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data[i]}\n")  
    fid.write("$EndNodeData")
    
    fid.close()
    return 
    
    
    
    
    
    
def write_node_data_2(fname, nodes , data1 , data2 , nombre_datos="\n"):
    
    Ndata = len(nodes)
    
    fid = open(fname,"w")
    
    
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")   
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")
    fid.write("3\n")
    fid.write("0\n")
    fid.write("3\n")            #DIMENSION DE LOS DATOS
    fid.write(f"{Ndata}\n")     #NUMERO DE DATOS
    
     
    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data1[i]} {data2[i]} 0.0 \n")  
    fid.write("$EndNodeData")
    
    fid.close()
    return     
    
    

def write_element_data(fname, elements , data , nombre_datos="\n"):
    
    Ndata = len(elements)
    
    fid = open(fname,"w")
    
    
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$ElementData\n")   
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")
    fid.write("3\n")
    fid.write("0\n")
    fid.write("1\n")            #DIMENSION DE LOS DATOS
    fid.write(f"{Ndata}\n")     #NUMERO DE DATOS
    
     
    for i in range(Ndata):
        fid.write(f"{elements[i]} {data[i]}\n")  
    fid.write("$EndElementData")
    
    fid.close()
    return     
    


    

def write_tension_Nodes_data(fname, elements , data , nombre_datos="\n"):
    
    Ndata = len(elements)
    
    fid = open(fname,"w")
    
    
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$ElementData\n")   
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")
    fid.write("3\n")
    fid.write("0\n")
    fid.write("1\n")            #DIMENSION DE LOS DATOS
    fid.write(f"{Ndata}\n")     #NUMERO DE DATOS
    
     
    for i in range(Ndata):
        fid.write(f"{i+1} {data[i]}\n")  
    fid.write("$EndElementData")
    
    fid.close()
    return   
    
   
    
    
        

def write_tension_Nodes_data1(fname,nodes,data,nombre_datos="\n"):
    
    Ndata = len(nodes)
    
    fid = open(fname,"w")
    
    
    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")   
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")
    fid.write("3\n")
    fid.write("0\n")
    fid.write("1\n")            #DIMENSION DE LOS DATOS
    fid.write(f"{Ndata}\n")     #NUMERO DE DATOS
    
     
    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data[i]}\n")  
    fid.write("$EndNodeData")
    
    fid.close()
    return 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    