# FE2021-10-HW3.3
### [ Group Members: Maximiliano Poblete , Lucas Raggio ] 
## Lastest Update -> For Homework 3.3, check below!

# FE2021-10-HW3.1a

Group Members: Maximiliano Poblete , Lucas Raggio

Hello, this is my HW3.1 , I used quad elements for the mesh, the dimensions of the plate are the requested. The First image is the 2D mesh, and the other images are just to see how the plate will look like (extrude not included in the .geo and .msh files)


![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Mesh2.png  )

![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Vista1.png  )

![alt text](   https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Vista%20thicc.png  )

![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Vista%203.png   )

# FE2021-10-HW3.1b

In part b of homework 3 part 1, FEM analysis was performed on a plate with different thicknesses. In the natural boundary (at the left), a force of 1 KN was applied to each node of the mesh. 

![alt text](   https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HOMEWORK_3.1b/Displacements.png  )

![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HOMEWORK_3.1b/Stresses_x.png   )





# FE2021-10-HW3.2

In this part of the task, the objective is to make 3 meshes of Quad 4 elements, coarse mesh (80 nodes), medium (950 nodes) and fine (3372 nodes), which are presented below.

![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Meshes.png   )


Below are the displacements for the 3 meshes calculated with the python script in the folder HW3.2. The units shown in the lower bar are in [m]. The amplifcation factors are 1000, 1200, and 1200 respectively

![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Displacements.png  )


Once the displacements are calculated, the x, y and xy stresses can be calculated. Below are the element stresses for the three types of meshes. The units of the stresses in the lower bar are in [Pa]. It is important to mention that in each of the stress graphs 25 intervals of iso-values were used to fill in and improve the visibility of the stresses. Unlike the previous installment, in which smoothing plug-ins were used.

<br>
<br>
<br>

### Sigma X
![alt text](  https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Sigma_X.png  )

<br>
<br>
<br>

### Sigma Y

![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Sigma_Y.png  )

<br>
<br>
<br>

### Sigma XY

![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/Sigma_XY.png  )

<br>
<br>
<br>


The maximum absolute stress component in each mesh was calculated and plotted according to the number of nodes for each one.

<br>



![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/NodesVsStresses.png  )

In the previous graph it can be seen that the more nodes the mesh has, the higher the maximum absolute stresses. Furthermore, a convergence in the maximum stresses can be seen as a function of the number of nodes. In order to estimate the convergence of the maximum stresses, it is necessary to make more types of meshes in order to have more points on the curve and the trend of the data can be analyzed with greater accuracy. Finally, the maximum values are found around the hole, therefore, a finer meshing implies that a greater detail of the stresses that are occurring in the hole due to the distributed load of 1 kN is obtained.

# FE2021-10-HW3.3
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_01.png  )

![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_02.png  )

![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_03.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_04.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_05.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_06.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_07.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_08.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_09.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_10.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_11.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_12.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_13.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_14.png  )
![alt text]( https://github.com/maxipoblete/FE2021-10-HW3/blob/main/HW3.3/Images/Homework%203%20-%20Finite%20Elements_Pa%CC%81gina_15.png  )
















