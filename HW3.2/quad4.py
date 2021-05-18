from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv,norm

def quad4(xy, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]

	ke = zeros((8,8))
	fe = zeros((8,1))

	#Primer punto de Gauss de la regla 2x2
	# xi = 1.0 / sqrt(3)
	# eta = -1.0 / sqrt(3)
	# wi = 1.0
	# wj = 1.0

	gauss_rule = [
		(-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
		( 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
		( 1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
		(-1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
	]

	for xi, eta, wi, wj in gauss_rule:

		# print(f"xi = {xi} eta = {eta}")

		x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
		y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
		dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
		dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
		dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
		dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
		
		dN0_dxi = eta/4. - 1/4.
		dN0_deta = xi/4. - 1/4.
		dN1_dxi = 1/4. - eta/4.
		dN1_deta = -xi/4. - 1/4.
		dN2_dxi = eta/4. + 1/4.
		dN2_deta = xi/4. + 1/4.
		dN3_dxi = -eta/4. - 1/4.
		dN3_deta = 1/4. - xi/4.

		# print(f"x = {x} y = {y}")

		J = array([
		[dx_dxi, dx_deta],
		[dy_dxi, dy_deta]
		]).T

		detJ = det(J)

		if detJ <= 0.:
			print(f"FATAL! detJ <= 0...")
			exit(-1)

		Jinv = inv(J)

		# print(f"J = {J}")
		# print(f"detJ = {detJ}")

		dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
		dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
		dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
		dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

		# ε = B ue
		B = zeros((3, 8))
		B[0,0] = dN0_dxy[0]
		B[1,1] = dN0_dxy[1]
		B[2,0] = dN0_dxy[1]
		B[2,1] = dN0_dxy[0]
		B[0,2] = dN1_dxy[0]
		B[1,3] = dN1_dxy[1]
		B[2,2] = dN1_dxy[1]
		B[2,3] = dN1_dxy[0]
		B[0,4] = dN2_dxy[0]
		B[1,5] = dN2_dxy[1]
		B[2,4] = dN2_dxy[1]
		B[2,5] = dN2_dxy[0]
		B[0,6] = dN3_dxy[0]
		B[1,7] = dN3_dxy[1]
		B[2,6] = dN3_dxy[1]
		B[2,7] = dN3_dxy[0]

		# print(f"B = {B}")


		ke += t * wi * wj * B.T @ Eσ @ B * detJ


	return ke, fe










def quad4_post(xy, u_e, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]


	σBIG = []
	εBIG = []

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]

	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]



	gauss_rule = [
		(-1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
		( 1.0 / sqrt(3), -1.0 / sqrt(3), 1.0, 1.0),
		( 1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
		(-1.0 / sqrt(3),  1.0 / sqrt(3), 1.0, 1.0),
	]

	for xi, eta, wi, wj in gauss_rule:



		x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
		y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
		dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
		dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
		dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
		dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
		
		dN0_dxi = eta/4. - 1/4.
		dN0_deta = xi/4. - 1/4.
		dN1_dxi = 1/4. - eta/4.
		dN1_deta = -xi/4. - 1/4.
		dN2_dxi = eta/4. + 1/4.
		dN2_deta = xi/4. + 1/4.
		dN3_dxi = -eta/4. - 1/4.
		dN3_deta = 1/4. - xi/4.


		J = array([
		[dx_dxi, dx_deta],
		[dy_dxi, dy_deta]
		]).T

		detJ = det(J)

		if detJ <= 0.:
			print(f"FATAL! detJ <= 0...")
			exit(-1)

		Jinv = inv(J)


		dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
		dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
		dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
		dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

		# ε = B ue
		B = zeros((3, 8))
		B[0,0] = dN0_dxy[0]
		B[1,1] = dN0_dxy[1]
		B[2,0] = dN0_dxy[1]
		B[2,1] = dN0_dxy[0]
		B[0,2] = dN1_dxy[0]
		B[1,3] = dN1_dxy[1]
		B[2,2] = dN1_dxy[1]
		B[2,3] = dN1_dxy[0]
		B[0,4] = dN2_dxy[0]
		B[1,5] = dN2_dxy[1]
		B[2,4] = dN2_dxy[1]
		B[2,5] = dN2_dxy[0]
		B[0,6] = dN3_dxy[0]
		B[1,7] = dN3_dxy[1]
		B[2,6] = dN3_dxy[1]
		B[2,7] = dN3_dxy[0]


		ε = B @ u_e
		σ = Eσ @ ε

		σBIG.append(σ)
		εBIG.append(ε)
		
		
		
	σxx = zeros(4)
	σxx[0]= σBIG[0][0]
	σxx[1]= σBIG[1][0]
	σxx[2]= σBIG[2][0]
	σxx[3]= σBIG[3][0]
	
	σyy = zeros(4)
	σyy[0]= σBIG[0][1]
	σyy[1]= σBIG[1][1]
	σyy[2]= σBIG[2][1]
	σyy[3]= σBIG[3][1]

	σxy = zeros(4)
	σxy[0]= σBIG[0][2]
	σxy[1]= σBIG[1][2]
	σxy[2]= σBIG[2][2]
	σxy[3]= σBIG[3][2]
	
	
	NN = zeros((4, 4))
	
	NN[0,0] = 1 + 0.5*(sqrt(3.))
	NN[0,1] = -0.5
	NN[0,2] = 1 - 0.5*(sqrt(3.))
	NN[0,3] = -0.5
	
	NN[1,0] =  -0.5
	NN[1,1] =  1 + 0.5*(sqrt(3.))
	NN[1,2] =  -0.5
	NN[1,3] =  1 - 0.5*(sqrt(3.))
	
	NN[2,0] = 1 - 0.5*(sqrt(3.))
	NN[2,1] = -0.5
	NN[2,2] = 1 + 0.5*(sqrt(3.))
	NN[2,3] = -0.5

	NN[3,0] =  -0.5
	NN[3,1] =  1 - 0.5*(sqrt(3.))
	NN[3,2] =  -0.5
	NN[3,3] =  1 + 0.5*(sqrt(3.))
	

	
	
	σXX =  NN @ σxx 
	σYY =  NN @ σyy 
	σXY =  NN @ σxy 
	
	
	σBIG = []
	
	for i in range(4):
		σBIG.append([σXX [i],σYY [i],σXY [i]])
	
	
	return εBIG, σBIG









def quad4_line_load(xy,properties_load, ni, nj):
    
    tx = properties_load["tx"]
    ty = properties_load["ty"]
    t  = properties_load["t"]
    #print(f" t ={t}")
    #print(f" tx ={tx}")

    xi = xy[ni,:]
    xj = xy[nj,:]

    L = norm(xj-xi)
    #print(f" L ={L}")
    #print(f" xi ={xi}")
    #print(f" xj ={xj}")
    
    fe = zeros((4,1))
    
    fe[0] = (t*tx*L)/2
    fe[1] = (t*ty*L)/2
    fe[2] = (t*tx*L)/2
    fe[3] = (t*ty*L)/2
    
    Suma = (t*tx*L)/2  + (t*tx*L)/2 
    #print(f" Suma ={Suma}")
    #print("")
    #print("")
    
    ke = 0
    
    # print (fe)
    
    return ke,fe
    




































#-------------------- EJEMPLO ---------------------

# xy = array([
# [1,1],
# [0,0]])

# properties = {}
# properties["E"] = 1.
# properties["nu"] = 0.25
# properties["bx"] = 0
# properties["by"] = 1.
# properties["t"] = 1.

# properties_load       = {}
# properties_load["t"]  = properties["t"]
# properties_load["tx"] = 1000 / (properties_load["t"]*4)
# properties_load["ty"] = 0


# ke, fe = quad4_line_load(xy,properties_load)




# print(f"ke = {ke}")

# fixed_dofs = [0, 1, 2, 3]
# free_dofs = [4, 5, 6, 7]

# ke_ff = ke[ix_(free_dofs, free_dofs)]
# fe_ff = array([0, -1, 0, -1])

# print(f"ke_ff = {ke_ff}")

# from scipy.linalg import solve

# u = zeros((8,1))
# uf = solve(ke_ff, fe_ff)

# print(f"uf = {uf}")











#-------------- CODIGO DE RESPALDO ---------------------




# def quad4_post_save(xy, u_e, properties):

# 	E = properties["E"]
# 	ν = properties["nu"]
# 	bx = properties["bx"]
# 	by = properties["by"]
# 	t = properties["t"]

# 	#Podemos pasarle otros valores de xi y eta
# 	if "xi" in properties:
# 		xi = properties["xi"]
# 	else:
# 		xi = 0.0

# 	if "eta" in properties:
# 		eta = properties["eta"]
# 	else:
# 		eta = 0.0

# 	Eσ = E / (1-ν**2) * array(
# 		[
# 		[1 , ν , 0       ]       ,
# 		[ν , 1 , 0       ]       ,
# 		[0 , 0 , (1-ν)/2 ]
# 		])

# 	x0 = xy[0,0]
# 	x1 = xy[1,0]
# 	x2 = xy[2,0]
# 	x3 = xy[3,0]

# 	y0 = xy[0,1]
# 	y1 = xy[1,1]
# 	y2 = xy[2,1]
# 	y3 = xy[3,1]

# 	x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
# 	y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
# 	dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
# 	dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
# 	dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
# 	dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
# 	
# 	dN0_dxi = eta/4. - 1/4.
# 	dN0_deta = xi/4. - 1/4.
# 	dN1_dxi = 1/4. - eta/4.
# 	dN1_deta = -xi/4. - 1/4.
# 	dN2_dxi = eta/4. + 1/4.
# 	dN2_deta = xi/4. + 1/4.
# 	dN3_dxi = -eta/4. - 1/4.
# 	dN3_deta = 1/4. - xi/4.

# 	# print(f"x = {x} y = {y}")

# 	J = array([
# 	[dx_dxi, dx_deta],
# 	[dy_dxi, dy_deta]
# 	]).T

# 	detJ = det(J)

# 	if detJ <= 0.:
# 		print(f"FATAL! detJ <= 0...")
# 		exit(-1)

# 	Jinv = inv(J)

# 	# print(f"J = {J}")
# 	# print(f"detJ = {detJ}")

# 	dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
# 	dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
# 	dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
# 	dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

# 	# ε = B ue
# 	B = zeros((3, 8))
# 	B[0,0] = dN0_dxy[0]
# 	B[1,1] = dN0_dxy[1]
# 	B[2,0] = dN0_dxy[1]
# 	B[2,1] = dN0_dxy[0]
# 	B[0,2] = dN1_dxy[0]
# 	B[1,3] = dN1_dxy[1]
# 	B[2,2] = dN1_dxy[1]
# 	B[2,3] = dN1_dxy[0]
# 	B[0,4] = dN2_dxy[0]
# 	B[1,5] = dN2_dxy[1]
# 	B[2,4] = dN2_dxy[1]
# 	B[2,5] = dN2_dxy[0]
# 	B[0,6] = dN3_dxy[0]
# 	B[1,7] = dN3_dxy[1]
# 	B[2,6] = dN3_dxy[1]
# 	B[2,7] = dN3_dxy[0]


# 	ε = B @ u_e
# 	σ = Eσ @ ε
# # 	print(f"B={B}")
# # 	print (f"Eσ")
# 	return ε, σ




# 	x = x0*(1 - eta)*(1 - xi)/4 + x1*(1 - eta)*(xi + 1)/4 + x2*(eta + 1)*(xi + 1)/4 + x3*(1 - xi)*(eta + 1)/4
# 	y = y0*(1 - eta)*(1 - xi)/4 + y1*(1 - eta)*(xi + 1)/4 + y2*(eta + 1)*(xi + 1)/4 + y3*(1 - xi)*(eta + 1)/4
# 	dx_dxi = -x0*(1 - eta)/4 + x1*(1 - eta)/4 + x2*(eta + 1)/4 - x3*(eta + 1)/4
# 	dx_deta = -x0*(1 - xi)/4 - x1*(xi + 1)/4 + x2*(xi + 1)/4 + x3*(1 - xi)/4
# 	dy_dxi = -y0*(1 - eta)/4 + y1*(1 - eta)/4 + y2*(eta + 1)/4 - y3*(eta + 1)/4
# 	dy_deta = -y0*(1 - xi)/4 - y1*(xi + 1)/4 + y2*(xi + 1)/4 + y3*(1 - xi)/4
# 	
# 	dN0_dxi = eta/4. - 1/4.
# 	dN0_deta = xi/4. - 1/4.
# 	dN1_dxi = 1/4. - eta/4.
# 	dN1_deta = -xi/4. - 1/4.
# 	dN2_dxi = eta/4. + 1/4.
# 	dN2_deta = xi/4. + 1/4.
# 	dN3_dxi = -eta/4. - 1/4.
# 	dN3_deta = 1/4. - xi/4.

# 	# print(f"x = {x} y = {y}")

# 	J = array([
# 	[dx_dxi, dx_deta],
# 	[dy_dxi, dy_deta]
# 	]).T

# 	detJ = det(J)

# 	if detJ <= 0.:
# 		print(f"FATAL! detJ <= 0...")
# 		exit(-1)

# 	Jinv = inv(J)

# 	# print(f"J = {J}")
# 	# print(f"detJ = {detJ}")

# 	dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
# 	dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
# 	dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
# 	dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])

# 	# ε = B ue
# 	B = zeros((3, 8))
# 	B[0,0] = dN0_dxy[0]
# 	B[1,1] = dN0_dxy[1]
# 	B[2,0] = dN0_dxy[1]
# 	B[2,1] = dN0_dxy[0]
# 	B[0,2] = dN1_dxy[0]
# 	B[1,3] = dN1_dxy[1]
# 	B[2,2] = dN1_dxy[1]
# 	B[2,3] = dN1_dxy[0]
# 	B[0,4] = dN2_dxy[0]
# 	B[1,5] = dN2_dxy[1]
# 	B[2,4] = dN2_dxy[1]
# 	B[2,5] = dN2_dxy[0]
# 	B[0,6] = dN3_dxy[0]
# 	B[1,7] = dN3_dxy[1]
# 	B[2,6] = dN3_dxy[1]
# 	B[2,7] = dN3_dxy[0]











