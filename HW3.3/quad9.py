from matplotlib.pylab import *
from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv,norm



def quad9(xy, properties):
    
#      E = properties["E"]
#     ν = properties["nu"]
#     bx = properties["bx"]
#     by = properties["by"]
#     t = properties["t"]
    
    
#     α = 1     
#     E_L = E
#     E_T = E*α    
#     ν_LT = ν
#     ν_TL = ν*α   
#     G_LT = E_L/(2*(1+ν_LT))
    
    
    
    

#     Eσ = 1 / (1-ν_LT*ν_TL) * array(
# 		[
# 		[E_L        ,   ν_TL*E_L   , 0       ]       ,
# 		[ν_LT*E_T   ,     E_T      , 0       ]       ,
# 		[0          , 0            , G_LT/(1-ν_LT*ν_TL) ]
# 		])

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
    x4 = xy[4,0]
    x5 = xy[5,0]
    x6 = xy[6,0]
    x7 = xy[7,0]
    x8 = xy[8,0]

    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]
    y4 = xy[4,1]
    y5 = xy[5,1]
    y6 = xy[6,1]
    y7 = xy[7,1]
    y8 = xy[8,1]

    ke = zeros((18,18))
    fe = zeros((18,1))
    
    
    gauss_rule = [
        ( -sqrt(3/5) , -sqrt(3/5) , 1.0, 25/81),
        ( 0          , -sqrt(3/5) , 1.0, 40/81),
        (  sqrt(3/5) , -sqrt(3/5) , 1.0, 25/81),
        ( -sqrt(3/5) , 0          , 1.0, 40/81),
        ( 0          , 0          , 1.0, 64/81),
        (  sqrt(3/5) , 0          , 1.0, 40/81),
        ( -sqrt(3/5) , +sqrt(3/5) , 1.0, 25/81),
        ( 0          , +sqrt(3/5) , 1.0, 40/81),
        (  sqrt(3/5) ,  sqrt(3/5) , 1.0, 25/81),    
    ]
    
    
    # N0 = (Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4
    # N1 = (Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4
    # N2 = (Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4
    # N3 = (Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4
    
    # N4 = (1 - Ξ**2 )*(ζ - 1) * ζ / 2
    # N5 = (1 - ζ**2 )*(Ξ + 1) * Ξ / 2
    # N6 = (1 - Ξ**2 )*(ζ + 1) * ζ / 2
    # N7 = (1 - ζ**2 )*(Ξ - 1) * Ξ / 2
    # N8 = (1 - ζ**2 )*(1 - Ξ**2 )
    
    
    for Ξ, ζ, wi, wj in gauss_rule:
    
        x = x0*((Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4)    + x1*((Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4)  + x2*((Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4) + x3*((Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4) + x4*((1 - Ξ**2 )*(ζ - 1) * ζ / 2) + x5*((1 - ζ**2 )*(Ξ + 1) * Ξ / 2) + x6*((1 - Ξ**2 )*(ζ + 1) * ζ / 2) + x7*((1 - ζ**2 )*(Ξ - 1) * Ξ / 2) + x8*((1 - ζ**2 )*(1 - Ξ**2 ))     
        y = y0*((Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4)    + y1*((Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4)  + y2*((Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4) + y3*((Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4) + y4*((1 - Ξ**2 )*(ζ - 1) * ζ / 2) + y5*((1 - ζ**2 )*(Ξ + 1) * Ξ / 2) + y6*((1 - Ξ**2 )*(ζ + 1) * ζ / 2) + y7*((1 - ζ**2 )*(Ξ - 1) * Ξ / 2) + y8*((1 - ζ**2 )*(1 - Ξ**2 ))    
        
        
        dN0_dΞ = Ξ*ζ*(ζ - 1)/4 + ζ*(Ξ - 1)*(ζ - 1)/4
        dN0_dζ = Ξ*ζ*(Ξ - 1)/4 + Ξ*(Ξ - 1)*(ζ - 1)/4
        
        dN1_dΞ = Ξ*ζ*(ζ - 1)/4 + ζ*(Ξ + 1)*(ζ - 1)/4
        dN1_dζ = Ξ*ζ*(Ξ + 1)/4 + Ξ*(Ξ + 1)*(ζ - 1)/4
        
        dN2_dΞ = Ξ*ζ*(ζ + 1)/4 + ζ*(Ξ + 1)*(ζ + 1)/4
        dN2_dζ = Ξ*ζ*(Ξ + 1)/4 + Ξ*(Ξ + 1)*(ζ + 1)/4
        
        dN3_dΞ = Ξ*ζ*(ζ + 1)/4 + ζ*(Ξ - 1)*(ζ + 1)/4
        dN3_dζ = Ξ*ζ*(Ξ - 1)/4 + Ξ*(Ξ - 1)*(ζ + 1)/4
        
        dN4_dΞ = -Ξ*ζ*(ζ - 1)
        dN4_dζ = ζ*(1 - Ξ**2)/2 + (1/2 - Ξ**2/2)*(ζ - 1)
        
        dN5_dΞ = Ξ*(1 - ζ**2)/2 + (1 - ζ**2)*(Ξ/2 + 1/2)
        dN5_dζ = -Ξ*ζ*(Ξ + 1)
        
        dN6_dΞ = -Ξ*ζ*(ζ + 1)
        dN6_dζ = ζ*(1 - Ξ**2)/2 + (1 - Ξ**2)*(ζ/2 + 1/2)
        
        dN7_dΞ = Ξ*(1 - ζ**2)/2 + (1/2 - ζ**2/2)*(Ξ - 1)
        dN7_dζ = -Ξ*ζ*(Ξ - 1)
        
        dN8_dΞ = -2*Ξ*(1 - ζ**2)
        dN8_dζ = -2*ζ*(1 - Ξ**2)
    
    
        dx_dΞ = x0*(dN0_dΞ) + x1*(dN1_dΞ) + x2*(dN2_dΞ) + x3*(dN3_dΞ) + x4*(dN4_dΞ) + x5*(dN5_dΞ) + x6*(dN6_dΞ) + x7*(dN7_dΞ) +x8*(dN8_dΞ)
        dx_dζ = x0*(dN0_dζ) + x1*(dN1_dζ) + x2*(dN2_dζ) + x3*(dN3_dζ) + x4*(dN4_dζ) + x5*(dN5_dζ) + x6*(dN6_dζ) + x7*(dN7_dζ) +x8*(dN8_dζ)
        
        dy_dΞ = y0*(dN0_dΞ) + y1*(dN1_dΞ) + y2*(dN2_dΞ) + y3*(dN3_dΞ) + y4*(dN4_dΞ) + y5*(dN5_dΞ) + y6*(dN6_dΞ) + y7*(dN7_dΞ) +y8*(dN8_dΞ)
        dy_dζ = y0*(dN0_dζ) + y1*(dN1_dζ) + y2*(dN2_dζ) + y3*(dN3_dζ) + y4*(dN4_dζ) + y5*(dN5_dζ) + y6*(dN6_dζ) + y7*(dN7_dζ) +y8*(dN8_dζ)
        
        
        J = array([
        [dx_dΞ, dx_dζ],
        [dy_dΞ, dy_dζ]
        ]).T

        detJ = det(J)
        
        if detJ <= 0.:
            print(f"FATAL! detJ <= 0...")
            break
    
        Jinv = inv(J)
    
        dN0_dxy = Jinv @ array([ dN0_dΞ, dN0_dζ ])
        dN1_dxy = Jinv @ array([ dN1_dΞ, dN1_dζ ])
        dN2_dxy = Jinv @ array([ dN2_dΞ, dN2_dζ ])
        dN3_dxy = Jinv @ array([ dN3_dΞ, dN3_dζ ])
        dN4_dxy = Jinv @ array([ dN4_dΞ, dN4_dζ ])
        dN5_dxy = Jinv @ array([ dN5_dΞ, dN5_dζ ])
        dN6_dxy = Jinv @ array([ dN6_dΞ, dN6_dζ ])
        dN7_dxy = Jinv @ array([ dN7_dΞ, dN7_dζ ])
        dN8_dxy = Jinv @ array([ dN8_dΞ, dN8_dζ ])

        B = zeros((3,18))
        
        B[0,0]  = dN0_dxy[0]
        B[1,1]  = dN0_dxy[1]
        B[2,0]  = dN0_dxy[1]
        B[2,1]  = dN0_dxy[0]
        
        B[0,2]  = dN1_dxy[0]
        B[1,3]  = dN1_dxy[1]
        B[2,2]  = dN1_dxy[1]
        B[2,3]  = dN1_dxy[0]
        
        B[0,4]  = dN2_dxy[0]
        B[1,5]  = dN2_dxy[1]
        B[2,4]  = dN2_dxy[1]
        B[2,5]  = dN2_dxy[0]
        
        B[0,6]  = dN3_dxy[0]
        B[1,7]  = dN3_dxy[1]
        B[2,6]  = dN3_dxy[1]
        B[2,7]  = dN3_dxy[0]
        
        B[0,8]  = dN4_dxy[0]
        B[1,9]  = dN4_dxy[1]
        B[2,8]  = dN4_dxy[1]
        B[2,9]  = dN4_dxy[0]
        
        B[0,10] = dN5_dxy[0]
        B[1,11] = dN5_dxy[1]
        B[2,10] = dN5_dxy[1]
        B[2,11] = dN5_dxy[0]

        B[0,12] = dN6_dxy[0]
        B[1,13] = dN6_dxy[1]
        B[2,12] = dN6_dxy[1]
        B[2,13] = dN6_dxy[0]

        B[0,14] = dN7_dxy[0]
        B[1,15] = dN7_dxy[1]
        B[2,14] = dN7_dxy[1]
        B[2,15] = dN7_dxy[0]
        
        B[0,16] = dN8_dxy[0]
        B[1,17] = dN8_dxy[1]
        B[2,16] = dN8_dxy[1]
        B[2,17] = dN8_dxy[0]
        
        
        
        ke += t * wi * wj * B.T @ Eσ @ B * detJ
        
        # print (f"\n\n ke={ke}")
        # print (f"\n\n B={B}")

        # print(f"LARGO dN0_dxy  = {len(dN0_dxy)}" )
    
    return ke, fe





def quad9_post(xy, u_e, properties):

    
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
    x4 = xy[4,0]
    x5 = xy[5,0]
    x6 = xy[6,0]
    x7 = xy[7,0]
    x8 = xy[8,0]

    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    y3 = xy[3,1]
    y4 = xy[4,1]
    y5 = xy[5,1]
    y6 = xy[6,1]
    y7 = xy[7,1]
    y8 = xy[8,1]

    
    gauss_rule = [
        ( -sqrt(3/5) , -sqrt(3/5) , 1.0, 25/81),
        ( 0          , -sqrt(3/5) , 1.0, 40/81),
        (  sqrt(3/5) , -sqrt(3/5) , 1.0, 25/81),
        ( -sqrt(3/5) , 0          , 1.0, 40/81),
        ( 0          , 0          , 1.0, 64/81),
        (  sqrt(3/5) , 0          , 1.0, 40/81),
        ( -sqrt(3/5) , +sqrt(3/5) , 1.0, 25/81),
        ( 0          , +sqrt(3/5) , 1.0, 40/81),
        (  sqrt(3/5) ,  sqrt(3/5) , 1.0, 25/81),    
    ]

    
    for Ξ, ζ, wi, wj in gauss_rule:
    
        x = x0*((Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4)    + x1*((Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4)  + x2*((Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4) + x3*((Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4) + x4*((1 - Ξ**2 )*(ζ - 1) * ζ / 2) + x5*((1 - ζ**2 )*(Ξ + 1) * Ξ / 2) + x6*((1 - Ξ**2 )*(ζ + 1) * ζ / 2) + x7*((1 - ζ**2 )*(Ξ - 1) * Ξ / 2) + x8*((1 - ζ**2 )*(1 - Ξ**2 ))     
        y = y0*((Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4)    + y1*((Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4)  + y2*((Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4) + y3*((Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4) + y4*((1 - Ξ**2 )*(ζ - 1) * ζ / 2) + y5*((1 - ζ**2 )*(Ξ + 1) * Ξ / 2) + y6*((1 - Ξ**2 )*(ζ + 1) * ζ / 2) + y7*((1 - ζ**2 )*(Ξ - 1) * Ξ / 2) + y8*((1 - ζ**2 )*(1 - Ξ**2 ))    
        
        
        dN0_dΞ = Ξ*ζ*(ζ - 1)/4 + ζ*(Ξ - 1)*(ζ - 1)/4
        dN0_dζ = Ξ*ζ*(Ξ - 1)/4 + Ξ*(Ξ - 1)*(ζ - 1)/4
        
        dN1_dΞ = Ξ*ζ*(ζ - 1)/4 + ζ*(Ξ + 1)*(ζ - 1)/4
        dN1_dζ = Ξ*ζ*(Ξ + 1)/4 + Ξ*(Ξ + 1)*(ζ - 1)/4
        
        dN2_dΞ = Ξ*ζ*(ζ + 1)/4 + ζ*(Ξ + 1)*(ζ + 1)/4
        dN2_dζ = Ξ*ζ*(Ξ + 1)/4 + Ξ*(Ξ + 1)*(ζ + 1)/4
        
        dN3_dΞ = Ξ*ζ*(ζ + 1)/4 + ζ*(Ξ - 1)*(ζ + 1)/4
        dN3_dζ = Ξ*ζ*(Ξ - 1)/4 + Ξ*(Ξ - 1)*(ζ + 1)/4
        
        dN4_dΞ = -Ξ*ζ*(ζ - 1)
        dN4_dζ = ζ*(1 - Ξ**2)/2 + (1/2 - Ξ**2/2)*(ζ - 1)
        
        dN5_dΞ = Ξ*(1 - ζ**2)/2 + (1 - ζ**2)*(Ξ/2 + 1/2)
        dN5_dζ = -Ξ*ζ*(Ξ + 1)
        
        dN6_dΞ = -Ξ*ζ*(ζ + 1)
        dN6_dζ = ζ*(1 - Ξ**2)/2 + (1 - Ξ**2)*(ζ/2 + 1/2)
        
        dN7_dΞ = Ξ*(1 - ζ**2)/2 + (1/2 - ζ**2/2)*(Ξ - 1)
        dN7_dζ = -Ξ*ζ*(Ξ - 1)
        
        dN8_dΞ = -2*Ξ*(1 - ζ**2)
        dN8_dζ = -2*ζ*(1 - Ξ**2)
    
    
        dx_dΞ = x0*(dN0_dΞ) + x1*(dN1_dΞ) + x2*(dN2_dΞ) + x3*(dN3_dΞ) + x4*(dN4_dΞ) + x5*(dN5_dΞ) + x6*(dN6_dΞ) + x7*(dN7_dΞ) +x8*(dN8_dΞ)
        dx_dζ = x0*(dN0_dζ) + x1*(dN1_dζ) + x2*(dN2_dζ) + x3*(dN3_dζ) + x4*(dN4_dζ) + x5*(dN5_dζ) + x6*(dN6_dζ) + x7*(dN7_dζ) +x8*(dN8_dζ)
        
        dy_dΞ = y0*(dN0_dΞ) + y1*(dN1_dΞ) + y2*(dN2_dΞ) + y3*(dN3_dΞ) + y4*(dN4_dΞ) + y5*(dN5_dΞ) + y6*(dN6_dΞ) + y7*(dN7_dΞ) +y8*(dN8_dΞ)
        dy_dζ = y0*(dN0_dζ) + y1*(dN1_dζ) + y2*(dN2_dζ) + y3*(dN3_dζ) + y4*(dN4_dζ) + y5*(dN5_dζ) + y6*(dN6_dζ) + y7*(dN7_dζ) +y8*(dN8_dζ)
        
        
        J = array([
        [dx_dΞ, dx_dζ],
        [dy_dΞ, dy_dζ]
        ]).T

        detJ = det(J)
        
        if detJ <= 0.:
            print(f"FATAL! detJ <= 0...")
            break
    
        Jinv = inv(J)

        dN0_dxy = Jinv @ array([ dN0_dΞ, dN0_dζ ])
        dN1_dxy = Jinv @ array([ dN1_dΞ, dN1_dζ ])
        dN2_dxy = Jinv @ array([ dN2_dΞ, dN2_dζ ])
        dN3_dxy = Jinv @ array([ dN3_dΞ, dN3_dζ ])
        dN4_dxy = Jinv @ array([ dN4_dΞ, dN4_dζ ])
        dN5_dxy = Jinv @ array([ dN5_dΞ, dN5_dζ ])
        dN6_dxy = Jinv @ array([ dN6_dΞ, dN6_dζ ])
        dN7_dxy = Jinv @ array([ dN7_dΞ, dN7_dζ ])
        dN8_dxy = Jinv @ array([ dN8_dΞ, dN8_dζ ])

        B = zeros((3,18))
        
        B[0,0]  = dN0_dxy[0]
        B[1,1]  = dN0_dxy[1]
        B[2,0]  = dN0_dxy[1]
        B[2,1]  = dN0_dxy[0]
        
        B[0,2]  = dN1_dxy[0]
        B[1,3]  = dN1_dxy[1]
        B[2,2]  = dN1_dxy[1]
        B[2,3]  = dN1_dxy[0]
        
        B[0,4]  = dN2_dxy[0]
        B[1,5]  = dN2_dxy[1]
        B[2,4]  = dN2_dxy[1]
        B[2,5]  = dN2_dxy[0]
        
        B[0,6]  = dN3_dxy[0]
        B[1,7]  = dN3_dxy[1]
        B[2,6]  = dN3_dxy[1]
        B[2,7]  = dN3_dxy[0]
        
        B[0,8]  = dN4_dxy[0]
        B[1,9]  = dN4_dxy[1]
        B[2,8]  = dN4_dxy[1]
        B[2,9]  = dN4_dxy[0]
        
        B[0,10] = dN5_dxy[0]
        B[1,11] = dN5_dxy[1]
        B[2,10] = dN5_dxy[1]
        B[2,11] = dN5_dxy[0]

        B[0,12] = dN6_dxy[0]
        B[1,13] = dN6_dxy[1]
        B[2,12] = dN6_dxy[1]
        B[2,13] = dN6_dxy[0]

        B[0,14] = dN7_dxy[0]
        B[1,15] = dN7_dxy[1]
        B[2,14] = dN7_dxy[1]
        B[2,15] = dN7_dxy[0]
        
        B[0,16] = dN8_dxy[0]
        B[1,17] = dN8_dxy[1]
        B[2,16] = dN8_dxy[1]
        B[2,17] = dN8_dxy[0]
        
        ε = B @ u_e
        σ = Eσ @ ε
        σBIG.append(σ)
        εBIG.append(ε)
        
    σxx = zeros(9)
    σyy = zeros(9)
    σxy = zeros(9)
    

    σxx[0]= σBIG[0][0]
    σxx[1]= σBIG[1][0]
    σxx[2]= σBIG[2][0]
    σxx[3]= σBIG[3][0]
    σxx[4]= σBIG[4][0]
    σxx[5]= σBIG[5][0]
    σxx[6]= σBIG[6][0]
    σxx[7]= σBIG[7][0]
    σxx[8]= σBIG[8][0]
    
    σyy[0]= σBIG[0][1]
    σyy[1]= σBIG[1][1]
    σyy[2]= σBIG[2][1]
    σyy[3]= σBIG[3][1]
    σyy[4]= σBIG[4][1]
    σyy[5]= σBIG[5][1]
    σyy[6]= σBIG[6][1]
    σyy[7]= σBIG[7][1]
    σyy[8]= σBIG[8][1]
    
    σxy[0]= σBIG[0][2]
    σxy[1]= σBIG[1][2]
    σxy[2]= σBIG[2][2]
    σxy[3]= σBIG[3][2]
    σxy[4]= σBIG[4][2]
    σxy[5]= σBIG[5][2]
    σxy[6]= σBIG[6][2]
    σxy[7]= σBIG[7][2]
    σxy[8]= σBIG[8][2]
    
    
    NN = zeros((9, 9))
    
    gauss_rule_inverse = [
        ( -sqrt(5/3) , -sqrt(5/3) , 1.0, 25/81),
        ( 0          , -sqrt(5/3) , 1.0, 40/81),
        (  sqrt(5/3) , -sqrt(5/3) , 1.0, 25/81),
        ( -sqrt(5/3) , 0          , 1.0, 40/81),
        ( 0          , 0          , 1.0, 64/81),
        (  sqrt(5/3) , 0          , 1.0, 40/81),
        ( -sqrt(5/3) , +sqrt(5/3) , 1.0, 25/81),
        ( 0          , +sqrt(5/3) , 1.0, 40/81),
        (  sqrt(5/3) ,  sqrt(5/3) , 1.0, 25/81),    
    ]
    
    
    
    o = 0
    for Ξ, ζ, wi, wj in gauss_rule_inverse:
        NN[0,o] = (Ξ - 1 )*(ζ - 1) * Ξ * ζ   / 4
        NN[1,o] = (Ξ + 1 )*(ζ - 1) * Ξ * ζ   / 4
        NN[2,o] = (Ξ + 1 )*(ζ + 1) * Ξ * ζ   / 4
        NN[3,o] = (Ξ - 1 )*(ζ + 1) * Ξ * ζ   / 4
        NN[4,o] = (1 - Ξ**2 )*(ζ - 1) * ζ / 2
        NN[5,o] = (1 - ζ**2 )*(Ξ + 1) * Ξ / 2
        NN[6,o] = (1 - Ξ**2 )*(ζ + 1) * ζ / 2
        NN[7,o] = (1 - ζ**2 )*(Ξ - 1) * Ξ / 2
        NN[8,o] = (1 - ζ**2 )*(1 - Ξ**2 )
        o+=1
        
    σXX =  NN @ σxx 
    σYY =  NN @ σyy 
    σXY =  NN @ σxy 

    σBIG = []
    for i in range(9):
        σBIG.append([σXX [i],σYY [i],σXY [i]])
    return εBIG, σBIG








































def quad9_line_load(xy,properties_load, ni, nj):
    
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
    
    
    
    # print (fe)
    
    return fe










