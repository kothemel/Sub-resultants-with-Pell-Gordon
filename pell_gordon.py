"""
Created on Thu Mar  6 11:40:43 2014

University of Thessaly
@author: kothemel
mail: kothemel@uth.gr
"""
# Subresultants PRS using Pell Gordon Theorem
# Until now = complete prs + pivot
# No usage of det of any matrix

import sympy  as sp
import giacpy as gp
i=sp.var('i')
x=sp.var('x')

def my_subresultants (p, q, x):
    """ So far we obtain the Sturm sequence, whether complete or incomplete.
        
        Next, we will obtain the subresultant prs 
    """    
    subresL = [p, q]

    my_poly = p
    p = my_poly
    if(sp.LC(my_poly) < 0):
        p = -p
        q = -q
    degree_rem_poly = 0
    poly1 = p
    
    # polynomial q is the first derivative of p
    # q = sp.simplify(sp.diff(p,x))
    poly2 = q
    ui = 0;
    vi = 0;
    rho_neg = sp.LC(poly1)  # is the LC of poly1 -> we'll use it in the rest of the cases...
    r0_i = rho_neg          # this will always be multiplied in the denominator
                            # of the formula         
    mass_of_u = pi = p0 = 1
    fd = sp.degree(poly1)
    while(1):
        rhoi  = sp.LC(poly2) # is the rho0, rho1 of the formula
        rem_poly = -sp.rem(poly1,poly2,x)
        p_plusplus = sp.degree(poly2) - sp.degree(rem_poly)
        pold = p_plusplus
        degree_rem_poly = sp.degree(rem_poly) # so as to know where to stop
        rem_LC = sp.LC(rem_poly)              # the LC of the new remnant
        poly1=poly2                           # refresh poly1 and poly2
        poly2=rem_poly
                                             # if we are in the first loop this is p1
                                                # if we are in the second loop this is p2
        fd = sp.LC(poly1)
        
        ui = sp.summation(i,(i,1,pi));
        vi = vi + pi;
        mass_of_u = (-1)**ui * mass_of_u;
        sign = mass_of_u * (-1)**vi;
        pi = p_plusplus;
        
        if(p_plusplus>1):
            r0_i = r0_i * rhoi**(1+p0) # the result of the Pell_Gordon formula
        else:
            r0_i = r0_i * rhoi**(p_plusplus+p0)
                            # multiply the denominator od the formula with
                            # the r0_i so as to find the det of the current submatrix

        LC = (rem_poly * r0_i)/sp.LC(my_poly)/sign
        # print(LC)
        subresL.append(LC)        
        
        # when the remnant is finally a number
        # LC and degree function throw an exception
        # so let's take a different case
        if(degree_rem_poly==1):            
            rhoi  = sp.LC(poly2)
            rem_poly = -sp.rem(poly1,poly2,x)
            p_plusplus = sp.degree(poly1)-1
            rem_LC = rem_poly
            
            ui = sp.summation(i,(i,1,pi))
            vi = vi + pi
            mass_of_u = (-1)**ui * mass_of_u
            sign = mass_of_u * (-1)**vi
            
            if(p_plusplus>1):
                
                if(rhoi<0):
                    r0_i = -r0_i*fd**(pold-p0) *(rhoi**(pold+p0))/sign # the result of the Pell_Gordon formula
                else:
                    r0_i =  r0_i*fd**(pold-p0) *(rhoi**(pold+p0))
            else:
                r0_i = r0_i * rhoi**(1+p0)
                
            LC = (rem_LC * r0_i)/sign
            #print(LC)
            subresL.append(LC)
            break
    #print("End of Current Computation")    
    return subresL
    
print("\n")
print('Example-1')
p = x**6+x**5-x**4-x**3+x**2-x+1
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-2')
p = x**3-7*x+7
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-3')
p = x**5-3*x-1
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-4')
p = -2*x**5+7*x**3+9*x**2-3*x+1
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq((p))[1])
print("\n")
print('Example-5')
p = 7*x**5+3*x**4-8*x-1
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-6')
p = 2*x**5-4*x**4+7*x+7
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-7')
p = 4*x**5-3*x**4+7
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-8')
p = 9*x**5+5*x**3+9
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-9')
p = x**8+x**6-3*x**4-3*x**3+8*x**2+2*x-5
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")
print('Example-10')
p = 3*x**6+5*x**4-4*x**2-9*x+21
q = sp.diff(p, x, 1)
print(my_subresultants(p, q, x))
print(gp.sturmseq(p)[1])
print("\n")


print("End of Program")
