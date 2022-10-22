import numpy as np
from sympy import solve, Symbol
from sympy import interpolating_spline as inter
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import time
import heapq


###Cubo de datos

hdul = fits.open("6302_degraded_withnoise_stray_hinode_mu05.fits")### es  6302 angstrom la linea sintetizada "Hierro"
hdul.info()
image_data = hdul[0].data


###FUNCIONES

### Esta funcion acota el numero de puntos a interpolar
def escoge(perfil,perfil_referencia,n):
    
    diferencia = np.absolute((perfil-np.min(perfil))/(np.max(perfil)-np.min(perfil))-perfil_referencia)#distancia entre los puntos y la linea
    menores=np.array(heapq.nsmallest(6, diferencia))#escoge los valores menores, en este caso los 10 menores 
    indices = list(np.where(np.isin(diferencia, menores, assume_unique=True))[0])#toma el indice de los valores menores
    indices.sort()#organiza de manera creciente los datos
    final_list = list(np.array_split(indices,2))#divide en 2 el arreglo
    puntos_inter1=[]
    puntos_inter2=[]
    for k in range(len(final_list[0])-1):
        #los siguientes 2 condicionales sirven para tomar los valores consecutivos entre puntos        
        if abs(final_list[0][k]-final_list[0][k+1])==1:
            puntos_inter1.append(final_list[0][k])
            puntos_inter1.append(final_list[0][k+1])
        else:
            pass
    for k in range(len(final_list[1])-1):  
        if abs(final_list[1][k]-final_list[1][k+1])==1:
            puntos_inter2.append(final_list[1][k])
            puntos_inter2.append(final_list[1][k+1])
        else:
            pass
    puntos_inter1=list(set(puntos_inter1))  #elimina datos repetidos
    puntos_inter2=list(set(puntos_inter2))
    puntos_inter1.sort()#vuelve a organizar los datos de manera creciente
    puntos_inter2.sort()
    return puntos_inter1, puntos_inter2 #retorna a listas con los valores consecutivos cercanos

def velocidad(i,j, perfil, perfil_referencia, c, lr, n,arreglo_correcto,indice_i,indice_j,sln_1,sln_2):
    
    x=Symbol("x")
    #arg_valor_menor=np.argmin(perfil) ###cambio perifil
    #valor_menor=arreglo_correcto[arg_valor_menor]
    primero, segundo= escoge(perfil, perfil_referencia, n)#usa la funcion anterior para acotar a los dos puntos
    primero_real = [arreglo_correcto[u] for u in primero]
    segundo_real = [arreglo_correcto[u] for u in segundo] 
    primero, segundo=np.array(primero), np.array(segundo)
    primero_real, segundo_real=np.array(primero_real), np.array(segundo_real)
    
    if len(primero)>2 and len(segundo)>2: #condicion por si no hay dos puntos entre la linea de ref
        perfil_editado_primero = image_data[j,i,0,primero]#primeros puntos que se van a interpolar
        perfil_editado_segundo = image_data[j,i,0,segundo]#los otros puntos que se van a interpolar
        
        
        primera_interpolacion = inter(2,x,primero_real,\
                                      ((perfil_editado_primero-np.min(perfil_editado_primero))/
                                       (np.max(perfil_editado_primero)-np.min(perfil_editado_primero))))#interpolacion
        segunda_interpolacion = inter(2,x,segundo_real,\
                                      ((perfil_editado_segundo-np.min(perfil_editado_segundo))/
                                       (np.max(perfil_editado_segundo)-np.min(perfil_editado_segundo))))

        a=solve([primera_interpolacion-perfil_referencia],[x])#sln de la interpolacion con la linea
        b=solve([segunda_interpolacion-perfil_referencia],[x])  
        
        if b!=[] and a==[]:
            if len(b[0])>=2:
                l=(b[0][0]+b[0][1])/2
                vel = float(((l-lr)/lr)*c) 
                return lr  
            else:
                pass
            
        if a!=[] and b==[]:    
            if len(a[0])>=2:
                l=(a[0][0]+a[0][1])/2
                vel = float(((l-lr)/lr)*c)
                return vel
            else: 
                pass
            
        if a!=[] and b!=[]:
            l=(a[0][0]+b[0][0])/2
            vel = float(((l-lr)/lr)*c) 
            return vel
            
        else:
            indice_i.append(i)
            indice_j.append(j)
            if len(a)!=0:
                sln_1.append(a)  
            else:
                sln_1.append(None)   
            if len(b)!=0:
                sln_2.append(b)
            else:
                sln_2.append(None)
                pass       
            return 0
        
    else:
        indice_i.append(i)
        indice_j.append(j)
        sln_1.append(None)   
        sln_2.append(None)
        return 0
        pass 
    
if __name__ == "__main__":

       print("""Código hecho por Óscar Calvo, establece las velocidad
en la línea de la visual por medio del método de bisectores
con bisector en I=0.7""")



