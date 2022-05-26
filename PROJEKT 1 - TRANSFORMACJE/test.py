# -*- coding: utf-8 -*-
"""
Created on Mon May 16 21:14:50 2022

@author: moja piękna
"""

#testowanie aplikacji 

#połączenie z testowaną aplikacją
from transformacje import Transformacje

#testowanie funkcji 
if __name__ == "__main__":
    test = Transformacje(model = "wgs84")
    x = 3664940.5
    y = 1409153.59
    z = 5009571.17
    #testowanie zamiana x,y,z na fi,lam,h
    p,l,h = test.xyz2plh(x, y, z)
    print('fi, lam, h =',round(p,5),round(l,5),round(h,3))
    #testowanie zamiana fi,lam,h na x,y,z
    X,Y,Z = test.plh2xyz(p, l, h)
    print('X,Y,Z =',round(X,2),round(Y,2),round(Z,2))
    #wyznaczenie współrzędnych w układzie 2000
    X2000, Y2000 = test.uklad2000(p,l,h)
    print('x2000, y2000 =', round(X2000,2),round(Y2000,2))
    #wyznaczenie współrzędnych w układzie 1992
    X92, Y92 = test.uklad1992(p,l,h)
    print('x1992, y1992 =',round(X92,2),round(Y92,2))
    
    
    
    
    
    
    
