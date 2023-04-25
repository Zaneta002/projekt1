# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 21:46:34 2023

@author: Zaneta
"""

import math
from math import *
import numpy as np


class Transformacje():
    
    def __init__(self, elipsoida: str ='GRS80'):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            e2 - mimośród^2
            
        """
        self.parametry_elipsoid = {
            "GRS80": (6378137.0, 6356752.31424518),
            "WGS84": (6378137.0, 6356752.31414036),
            "Krasowskiego": (6378245.0, 6356863.019),
        }
        try:
            self.a, self.b = self.parametry_elipsoid[elipsoida]
        except KeyError:
            raise NotImplementedError(f"{elipsoida} elipsoida nie została zaimplementowana")
        
        self.flat = (self.a - self.b) / self.a
        self.e = np.sqrt(2 * self.flat - self.flat ** 2)
        self.e2 = (2 * self.flat - self.flat ** 2)
            

    
    
    def hirvonen(self, X, Y, Z, a, e2):
        """
        DODAC DOKUMENTACJE 
        
        """
    
        p = np.sqrt(X**2 + Y**2)        #promień równoleżnika 
        f = np.arctan(Z/(p*(1 - e2)))     #przybliżona wartosc fi
        while True:
            N = a / np.sqrt(1-e2 * np.sin(f)**2)
            h = (p/np.cos(f)) - N
            fpop = f
            f = np.arctan(Z/(p*(1-(e2*N)/(N+h))))
            if abs(fpop - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)
        return(f,l,h) 
    
    
    def flh2XYZ(self, f, l, h, a, e2):
        """
        DODAC DOKUMENTACJE
        
        """
        N = a / np.sqrt(1-e2 * np.sin(f)**2)   #promień krzywizny w I wertykale
        Xk = (N + h)*np.cos(f)*np.cos(l)
        Yk = (N + h)*np.cos(f)*np.sin(l)
        Zk = (N*(1-e2)+h)*np.sin(f)
        return(Xk, Yk, Zk)
    
   
    
        
# XYZ2neu to nie wiem czy git jest   
    
    def XYZ2neu(self, X, Y, Z, f, l):
        """
        DOKUMEN
        
        """

        dX = np.array([X, Y, Z])
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                     [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                     [np.cos(f), 0, np.sin(f)]])
        dneu = R.T @ dX
        s = np.sqrt(dneu @ dneu)
        alfa = np.arctan2(dneu[1], dneu[0])
        z = np.arccos(dneu[2]/s)
        return (s, alfa, z)
    
# inna wersja luknij Madzia
   
    def xyz_to_neu(self, X, Y, Z, f, l):
        f = radians(f)
        l = radians(l)
        R = a / np.sqrt(1-e2 * np.sin(f)**2)
        r = sqrt(X**2 + Y**2 + Z**2)
        f_g = asin(Z / r)    #geocentryczne fi i lambda punktu XYZ
        l_g = atan2(Y, X)
        N = -sin(f_g) * cos(l_g) * X - sin(f_g) * sin(l_g) * Y + cos(f_g) * Z
        E = -sin(l) * X + cos(l) * Y
        U = cos(f_g) * cos(l_g) * X + cos(f_g) * sin(l_g) * Y + sin(f_g) * Z - R
        return (N, E, U)
    
 
    
 

    def fl_to_uk2000(self, f, l, L0, n, a, e2):  #L0 - południk0 (15,18,21,24); n - strefa(5,6,7,8); nie wiem jak to ustawić
        """
        DOKUM
        
        """
        b2 = (a**2)*(1-e2)
        e22 = ((a**2) - b2) / b2
        dl = l - L0
        t = tan(f)
        eta2 = e22 * ((cos(f))**2)
        N = a / np.sqrt(1-e2 * np.sin(f)**2)
        
        A0 = 1 - (e2/4) - (3 * e2**2/64) - (5 * e2**3/256)
        A2 = 3/8 *((e2) + ( e2**2/4) + (15 * e2**3/128))
        A4 = 15/256 *((e2**2) + (3 * e2**3/4))
        A6 = 35 * e2**3/3072
        sig = a * (A0 * f - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        
        Xgk = sig + ((dl**2)/2)*N*(sin(f))*(cos(f))*(1+((dl**2)/12)*((cos(f))**2)*(5-(t**2)+9*(eta2)+4*((eta2)**2))+((dl**4)/360)*((cos(f))**4)*(61-(58*(t**2))+(t**4)+(270*(eta2))-(330*(eta2)*(t**2))))                         
        Ygk = dl*N*(cos(f))*(1+((dl**2)/6)*((cos(f))**2)*(1-(t**2)+(eta2))+((dl**4)/120)*((cos(f))**4)*(5-(18*(t**2))+(t**4)+(14*eta2)-(58*(eta2)*(t**2))))
        
        X_2000 = Xgk * 0.999923
        Y_2000 = Ygk * 0.999923 + (n*1000000) + 500000
        return(X_2000, Y_2000)
    
    
    def fl_to_uk1992(self, f, l, L0, a, e2):  #nie wiem jak ustawić zeby L0 było zawsze 19*pi/180
        """
        DOKUM
        
        """
        b2 = (a**2)*(1-e2)
        e22 = ((a**2) - b2) / b2
        dl = l - L0
        t = tan(f)
        eta2 = e22 * ((cos(f))**2)
        N = a / np.sqrt(1-e2 * np.sin(f)**2)
        
        A0 = 1 - (e2/4) - (3 * e2**2/64) - (5 * e2**3/256)
        A2 = 3/8 *((e2) + ( e2**2/4) + (15 * e2**3/128))
        A4 = 15/256 *((e2**2) + (3 * e2**3/4))
        A6 = 35 * e2**3/3072
        sig = a * (A0 * f - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        
        Xgk = sig + ((dl**2)/2)*N*(sin(f))*(cos(f))*(1+((dl**2)/12)*((cos(f))**2)*(5-(t**2)+9*(eta2)+4*((eta2)**2))+((dl**4)/360)*((cos(f))**4)*(61-(58*(t**2))+(t**4)+(270*(eta2))-(330*(eta2)*(t**2))))                         
        Ygk = dl*N*(cos(f))*(1+((dl**2)/6)*((cos(f))**2)*(1-(t**2)+(eta2))+((dl**4)/120)*((cos(f))**4)*(5-(18*(t**2))+(t**4)+(14*eta2)-(58*(eta2)*(t**2))))
        
        X_1992 = (Xgk * 0.9993) - 5300000
        Y_1992 = (Ygk * 0.9993) + 500000
        return(X_1992, Y_1992)


