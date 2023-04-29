# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 21:46:34 2023

@author: Zaneta
"""

import math
from math import *
import numpy as np
import argparse

###### w tym przykładzie to było zdefiniowane jako objekt wiec dodałam ale nw czy git
o = object()


class Transformacje():
    
    def __init__(self, elipsoida: str ='GRS80'):
        """
        Parametry elipsoidy:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            e2 - mimośród^2
            
        """
        self.parametry_elipsoidy = {
            "GRS80": (6378137.0, 6356752.31424518),
            "WGS84": (6378137.0, 6356752.31414036),
            "Krasowskiego": (6378245.0, 6356863.019),
        }
        try:
            self.a, self.b = self.parametry_elipsoidy[elipsoida]
        except KeyError:
            raise NotImplementedError(f"{elipsoida} elipsoida nie została zaimplementowana")
        
        self.flat = (self.a - self.b) / self.a
        self.e = np.sqrt(2 * self.flat - self.flat ** 2)
        self.e2 = (2 * self.flat - self.flat ** 2)
            
        
    
    def xyz2flh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena 
        Funkcja przelicza współrzędne ortokartezjańskie (XYZ) na współrzędne geodezyjne (długość,
        szerokość i wysokość elipsoidalna) - (phi, lambda, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia współrzędnej phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        
        Parametry:
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie ortokartezjańskim, 

        Zwraca:
        -------
        f (phi)
            [stopnie dziesiętne] - szerokość geodezyjna
        l (lambda)
            [stopnie dziesiętne] - długośc geodezyjna
        h : TYPE
            [metry] - wysokość elipsoidalna
            
        output [STR] - opcjonalne, domyślne 
            dec_degree - stopnie dziesiętne
            dms - stopnie, minuty sekundy
        
        """
        p = np.sqrt(X**2 + Y**2)        #promień równoleżnika 
        f = np.arctan(Z/(p*(1 - self.e2)))     #przybliżona wartosc fi
        while True:
            N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
            h = (p/np.cos(f)) - N
            fpop = f
            f = np.arctan(Z / (p * (1 - (self.e2 * N) / (N + h))))
            if abs(fpop - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)
      
  #### o tu cos z tym trzeba zrobic chyba, zeby przeliczało na stopnie, minuty, sekundy 
 # tak analogicznie do tego co dał w przykładzie - dodałam funkcje deg2dms i  te self w przeliczaniu stopni na stopnie, minuty i sekundy jak cos
        def dms(x):
            sig = ''
            if x < 0:
                sig = '-'
                x = abs(x)  
            x = x * 180/pi
            d = int(x)
            m = int((x-d) * 60)
            s = (x - d - m/60) * 3600
            return(sig,"%3d %2d %7.5f" %(d,m,s))
        def deg2dms(degrees):
            # przelicza stopnie na stopnie, minuty i sekundy
            d = int(degrees)
            m = int((degrees - d) * 60)
            s = (degrees - d - m/60) * 3600
            # Zwracamy wynik jako krotkę
            return (d, m, s)
        if output == "dec_degree":
            return degrees(f), degrees(l), h 
        elif output == "dms":
            f = self.deg2dms(degrees(f))
            l = self.deg2dms(degrees(l))
            return f"{f[0]:02d}:{f[1]:02d}:{f[2]:.2f}", f"{l[0]:02d}:{l[1]:02d}:{l[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - jednostka wyjsciowa nie została zdefiniowana")
              
                    
    
    def flh2xyz(self, f, l, h, output = 'metry'):
        """
        Funkcja odwrotna do algorytmu Hirvonena. Przelicza współrzędne geodezyjne (phi, lambda, h) na 
        współrzędne ortokartezjańskie (XYZ).
        Dane:
        ----------
        f, l, h : FLOAT
             współrzędne w układzie geodezyjnym

        Zwraca:
        -------
        X : TYPE
            [metry]
      
        Y : TYPE
            [metry]
            
        Z: TYPE
            [metry]
            
        """
        f = radians(f)
        l = radians(l)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)   #promień krzywizny w I wertykale
        Xk = (N + h) * np.cos(f) * np.cos(l)
        Yk = (N + h) * np.cos(f )* np.sin(l)
        Zk = (N * (1 - self.e2) +h) * np.sin(f)
        return(Xk, Yk, Zk)
 
        
# XYZ2neu to nie wiem czy git jest   
    
    def xyz2neu(self, X, Y, Z, f, l):
       """
       Funkcja odwrotna do algorytmu Hirvonena. Przelicza współrzędne geodezyjne (phi, lambda, h) na 
       współrzędne ortokartezjańskie (XYZ).
       Dane:
       ----------
       f, l, h : FLOAT
            współrzędne w układzie geodezyjnym

       Zwraca:
       -------
       X : TYPE
           [metry]
     
       Y : TYPE
           [metry]
           
       Z: TYPE
           [metry]
           
       """
       
       f = radians(f)
       l = radians(l)
       
       dX = np.array([X, Y, Z])
       R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                     [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
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
        R = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        r = sqrt(X**2 + Y**2 + Z**2)
        f_g = asin(Z / r)    #geocentryczne fi i lambda punktu XYZ
        l_g = atan2(Y, X)
        N = -sin(f_g) * cos(l_g) * X - sin(f_g) * sin(l_g) * Y + cos(f_g) * Z
        E = -sin(l) * X + cos(l) * Y
        U = cos(f_g) * cos(l_g) * X + cos(f_g) * sin(l_g) * Y + sin(f_g) * Z - R
        return (N, E, U)
    
    

    def fl_to_uk2000(self, f, l):  
        """
        Funkcja przelicza współrzędne geodezyjne elipsoidalne (fi,lambda) na współrzędne płaskie (X,Y) w układzie PL-2000.
        
        Układ Pl_2000 powstał w wyniku zastosowania odwzorowania Gaussa-Krügera dla elipsoidy GRS'80 w czterech 
        trzystopniowych strefach o południkach osiowych 15°E, 18°E, 21°E i 24°E. 
        Strefy oznaczone są kolejno numerami: 5,6,7,8.                                             
        Skala długości odwzorowania na południkach osiowych wynosi m0 = 0,999923.
        
        Parametry:
        ----------
        f, l : FLOAT
             współrzędne geodezyjne elipsoidalne

        Zwraca:
        -------
        X_2000 : TYPE
                    [metry]
      
        Y_2000 : TYPE
                   [metry]
                   
        """
        f = radians(f)
        l = radians(l)
        
        L0 = 0
        n = 0
        if l > radians(13.5) and l < radians(16.5):
            L0 = L0 + radians(15)
            n = n + 5
        if l > radians(16.5) and l < radians(19.5): 
            L0 = L0 + radians(18)
            n = n + 6
        if l > radians(19.5) and l < radians(22.5): 
            L0 = L0 + radians(21)
            n = n + 7
        if l > radians(22.5) and l < radians(25.5): 
            L0 = L0 + radians(24)
            n = n + 8           
        
        b2 = (self.a**2) * (1 - self.e2)
        e22 = ((self.a**2) - b2) / b2
        dl = l - L0
        t = tan(f)
        eta2 = e22 * ((cos(f))**2)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        
        A0 = 1 - (self.e2/4) - (3 * self.e2**2/64) - (5 * self.e2**3/256)
        A2 = 3/8 *((self.e2) + (self.e2**2/4) + (15 * self.e2**3/128))
        A4 = 15/256 * ((self.e2**2) + (3 * self.e2**3/4))
        A6 = 35 * self.e2**3/3072
        sig = self.a * (A0 * f - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        
        Xgk = sig + ((dl**2)/2)*N*(sin(f))*(cos(f))*(1+((dl**2)/12)*((cos(f))**2)*(5-(t**2)+9*(eta2)+4*((eta2)**2))+((dl**4)/360)*((cos(f))**4)*(61-(58*(t**2))+(t**4)+(270*(eta2))-(330*(eta2)*(t**2))))                         
        Ygk = dl*N*(cos(f))*(1+((dl**2)/6)*((cos(f))**2)*(1-(t**2)+(eta2))+((dl**4)/120)*((cos(f))**4)*(5-(18*(t**2))+(t**4)+(14*eta2)-(58*(eta2)*(t**2))))
        
        X_2000 = Xgk * 0.999923
        Y_2000 = Ygk * 0.999923 + (n*1000000) + 500000
        return(X_2000, Y_2000)

    
    
    def fl_to_uk1992(self, f, l):  
        """
        Funkcja przelicza współrzędne geodezyjne elipsoidalne (fi,lambda) na współrzędne płaskie (X,Y) w układzie PL-1992.
        
        Układ Pl-1992 jest oparty na odwzorowaniu Gaussa-Krügera dla elipsoidy GRS'80 w jednej dziesięciostopniowej 
        strefie. Początkiem układu jest punkt przecięcia południka zerowego 19°E z obrazem równika. Południk zerowy
        odwzorowuje się na linię prostą w skali m0 = 0,9993.
             
        Parametry:
        ----------
        f, l : FLOAT
             współrzędne geodezyjne elipsoidalne

        Zwraca:
        -------
        X_1992 : TYPE
                    [metry]
      
        Y_1992 : TYPE
                    [metry]
            
        """
        L0_92 = 19*pi/180
        f = radians(f)
        l = radians(l)
        
        b2 = (self.a**2) * (1 - self.e2)
        e22 = ((self.a**2) - b2) / b2
        dl = l - L0_92
        t = tan(f)
        eta2 = e22 * ((cos(f))**2)
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        
        A0 = 1 - (self.e2/4) - (3 * self.e2**2/64) - (5 * self.e2**3/256)
        A2 = 3/8 *((self.e2) + (self.e2**2/4) + (15 * self.e2**3/128))
        A4 = 15/256 *((self.e2**2) + (3 * self.e2**3/4))
        A6 = 35 * self.e2**3/3072
        sig = self.a * (A0 * f - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        
        Xgk = sig + ((dl**2)/2)*N*(sin(f))*(cos(f))*(1+((dl**2)/12)*((cos(f))**2)*(5-(t**2)+9*(eta2)+4*((eta2)**2))+((dl**4)/360)*((cos(f))**4)*(61-(58*(t**2))+(t**4)+(270*(eta2))-(330*(eta2)*(t**2))))                         
        Ygk = dl*N*(cos(f))*(1+((dl**2)/6)*((cos(f))**2)*(1-(t**2)+(eta2))+((dl**4)/120)*((cos(f))**4)*(5-(18*(t**2))+(t**4)+(14*eta2)-(58*(eta2)*(t**2))))
        
        X_1992 = (Xgk * 0.9993) - 5300000
        Y_1992 = (Ygk * 0.9993) + 500000
        return(X_1992, Y_1992)




if __name__ == "__main__":
    #1
    tr1 = Transformacje(elipsoida = "WGS84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    f, l, h = tr1.xyz2flh(X, Y, Z)
    print('f = '"%7.5f" % f, 'l = ' "%7.5f" % l, 'h = '"%7.5f" % h)
     
    #2
    tr2 = Transformacje(elipsoida = "WGS84")
    # dane flh
    f = 52.0972722; l = 21.0315333279; h = 141.3986623911
    X, Y, Z = tr2.flh2xyz(f,l,h)
    print("%11.5f" % X, "%11.5f" % Y, "%11.5f" % Z)
    
    #3
    tr3 = Transformacje(elipsoida = "WGS84")
    # dane flh
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170; f = 52.0972722; l = 21.0315333279
    n, e, u = tr3. xyz2neu(X, Y, Z, f,l)
    print(n,e,u)
    
    #4
    tr4 = Transformacje(elipsoida = "WGS84")
    # dane flh
    f = 52.0972722; l = 21.0315333279
    X_2000,Y_2000 = tr4.fl_to_uk2000(f,l)
    print("%11.5f" % X_2000, "%11.5f" % Y_2000)
    
    #5
    tr5 = Transformacje(elipsoida = "WGS84")
    # dane flh
    f = 52.0972722; l = 21.0315333279
    X_1992,Y_1992 = tr5.fl_to_uk1992(f,l)
    print("%11.5f" % X_1992, "%11.5f" % Y_1992)
  
    
"""  
 #===============================================
# próbuje tu tą biblioteke arparse czy cos   
    
    parser = argparse.ArgumentParser()
    parser.add_argument('x', type = float, help ='współrzędna x punktu')
    parser.add_argument('y', type = float, help ='współrzędna y punktu')
    parser.add_argument('z', type = float, help ='współrzędna z punktu')
    args = parser.parse_args()
    xyz2flh(args.x, args.y, args.z)  
    print(f'Uzyskane współrzędne punktu to ({args.x}, {args.y}, {args.z})')

####========================================================

# funkcja na transformowanie danych z pliku(on tam dał w przykładzie taki plik - wsp_inp, wiec tego chyba uzyje)
# Nie wiem czy ona działa bo sie nic nie dzieje(nie ma nawet błędu) jak ją wywoluje i nw dlaczego...
    def przelicz_dane_plik(input_file, output_file):
        
#       Funkcja ta przyjmuje dwa argumenty: input_file to nazwa pliku wejściowego, a output_file to nazwa pliku wynikowego. 
        Funkcja korzysta z obiektu Transformacje utworzonego wcześniej i przelicza współrzędne geocentryczne dla każdej linii z pliku wejściowego.
        
        t = Transformacje()
        with open(input_file, 'r') as f:
            lines = f.readlines()
        with open(output_file, 'w') as f:
            for line in lines:
                if line.startswith('#'):
                    continue
                x, y, z = [float(x.strip()) for x in line.split(',')]
                f.write(f"{t.xyz2flh(x, y, z)}\n")

# PRZYKLAD UZYCIA: Po uruchomieniu funkcji przelicz_dane_plik("wsp_inp.txt", "wyniki.txt") zostanie stworzony plik 'wyniki.txt' z wynikami obliczeń dla danych z pliku 'wsp_inp.txt'.

WYNIKI = przelicz_dane_plik("wsp_inp.txt", "wyniki.txt")
        
"""
def transformacje_plik(X, Y, Z):
    trans = Transformacje()
    with open('wyniki.txt', 'w') as wyniki:
        wyniki.write('Wyniki transformacji współrzędnych ECEF na geodezyjne\n\n')
        wyniki.write('Nazwa elipsoidy: GRS80\n')
        wyniki.write('Metoda transformacji: Algorytm Hirvonena\n\n')
        wyniki.write('Współrzędne ECEF [m]:\nX[m]         Y[m]        Z[m]\n')
        for x, y, z in zip(X, Y, Z):
            f, l, h = trans.xyz2flh(x, y, z, output='dec_degree')
            wyniki.write(f'{x:.3f}, {y:.3f}, {z:.3f} -> {f:.10f}, {l:.10f}, {h:.3f}\n')

with open('wsp_in.txt', 'r') as f:
    next(f)  # pomijamy nagłówek
    X, Y, Z = [], [], []
    for line in f:
        x, y, z = map(float, line.strip().split(','))
        X.append(x)
        Y.append(y)
        Z.append(z)
    transformacje_plik(X, Y, Z)
