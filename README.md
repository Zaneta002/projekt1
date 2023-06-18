## Funkcjonalności programu
Program służy do transformacji współrzędnych między poszczególnymi układami.
Umożliwia:
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne geodezyjne (fi,lambda,h), 
* zamianę odwrotną do powyższej (fi,lambda,h) -> (X,Y,Z),
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne topocentryczne NEU (northing,easting,up),
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 2000,
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 1992.
  
Program obsługuje elipsoidy GRS80, WGS84 oraz elipsoidę Krasowskiego.

## System operacyjny i programy
Program został napisany dla systemu operacyjnego Microsoft Windows 10 w języku programowania Python w środowisku programistycznym Spyder.

## Wymagania techniczne 
Aby program poprawnie działał należy mieć na danym komputerze zainstalowaną aplikację Python w wersji 3.10.10, 
a także zainstalowane biblioteki - NumPy, Math oraz Argparse.

## Działanie programu
Program można uruchomić przy pomocy Wiersza poleceń, co jest możliwe dzięki zastosowaniu biblioteki Argparse.
Wówczas przekazujemy jako argumenty wartości współrzędnych do przeliczenia oraz parametry elipsoidy. 
Opcjonalnie możemy określić jednostkę wyjściową przez podanie wartości opcji '--output'.
Przykładowe polecenie, które uruchomi program wygląda następująco:
"python skrypt.py --elipsoida WGS84 --X 54623.34524 --Y 64697.89021 --Z 12523.45236 --output dec_degree".

## Sposób przekazania danych
W celu zaimportowania danych do programu, należy wpisać poszczególne wartości współrzędnych w wierszu poleceń z odpowiednim wywołaniem skryptu.
Współrzędne należy wpisywać w metrach (w przypadku transformacji XYZ -> BLH), stopniach dziesiętnych (w przypadku transformacji: BLH -> XYZ, 
BL -> 2000, BL -> 1992) oraz w metrach (XYZ) i stopniach dziesiętnych (BL) w przypadku transformacji XYZ -> neu.

## Znane błędy, które nie zostały jeszcze naprawione
Transformacje:
* BL (Krasowski) -> XY 2000
* BL (Krasowski) -> XY 1992
  
dają błedne rezultaty i nie powinny być używane.




