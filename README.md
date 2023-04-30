Program został napisany w języku programowania Python. Służy do transformacji współrzędnych między poszczególnymi układami.
Umożliwia:
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne geodezyjne (fi,lambda,h), 
* zamianę odwrotną do powyższej (fi,lambda,h) -> (X,Y,Z),
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne topocentryczne NEU (northing,easting,up),
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 2000,
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 1992.

Program obsługuje elipsoidy GRS80, WGS84 oraz elipsoidę Krasowskiego. 
Aby program poprawnie działał należy mieć na danym komputerze zainstalowaną aplikację Python w wersji 3.10.10, 
a także zainstalowane biblioteki - NumPy, Math oraz Argparse.
Skrypt został napisany przy użyciu środowiska programistycznego Spyder w systemie operacyjnym Microsoft Windows 10.

Program można uruchomić przy pomocy Wiersza poleceń, co jest możliwe dzięki zastosowaniu biblioteki Argparse.
Wówczas przekazujemy jako argumenty wartości współrzędnych do przeliczenia oraz parametry elipsoidy. 
Opcjonalnie możemy określić jednostkę wyjściową przez podanie wartości opcji '--output'.
Przykładowe polecenie, które uruchomi program wygląda następująco:
"python skrypt.py --elipsoida WGS84 --X 54623.34524 --Y 64697.89021 --Z 12523.45236 --output dec_degree"





