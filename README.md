Stworzony program służy do transformacji współrzędnych między poszczególnymi układami.
Umożliwia:
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne geodezyjne (fi,lambda,h), 
* zamianę odwrotną do powyższej (fi,lambda,h) -> (X,Y,Z),
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne topocentryczne NEU (northing,easting,up),
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 2000,
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 1992.

Program obsługuje elipsoidy GRS80, WGS84 oraz elipsoidę Krasowskiego. 
Aby program poprawnie działał należy mieć na danym komputerze zainstalowaną aplikację Python w wersji Python 3.10.10 oraz zainstalowane biblioteki - numpy, math oraz moduł argparse.
Program został napisany dla systemu operacyjnego Microsoft Windows.