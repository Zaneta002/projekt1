# Funkcjonalności programu
Program służy do transformacji współrzędnych między poszczególnymi układami.
Umożliwia:
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne geodezyjne (fi,lambda,h), 
* zamianę odwrotną do powyższej (fi,lambda,h) -> (X,Y,Z),
* zamianę współrzędnych ortokartezjańskich (X,Y,Z) na współrzędne topocentryczne NEU (northing,easting,up),
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 2000,
* zamianę współrzędnych geodezyjnych (fi,lambda) na współrzędne w układzie 1992.
  
Program obsługuje elipsoidy GRS80, WGS84 oraz elipsoidę Krasowskiego.

# System operacyjny i programy
Program został napisany dla systemu operacyjnego Microsoft Windows 10 w języku programowania Python w środowisku programistycznym Spyder.

# Wymagania techniczne 
Aby program poprawnie działał należy mieć na danym komputerze zainstalowaną aplikację Python w wersji 3.10.10, 
a także zainstalowane biblioteki - NumPy, Math oraz Argparse.

# Działanie programu
Program można uruchomić przy pomocy Wiersza poleceń, co jest możliwe dzięki zastosowaniu biblioteki Argparse.
Wówczas przekazujemy jako argumenty wartości współrzędnych do przeliczenia oraz parametry elipsoidy. 
Opcjonalnie możemy określić jednostkę wyjściową przez podanie wartości opcji '--output'.
Przykładowe polecenie, które uruchomi program wygląda następująco:
"python skrypt.py --elipsoida WGS84 --X 54623.34524 --Y 64697.89021 --Z 12523.45236 --output dec_degree".

Możliwe jest również użycie skryptu dla transformacji wielu współrzędnych jednocześnie. Wówczas argumentem wejściowym może być plik z zestawem danych do transformacji, a wynikiem działania programu jest plik z przetransformowanymi współrzędnymi. 
Przykład takiego użycia przedstawiony jest poniżej:

Plik wejściowy:

Współrzedne geocentryczny ECEF stacji pemanentnej GNSS
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu

  X[m]         Y[m]        Z[m]
# -----------------------------------------------------
3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
3664940.520,1409153.570,5009571.167
3664940.530,1409153.560,5009571.168
3664940.520,1409153.590,5009571.170
3664940.514,1409153.584,5009571.166
3664940.525,1409153.575,5009571.166
3664940.533,1409153.564,5009571.169
3664940.515,1409153.590,5009571.170
3664940.514,1409153.584,5009571.169
3664940.515,1409153.595,5009571.169
3664940.513,1409153.584,5009571.171

Plik wyjściowyy:
Wyniki transformacji współrzędnych ECEF na geodezyjne

Nazwa elipsoidy: GRS80
Metoda transformacji: Algorytm Hirvonena

Współrzędne ECEF [m]:

X[m]         Y[m]        Z[m]

3664940.500, 1409153.590, 5009571.170 -> 52.0972722184, 21.0315333328, 141.399
3664940.510, 1409153.580, 5009571.167 -> 52.0972721611, 21.0315331442, 141.400
3664940.520, 1409153.570, 5009571.167 -> 52.0972721204, 21.0315329557, 141.403
3664940.530, 1409153.560, 5009571.168 -> 52.0972720852, 21.0315327671, 141.408
3664940.520, 1409153.590, 5009571.170 -> 52.0972720860, 21.0315332281, 141.410
3664940.514, 1409153.584, 5009571.166 -> 52.0972721189, 21.0315331778, 141.402
3664940.525, 1409153.575, 5009571.166 -> 52.0972720690, 21.0315329976, 141.406
3664940.533, 1409153.564, 5009571.169 -> 52.0972720606, 21.0315328059, 141.411
3664940.515, 1409153.590, 5009571.170 -> 52.0972721191, 21.0315332542, 141.407
3664940.514, 1409153.584, 5009571.169 -> 52.0972721355, 21.0315331778, 141.405
3664940.515, 1409153.595, 5009571.169 -> 52.0972721009, 21.0315333223, 141.408
3664940.513, 1409153.584, 5009571.171 -> 52.0972721532, 21.0315331830, 141.406

# Znane błędy, które nie zostały jeszcze naprawione
Transformacje:
* Krasowski (fi,lambda) -> (XY) 2000
* Krasowski (fi,lambda) -> (XY) 1992
  
dają błedne rezultaty i nie powinny być używane.




