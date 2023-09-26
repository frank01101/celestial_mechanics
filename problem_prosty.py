#!/usr/bin/env python
#-*- coding: utf-8 -*-
# mechanika nieba I: problem prosty
# autor: Franciszek Humieja
# wersja 1.1 (26.I.2017)

from sys import argv
import opracuj_dane as opr
import math as m
from scipy.optimize import fsolve
import znajdz_kat

def wczytaj_elementy_orbitalne(nazwaPliku):
	dElementyOrbitalne_linie = opr.wczytaj_dane(nazwaPliku)					#Wczytuje linie tekstu z pliku wejściowego
	dElementyOrbitalne_bez_komentarzy = opr.usun_komentarze(dElementyOrbitalne_linie)	#Usuwa komentarze
	dElementyOrbitalne = opr.uporzadkuj_z_rozdz_tekst(dElementyOrbitalne_bez_komentarzy)	#Rozdziela dane na dwie listy: pierwszą z nazwami ciał niebieskich, a drugą z ich el. orbitalnymi (jeszcze jako stringi)
	return dElementyOrbitalne

def wczytaj_nazwy_obiektow(nazwaPliku):
	nazwyObiektow = wczytaj_elementy_orbitalne(nazwaPliku)[0]
	return nazwyObiektow

def wczytaj_wartosci_elementow_orbitalnych(nazwaPliku):
	dElementyOrbitalne = wczytaj_elementy_orbitalne(nazwaPliku)[1]
	n = len(dElementyOrbitalne)
	elementyOrbitalne = []
	for i in range(n):
		elementyOrbitalne[i:i] = [[float(k) for k in dElementyOrbitalne[i]]]
	return elementyOrbitalne

def przypisz_czas_obserwacji(czas):
	"""
	Zapisuje żądany czas UTC obserwacji w formie listy [r, m, d, g, m, s], w której liczba oznaczająca sekundę jest typu float, natomiast pozostałe są typu integer.
	"""
	t_UTC_int = czas[0:5]
	t_UTC_float = czas[5]
	t_UTC = [int(k) for k in t_UTC_int]
	t_UTC[6:6] = [t_UTC_float]
	return t_UTC

def przelicz_deg_na_rad(kat):
	kat_rad = kat/180.0*m.pi
	return kat_rad

def przelicz_rad_na_deg(kat):
	kat_deg = kat/m.pi*180
	return kat_deg

def przelicz_rad_na_h(kat):
	kat_h = kat/m.pi*12
	return kat_h

def przelicz_czas_na_JD(czas):
	"""
	Przelicza datę i godzinę z kalendarza gregoriańskiego na datę juliańską. Argument funkcji musi być listą o strukturze [rok, miesiąc, dzień, godzina, minuta, sekunda]. Źródło algorytmu: angielska Wikipedia pod hasłem "Julian day".
	"""
	rok = czas[0]
	miesiac = czas[1]
	dzien = czas[2]
	godzina = czas[3]
	minuta = czas[4]
	sekunda = czas[5]
	a = int((14 - miesiac)/12.)	#Styczeń i luty to w JD miesiące 13. i 14. poprzedniego roku; zmienna jest przełącznikiem 0-1 w przypadku wystąpienia tych dwóch miesięcy
	y = rok + 4800 - a		#Numer roku od początku rachuby
	m = miesiac + 12*a - 3		#Numer miesiąca od 0 (marzec) do 11 (luty)
	JDN = dzien + int((153*m + 2)/5.) + 365*y + int(y/4.) - int(y/100.) + int(y/400.) - 32045	#Drugi od lewej składnik sumy bierze się z różnej ilości dni w poszczególnych miesiącach
	JD = JDN + (godzina - 12)/24. + minuta/1440. + sekunda/86400.
	return JD

def przelicz_JD_na_JED(JD):
	"""
	Przelicza czas JD (Julian Date), odpowiadający czasowi uniwersalnemu UTC, na czas JED (Julian Ephemeris Date), odpowiadający czasowi ziemskiemu TT.
	"""
	Dt = 69.184		#Różnica w sekundach pomiędzy czasami TT i UTC na 2017 rok
	JED = JD + Dt/86400.
	return JED

def wczytaj_pozycje_ziemi(nazwa_pliku):
	dPozycjeZiemi_linie = opr.wczytaj_dane(nazwa_pliku)
	dPozycjeZiemi_bez_komentarzy = opr.usun_komentarze(dPozycjeZiemi_linie)
	dPozycjeZiemi = opr.uporzadkuj(dPozycjeZiemi_bez_komentarzy)
	return dPozycjeZiemi

def interpoluj_wspolrzedna_ziemi(nazwa_pliku, t, nr_wspolrzednej):
	"""
	Interpoluje jedną, wybraną współrzędną położenia Ziemi w układzie barycentrycznym.
	"""
	dPozycjeZiemi = wczytaj_pozycje_ziemi(nazwa_pliku)
	t_min = int(t) + 0.5	#Zakładam tutaj, że dane w pliku z położeniami Ziemi wyznaczone są na godzinę 00:00 każdego kolejnego dnia!
	t_max = int(t) + 1.5	#j.w.
	wsp = opr.znajdz_wartosci_po_wierszach(dPozycjeZiemi, [t_min, t_max], nr_wspolrzednej)
	a = (wsp[1] - wsp[0])/(t_max - t_min)
	b = (wsp[0]*t_max - wsp[1]*t_min)/(t_max - t_min)
	wspolrzedna = a*t + b
	return wspolrzedna

def interpoluj_pozycje_ziemi(nazwa_pliku, t):
	"""
	Wywołuje interpolację dla wszystkich trzech współrzędnych położenia Ziemi w układzie barycentrycznym.
	"""
	X_z = interpoluj_wspolrzedna_ziemi(nazwa_pliku, t, 1)
	Y_z = interpoluj_wspolrzedna_ziemi(nazwa_pliku, t, 2)
	Z_z = interpoluj_wspolrzedna_ziemi(nazwa_pliku, t, 3)
	return [X_z, Y_z, Z_z]

def oblicz_okres(a):
	d = 2*m.pi/0.01720209895	#Ilość dób słonecznych w roku gwiazdowym (rok zwrotnikowy, ze względu na precesję Ziemi, jest krótszy od gwiazdowego o ok. 20.4 min)
	T = m.sqrt(a**3)*d		#Okres wyrażony w dniach; wzór ma postać T = 2*pi*sqrt(a**3/(G*M_s)), ale upraszcza się, gdy a wyrażone jest w AU, ponieważ G = 4*pi**2 AU**3 yr**(-2) M_s**(-1).
	return T

def rownanie_na_anomalie_mimosrodowa(E, e, t, t_p, T, h):
	rownanie = E - e*m.sin(E) - 2*m.pi*(t - h - t_p)/T
	return rownanie

def wyznacz_anomalie_mimosrodowa(e, t, t_p, T, h):
	"""
	Rozwiązuje numerycznie równanie na anomalię mimośrodową E. E powinno znajdować się w przedziale od -pi do +pi.
	"""
	E = fsolve(rownanie_na_anomalie_mimosrodowa, 0.0, (e, t, t_p, T, h))
	return E[0]

def wyznacz_anomalie_prawdziwa(E, e):
	phi = 2*m.atan(m.sqrt((1.0 + e)/(1.0 - e))*m.tan(E/2.0))	#Prosty wzór z arctg, ponieważ phi znajduje się w przedziale od -pi do +pi
	return phi

def wyznacz_odleglosc_od_barycentrum(a, e, E):
	R = a*(1 - e*m.cos(E))
	return R

def wyznacz_polozenie_barycentryczne(R, Omega, phi, omega, i):
	X = R*(m.cos(Omega)*m.cos(phi + omega) - m.sin(Omega)*m.cos(i)*m.sin(phi + omega))
	Y = R*(m.sin(Omega)*m.cos(phi + omega) + m.cos(Omega)*m.cos(i)*m.sin(phi + omega))
	Z = R*m.sin(i)*m.sin(phi + omega)
	return [X, Y, Z]

def wyznacz_polozenie_geocentryczne(X, X_z):
	epsilon_deg = 23.439292				#nachylenie ekliptyki do równika [deg]
	epsilon = przelicz_deg_na_rad(epsilon_deg)	#[rad]
	x = X[0] - X_z[0]
	y = (X[1] - X_z[1])*m.cos(epsilon) - (X[2] - X_z[2])*m.sin(epsilon)
	z = (X[1] - X_z[1])*m.sin(epsilon) + (X[2] - X_z[2])*m.cos(epsilon)
	return [x, y, z]

def wyznacz_odleglosc_od_ziemi(x):
	r = m.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
	return r

def wyznacz_opoznienie_sygnalu(r):
	c_ms = 299792458.0	#prędkość światła w próżni [m/s]
	d = 86400.0		#przelicznik s na d [s/d]
	a = 149597870700.0	#przelicznik m na AU [m/AU]
	c = c_ms*d/a		#[AU/d]
	h = r/c
	return h

def wyznacz_rektascensje_i_deklinacje(x, r):
	delta_rad = m.asin(x[2]/r)					#Deklinacja zmienia się w przedziale od -pi/2 do +pi/2, zatem można użyć zwykłej funkcji odwrotnej sinusa; [rad]
	cos_alpha = x[0]/(r*m.cos(delta_rad))
	sin_alpha = x[1]/(r*m.cos(delta_rad))
	alpha_rad = znajdz_kat.z_funkcji_tryg(sin_alpha, cos_alpha)	#Rektascensja zmienia się w przedziale od 0 do 2*pi, więc nie można użyć zwykłej funkcji odwrotnej; [rad]
	alpha = przelicz_rad_na_h(alpha_rad)				#[arch]
	delta = przelicz_rad_na_deg(delta_rad)				#[deg]
	return [alpha, delta]

def zapisz_efemerydy(nazwaPliku, nazwyObiektow, czasyObserwacji, odleglosciOdZiemi, wspolrzedneRownikowe):
	dZapis = opr.zlacz_dane_z_jednej_linii(nazwyObiektow, czasyObserwacji, odleglosciOdZiemi, wspolrzedneRownikowe)
	dZapis[0:0] = [["#nazwa", "UTC(r", "m", "d", "g", "m", "s)", "odleglosc [AU]", "RA [arch]", "DEC [deg]"]]
	opr.zapisz_dane(nazwaPliku, dZapis)

if __name__ == "__main__":
	dokl = 12							#Dokładność, z jaką mają być porównywane wartości kolejnych wyznaczonych anomalii mimośrodowych E
	fElementyOrbitalne = argv[1]					#Pobranie nazwy pliku wejściowego z elementami, która jest argumentem wywołania programu
	fPozycjaZiemi = "pozycja_ziemi.dat"				#Nazwa pliku wejściowego ze współrzędnymi Ziemi w układzie barycentrycznym
	fEfemerydy = "efemerydy.dat"			 		#Nazwa pliku wyjściowego z efemerydami
	nazwyObiektow = wczytaj_nazwy_obiektow(fElementyOrbitalne)
	elementyOrbitalne = wczytaj_wartosci_elementow_orbitalnych(fElementyOrbitalne)
	czasyObserwacji = []
	odleglosciOdZiemi = []
	wspolrzedneRownikowe = []
	n = len(elementyOrbitalne)
	for k in range(n):							#Pętla po obiektach z listy
		e = elementyOrbitalne[k][0]					#mimośród
		a = elementyOrbitalne[k][1]					#półoś wielka [AU]
		t_p = elementyOrbitalne[k][2]					#czas momentu przejścia przez perycentrum [JD]
		Omega = przelicz_deg_na_rad(elementyOrbitalne[k][3])		#długość węzła wstępującego [rad]
		i = przelicz_deg_na_rad(elementyOrbitalne[k][4])		#inklinacja [rad]
		omega = przelicz_deg_na_rad(elementyOrbitalne[k][5])		#argument perycentrum [rad]
		t_UTC = przypisz_czas_obserwacji(elementyOrbitalne[k][6:12])	#żądany czas UTC obserwacji [r, m, d, g, m, s]
		t_JD = przelicz_czas_na_JD(t_UTC)				#[d]
		t = przelicz_JD_na_JED(t_JD)					#[d]
		X_z = interpoluj_pozycje_ziemi(fPozycjaZiemi, t)		#[AU, AU, AU]
		T = oblicz_okres(a)						#[d]
		h = 0.0								#czas dotarcia sygnału świetlnego od ciała niebieskiego do Ziemi --- poprawka do wzoru na anomalię mimośrodową
		E_kontrolna = 0.0						#anomalia mimośrodowa
		koniec = False
		while not koniec:
			E = wyznacz_anomalie_mimosrodowa(e, t, t_p, T, h)		#[rad]
			phi = wyznacz_anomalie_prawdziwa(E, e)				#[rad]
			R = wyznacz_odleglosc_od_barycentrum(a, e, E)			#[AU]
			X = wyznacz_polozenie_barycentryczne(R, Omega, phi, omega, i)	#[AU, AU, AU]
			x = wyznacz_polozenie_geocentryczne(X, X_z)			#[AU, AU, AU]
			r = wyznacz_odleglosc_od_ziemi(x)				#[AU]
			h = wyznacz_opoznienie_sygnalu(r)				#[d]
			if round(E - E_kontrolna, dokl) == 0.0:				#Sprawdzenie, czy poprawki mają wpływ na wartość E
				koniec = True
			E_kontrolna = E
		wspolrzedneRownikowe_k = wyznacz_rektascensje_i_deklinacje(x, r)	#[arch, deg]
		czasyObserwacji[k:k] = [t_UTC]
		odleglosciOdZiemi[k:k] = [r]
		wspolrzedneRownikowe[k:k] = [wspolrzedneRownikowe_k]
	zapisz_efemerydy(fEfemerydy, nazwyObiektow, czasyObserwacji, odleglosciOdZiemi, wspolrzedneRownikowe)
