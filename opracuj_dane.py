#!/usr/bin/env python
#-*- coding: utf-8 -*-
# autor: Franciszek Humieja
# wersja 1.1 --- 2017-01-22

from sys import argv
from math import sqrt

#Określa, do którego miejsca po przecinku mają być zaokrąglane wyliczane warości.
zaokr = 6

def wczytaj_dane(nazwa_pliku):
   """
   Wczytuje dane z pliku wejściowego i zwraca listę, której każdy z elementow jest jedną linijką z pliku wejściowego.
   """
   plik = open(nazwa_pliku, "r")
   dane = plik.readlines()
   plik.close()
   return dane

def uporzadkuj(lista_danych):
   """
   Porządkuje dane: zwraca dwuwymiarową listę złożona z liczb o takim samym układzie jak w pliku wejściowym.
   """
   n = len(lista_danych)
   lista_uporz = []
   for i in range(n):
      lista_uporz[i:i] = [lista_danych[i].split()]
      lista_uporz[i] = [float(k) for k in lista_uporz[i]]
   return lista_uporz

def uporzadkuj_tekst(lista_danych):
   """
   Wersja funkcji uporzadkuj dla danych, w których może wystąpic tekst.
   Porządkuje dane: zwraca dwuwymiarową listę złożona ze słów o takim samym układzie jak w pliku wejściowym (nie konwertuje tekstu na liczbę).
   """
   n = len(lista_danych)
   lista_uporz = []
   for i in range(n):
      lista_uporz[i:i] = [lista_danych[i].split()]
   return lista_uporz

def uporzadkuj_z_rozdz(lista_danych, czy_trzecia_dana = False):
   """
   Porządkuje dane:
   1) Liczby z pierwszej kolumny pliku wejściowego umieszcza w liście d1.
   2) Jeżeli użytkownik oznajmił, że w ostatniej kolumnie znajduje się trzeci rodzaj wielkości to:
      2.1) Liczby z kolumn: od drugiej do przedostatniej pliku wejściowego umieszcza w dwuwymiarowej liście d2.
      2.2) Liczby z ostatniej kolumny pliku wejściowego umieszcza w liscie d3.
   3) Jeżeli użytkownik oznajmił, że w ostatniej kolumnie nie znajduje się trzeci rodzaj wielkości to:
      2.1) Liczby z kolumn: od drugiej do ostatniej pliku wejściowego umieszcza w dwuwymiarowej liście d2.
   4) Zwraca listy d1, d2, d3.
   """
   lista_uporz = uporzadkuj(lista_danych)
   n = len(lista_uporz)
   d1 = []
   d2 = []
   d3 = []
   for i in range(n):
      m = len(lista_uporz[i])
      d1[i:i] = [lista_uporz[i][0]]
      if czy_trzecia_dana:
         d2[i:i] = [lista_uporz[i][1:m-1]]
         d3[i:i] = [lista_uporz[i][m-1]]
      else:
         d2[i:i] = [lista_uporz[i][1:m]]
   return d1, d2, d3

def uporzadkuj_z_rozdz_tekst(lista_danych, czy_trzecia_dana = False):
   """
   Wersja funkcji uporzadkuj_z_rozdz dla danych, w których może wystąpić tekst.
   Porządkuje dane:
   1) Znaki z pierwszej kolumny pliku wejściowego umieszcza w liście d1.
   2) Jeżeli użytkownik oznajmił, że w ostatniej kolumnie znajduje się trzeci rodzaj wielkości to:
      2.1) Znaki z kolumn: od drugiej do przedostatniej pliku wejściowego umieszcza w dwuwymiarowej liście d2.
      2.2) Znaki z ostatniej kolumny pliku wejściowego umieszcza w liscie d3.
   3) Jeżeli użytkownik oznajmił, że w ostatniej kolumnie nie znajduje się trzeci rodzaj wielkości to:
      2.1) Znaki z kolumn: od drugiej do ostatniej pliku wejściowego umieszcza w dwuwymiarowej liście d2.
   4) Zwraca listy d1, d2, d3.
   """
   lista_uporz = uporzadkuj_tekst(lista_danych)
   n = len(lista_uporz)
   d1 = []
   d2 = []
   d3 = []
   for i in range(n):
      m = len(lista_uporz[i])
      d1[i:i] = [lista_uporz[i][0]]
      if czy_trzecia_dana:
         d2[i:i] = [lista_uporz[i][1:m-1]]
         d3[i:i] = [lista_uporz[i][m-1]]
      else:
         d2[i:i] = [lista_uporz[i][1:m]]
   return d1, d2, d3

def usun_komentarze(lista_danych, znak_komentarza = '#'):
   """
   Usuwa wszystkie znaki począwszy od znaku, podanego jako znak komentarza, aż do końca linii.
   """
   n = len(lista_danych)
   lista_bez_kom = []
   for i in range(n):
      ciecie = len(lista_danych[i])
      if znak_komentarza in lista_danych[i]:
         m = len(lista_danych[i])
         for k in range(m):
            if lista_danych[i][k] == znak_komentarza and ciecie >= k:
               ciecie = k
      lista_bez_kom[i:i] = [lista_danych[i][0:ciecie]]
   lista_wyczyszczona = [k for k in lista_bez_kom if k != '']
   return lista_wyczyszczona

def stworz_nazwe_pliku_out(nazwa_pliku_in, rozszerzenie_pliku_out = ".out"):
   """
   Zwraca nazwę pliku wyjściowego. W przypadku gdy nazwa pliku wejściowego posiada rozszerzenie to nazwa pliku wyjściowego jest taka sama z zamienionym rozszerzeniem na ciąg znaków rozszerzenie_pliku_out. Gdy nazwa pliku wejściowego nie posiada rozszerzenia (bądź posiada rozszerzenie takie samo jak rozszerzenie_pliku_out) to nazwa pliku wyjściowego jest taka sama z dodanym na końcu ciągiem znaków rozszerzenie_pliku_out.
   """
   n = len(nazwa_pliku_in)
   if "." in nazwa_pliku_in and rozszerzenie_pliku_out not in nazwa_pliku_in:
      for i in range(1, n+1):
         if nazwa_pliku_in[-i] == ".":
            nazwa_pliku_out = nazwa_pliku_in[0:n-i] + rozszerzenie_pliku_out
            return nazwa_pliku_out
   else:
      nazwa_pliku_out = nazwa_pliku_in + rozszerzenie_pliku_out
      return nazwa_pliku_out

def zlacz_dane_z_jednej_linii(*listy_danych):
   """
   Zwraca listę dwuwymiarową, której elementy i-tego wiersza będą zapisane w i-tej linijce w pliku wyjściowym.
   """
   ilosc_przekazanych_list = len(listy_danych)
   n = len(listy_danych[0])
   d = [listy_danych[i] for i in range(ilosc_przekazanych_list) if len(listy_danych[i]) == n]   
   m = len(d)
   wartosci_linia = []
   for i in range(n):
      wartosci_linia[i:i] = [[]]
      for k in range(m):
         if type(d[k][i]) != list:
            wartosci_linia[i] += [d[k][i]]
         else:
            wartosci_linia[i] += d[k][i]
      wartosci_linia[i] = [str(l) for l in wartosci_linia[i]]
   return wartosci_linia

def zapisz_dane(nazwa_pliku_out, lista_z_wierszami, separator = "\t"):
   """
   Łączy wszystkie elementy i-tego wiersza listy lista_z_wierszami w jeden ciąg znaków (elementy oddziela od siebie ustalonym separatorem) i zapisuje go do pliku. Powtarza procedurę dla wszystkich n wierszy listy.
   """
   n = len(lista_z_wierszami)
   plik_out = open(nazwa_pliku_out, "w")
   for i in range(n):
      linia = separator.join(lista_z_wierszami[i])
      plik_out.write(linia)
      plik_out.write("\n")
   plik_out.close()

def ustal_czy_wyst_opcja(tresc_opcji, od = 2, argumenty = argv):
   """
   Zwraca True, gdy jako od-ty lub dalszy argument wywołania programu podaliśmy opcję o treści tresc_opcji.
   """
   n = len(argumenty)
   czy_wyst = False
   for i in range(od, n):
      if argumenty[i] == tresc_opcji:
         czy_wyst = True
   return czy_wyst

def znajdz_argument_opcji(tresc_opcji, od = 2, argumenty = argv):
   """
   Zwraca wartość argumentu opcji tresc_opcji. Użytkownik powinien wpisywać argumenty po opcjach im odpowiadających (np: -t 0.95).
   """
   n = len(argumenty)
   for i in range(od, n-1):
      if argumenty[i] == tresc_opcji:
         arg_opcji = argumenty[i+1]
   return arg_opcji

def wyznacz_dlugosci_podlist(lista):
   """
   Zwraca listę, w której i-ty element jest długością podlisty lista[i].
   """
   n = len(lista)
   dlugosci_podlist = [len(lista[i]) for i in range(n)]
   return dlugosci_podlist

def oblicz_srednie(lista):
   """
   Zwraca listę, w której i-ty element jest średnią arytmetyczną liczb znajdujących się w podliście lista[i].
   """
   n = len(lista)
   srednie = []
   for i in range(n):
      m = len(lista[i])
      suma = 0.0
      for k in range(m):
         suma += lista[i][k]
      srednie[i:i] = [round(float(suma)/m, zaokr)]
   return srednie

def oblicz_odchylenia_std(lista, wartosci_srednie):
   """
   Zwraca listę, w której i-ty element jest odchyleniem standardowym średniej arytmetycznej liczb znajdujących się w podliście lista[i].
   """
   n = len(lista)
   odchylenia = []
   for i in range(n):
      m = len(lista[i])
      w = 0.0
      for k in range(m):
         w += (lista[i][k]-wartosci_srednie[i])**2
      odchylenia[i:i] = [round(sqrt(float(w)/(m*(m-1))), zaokr)]
   return odchylenia

def oblicz_niep(odchylenia_std, niep_graniczne = []):
   """
   Zwraca wartość niepewności systematycznej.
   """
   n = len(odchylenia_std)
   if len(niep_graniczne) != 0:
      niep = [round(sqrt(odchylenia_std[i]**2+niep_graniczne[i]**2/3), zaokr) for i in range(n)]
   else:
      niep = odchylenia_std
   return niep

def oblicz_niep_rozsz(wartosci_wsp_t, niep):
   """
   Zwraca wartość niepewności rozszerzonej.
   """
   n = len(niep)
   niep_rozsz = [round(wartosci_wsp_t[i]*niep[i], zaokr) for i in range(n)]
   return niep_rozsz

def znajdz_kolumne(lista, wartosc_szukana):
   """
   Zwraca numer kolumny, na początku której znajduje się wartość najbliższa zadanej wartości z zadanej jako argument listy.
   """
   n = len(lista[0])
   min_roznica = 1.0
   kol = n-1
   for i in range(n):
      if abs(round(wartosc_szukana-lista[0][i], zaokr)) <= min_roznica:
         min_roznica = abs(round(wartosc_szukana-lista[0][i], zaokr))
         kol = i
   return kol

def znajdz_wartosci_po_wierszach(lista, wiersze_szukane, nr_kol):
   """
   Szuka, na początku których wierszy zadanej listy znajdują się wartości z listy wiersz_szukany, a następnie zwraca listę elementów z tych wierszy, których druga współrzędna jest równa nr_kol.
   """
   n = len(wiersze_szukane)
   m = len(lista)
   wartosci = []
   for i in range(n):
      wartosci[i:i] = [0.0]
      for k in range(m):
         if wiersze_szukane[i] == lista[k][0]:
            wartosci[i] = lista[k][nr_kol]
   return wartosci

if __name__ == "__main__":
   nazwa_pliku_in = argv[1]
   nazwa_pliku_out = stworz_nazwe_pliku_out(nazwa_pliku_in)
   nazwa_pliku_rozkl_stud = "tstudent.dat"
   tresc_opcji_niep_graniczna = "-g"
   tresc_opcji_rozkl_stud = "-t"
   czy_niep_graniczna = ustal_czy_wyst_opcja(tresc_opcji_niep_graniczna)
   czy_rozkl_stud = ustal_czy_wyst_opcja(tresc_opcji_rozkl_stud)
   dane_in = wczytaj_dane(nazwa_pliku_in)
   y, x, dx = uporzadkuj_z_rozdz(dane_in, czy_niep_graniczna)
   x_sr = oblicz_srednie(x)
   sx = oblicz_odchylenia_std(x, x_sr)
   ux = oblicz_niep(sx, dx)
   uux = []
   if czy_rozkl_stud:
      podany_poz_ufnosci = float(znajdz_argument_opcji(tresc_opcji_rozkl_stud))
      dane_rozkl_stud = wczytaj_dane(nazwa_pliku_rozkl_stud)
      lista_rozkl_stud = uporzadkuj(dane_rozkl_stud)
      nr_kolumny_wsp_t = znajdz_kolumne(lista_rozkl_stud, podany_poz_ufnosci)
      ilosci_x = wyznacz_dlugosci_podlist(x)
      t = znajdz_wartosci_po_wierszach(lista_rozkl_stud, ilosci_x, nr_kolumny_wsp_t)
      uux = oblicz_niep_rozsz(t, ux)
   dane_out = zlacz_dane_z_jednej_linii(y, x, dx, x_sr, sx, ux, uux)
   zapisz_dane(nazwa_pliku_out, dane_out)
