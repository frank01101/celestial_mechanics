#!/usr/bin/env python
#-*- coding: utf-8 -*-
#autor: Franciszek Humieja
#wersja 1.0 (17.XI.2014)

import math as m

dokl_arg = 12

def walcz_z_dziadostwem(a, b):
  """
  Walczy z dziadostwem (gdy funkcje trygonometryczne przyjmują wartości większe niż 1 lub mniejsze niż -1).
  """
  if(round(a, dokl_arg) > 1):
    a = 1.0
  elif(round(a, dokl_arg) < -1):
    a = -1.0
  if(round(b, dokl_arg) > 1):
    b = 1.0
  elif(round(b, dokl_arg) < -1):
    b = -1.0
  return a, b

def z_funkcji_tryg(inc_sinus, inc_cosinus):
  """
  Znajduje kąt w przedziale [0,2*pi) gdy znane są funkcje sin oraz cos tego kąta.
  """
  sinus, cosinus = walcz_z_dziadostwem(inc_sinus, inc_cosinus)
  dokl = 6
  while(dokl >= 0):
    #W przedziale [0, 2*pi) istnieją dwa argumenty x, dla których funkcje sin oraz cos mają daną wartość. Poniżej program wyszukuje te dwie wartości używając funkcji trygonometrycznych odwrotnych.
    if sinus >= 0:
      x1 = set([round(m.asin(round(sinus, dokl_arg)), dokl), round(m.pi-m.asin(round(sinus, dokl_arg)), dokl)])
    if sinus < 0:
      x1 = set([round(m.asin(round(sinus, dokl_arg))+2*m.pi, dokl), round(m.pi-m.asin(round(sinus, dokl_arg)), dokl)])
    x2 = set([round(m.acos(round(cosinus, dokl_arg)), dokl), round(2*m.pi-m.acos(round(cosinus, dokl_arg)), dokl)])
    if not x1.isdisjoint(x2): #Jeżeli zbiory x1 i x2 mają wspólny element...
      x = x1 & x2 #Część wspólna zbiorów x1 i x2 to szukana wartość kąta.
      return list(x)[0]
    else:
      dokl -= 1 #Dokładność funkcji trygonometrycznych w Pythonie jest niska, dlatego gdy zbiory x1 i x2 nie posiadają wspólnego elementu należy zmniejszyć dokładność porównywania.
