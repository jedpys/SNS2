# -*- coding: utf-8 -*-
"""

@author: Maciek


"""
from readrnx_studenci import readrnxnav, date2tow, readrnxobs
import numpy as np
import math as m

nav_file = "BRDC00WRD_R_20240650000_01D_GN.rnx"

nav, inav = readrnxnav(nav_file)

zdrowe = nav[:, 30] == 0
inav = inav[zdrowe]
nav = nav[zdrowe, :]


sat = 2
indeks_satelity = inav == 2
nav_satelity = nav[indeks_satelity, :]

data_obliczen = [2024, 3, 5, 11, 15, 0]
week, tow, dow = date2tow(data_obliczen)

t_z_tygodniami = week * 7 * 86400 + tow


toe = nav_satelity[:, 17]

# gps_week = nav_satelity[:, 27]
# toe_z_tygodniami = gps_week

roznica = np.abs(tow - toe)

id_min_roznica = np.argmin(roznica)


nav_satelity_wybrana_epoka = nav_satelity[id_min_roznica]

tst = toe[7]

mi = 3.986005 * (10 ** 14)
omegaE = 7.2921151467 * (10 ** -5)
c = 299792458

def satpos(week, tow, nav_sat):
    t_full = tow + week * 7 * 24 * 60 * 60
    toe_full = nav_sat[17] + nav_sat[27] * 7 * 24 * 60 * 60
    tk = t_full - toe_full

    a = nav_sat[16] ** 2
    n0 = np.sqrt(mi/(a ** 3))
    n = n0 + nav_sat[11]
    Mk = nav_sat[12] + n * tk

    E1 = Mk
    E2 = Mk + nav_sat[14] * np.sin(Mk)
    while abs(E2 - E1) > pow(10, -12):
        E1 = E2
        E2 = Mk + nav_sat[14] * m.sin(E2)
    Ek = E2

    vk = np.arctan2((np.sqrt(1 - nav_sat[14] ** 2) * np.sin(Ek)), np.cos(Ek) - nav_sat[14])

    fik = vk + nav_sat[23]

    delta_uk = nav_sat[15] * np.sin(2 * fik) + nav_sat[13] * np.cos(2 * fik)
    delta_rk = nav_sat[10] * np.sin(2 * fik) + nav_sat[22] * np.cos(2 * fik)
    delta_ik = nav_sat[20] * np.sin(2 * fik) + nav_sat[18] * np.cos(2 * fik)


    uk = fik + delta_uk
    rk = a * (1 - nav_sat[14] * np.cos(Ek)) + delta_rk
    ik = nav_sat[21] + nav_sat[25] * tk + delta_ik

    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)


    Omega_k = nav_sat[19] + (nav_sat[24] - omegaE) * tk - omegaE * nav_sat[17]

    Xk = xk * np.cos(Omega_k) - yk * np.cos(ik) * np.sin(Omega_k)
    Yk = xk * np.sin(Omega_k) + yk * np.cos(ik) * np.cos(Omega_k)
    Zk = yk * np.sin(ik)

    #if np.abs(rk - np.sqrt(Xk ** 2 + Yk ** 2 + Zk ** 2)) >= 0.01:
    #    print("bład")

    delta_ts = nav_sat[6] + nav_sat[7] * (tow - nav_sat[17]) + nav_sat[8] * ((tow - nav_sat[17]) ** 2)
    delta_trel = (-2 * np.sqrt(mi) * nav_sat[14] * np.sqrt(a) * np.sin(Ek))/(c ** 2)
    delta_trels = delta_ts + delta_trel

    return ((Xk, Yk, Zk), delta_trels)


def uklad_chwilowy_trans(xyz, new_tau):
    tau = new_tau
    R = np.column_stack([[np.cos(omegaE*tau), -np.sin(omegaE*tau), 0],
                        [np.sin(omegaE*tau), np.cos(omegaE*tau), 0],
                        [0, 0, 1]])
    new_xyz = R * np.array(xyz)
    return new_xyz



# cieżka do pliku nawigacyjnego
nav_file = 'BRDC00WRD_R_20240650000_01D_GN.rnx'
# cieżka do pliku obserwacyjnego
obs_file = 'JOZ200POL_R_20240650000_01D_30S_MO.rnx'

# zdefiniowanie czasu obserwacji: daty początkowej i końcowej
# dla pierwszej epoki z pliku będzie to:
time_start =  [2024, 3, 5, 0, 0, 0]
time_end =    [2024, 3, 5, 23, 59, 59]

# data_obliczen = [2024]
# odczytanie danych z pliku obserwacyjnego
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
# odczytanie danych z pliku nawigacyjnego:
nav, inav = readrnxnav(nav_file)

# print(obs)
# print(iobs)


#%%
"""
zdefiniowanie współrzędnych przybliżonych odbiornika - mogą to być współrzędne z nagłówka 
pliku obserwacyjnego, skopiowane "z palca" lub pobierane automatycznie z treci nagłówka pliku Rinex
"""
xr0 = [3660000.,  1400000.,  5000000.]

"""
Wprowadzenie ustawień, takich jak maska obserwacji, czy typ poprawki troposferycznej
"""
el_mask = 10 # elevation mask/cut off in degrees

"""
Przeliczenie daty początkowej i końcowej do sekund tygodnia GPS - niezbędne w celu
poprawnej definicji pętli związanej z czasem obserwacji w ciągu całej doby
"""
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]
#%% Obliczenia



"""
Otwieramy dużą pętlę
dt = 30
for t in range(tow, tow_end+1, dt): gdzie dt równe 30
"""
t = 213300  #czas odbiornika - 11:15

index_t = iobs[:, 2] == t  #wybieramy z iobs tylko obserwacje zgodne z czasem odbiornika t
Pobs = obs[index_t, 0]  #w Pobs tylko interesujące nas pseudoodległości

sats = iobs[index_t, 0]  #punkt 2 instrukcja z kodu  ##wczesniej t zamaist index_t
print(sats)
dtr = 0
tau = 0.07

'''
Wewnątrz tej pętli, zajmujemy się obserwacjami wyłącznie dla jednej epoki (epoka t), zatem:
    1. Wybieramy obserwacje dla danej epoki, na podstawie tablicy iobs oraz naszej epoki t
    czyli, wybieramy te obserwacje z tablicy obs, dla których w tablicy iobs ostatnia kolumna 
    jest równa t - przypisujemy do zmiennej np. Pobs
    2. wybieramy satelity, obserwowane w danej epoce, na podstawie tablicy iobs - na podstawie 
    naszego t - przypisujemy do zmiennej np. sats
    3. Definiujemy wartości przybliżone błąd zegara odbiornika
    dtr = 0 oraz czasu propagacji sygnału tau = 0.07
    4. Najprawdopodobniej przyda się definicja pustych wektorów, np. zawierających 
    odległosci geometryczne (wartoci przybliżone na podstawie tau)     
    Przechodzimy do iteracyjnego obliczenia współrzędnych odbiornika - w pierwszych testach naszego programu, zróbmy obliczenia nieiteracyjnie, 
    ale pamiętajmy o tym, że będzie trzeba przygotować kod do działania w pętli:
        
        Po weryfikacji działania programu, można zamienić pętlę for na pętle while, dopisując
        warunek zbieżnoci kolejnych współrzędnych - skróci nam to czas obliczeń, ponieważ 
        najczęściej wystarcza mniej iteracji niż 5
'''
for i in range(1):  #potem zamień na 5
    for sat in range(len(sats)):
        print(sat)
        ts = t - tau + dtr  #czas emisji sygnału satelity
        xyzs, dts = satpos(week, ts, nav_satelity[sat])
        xyzs_chwil = uklad_chwilowy_trans(xyzs, tau)



'''
        
        for i in range(5):
            Wykonujemy kolejne obliczenia, niezależnie dla kolejnych satelitów, obserwowanych
            w danej epoce, czyli przechodzimy do pętli:
                for sat in sats: (przyda nam się tutaj również indeks satelity, np. for i, sat in enumerate(sats):)
                    Obliczamy czas emisji sygnału:
                        ts = t - tau + dtr
                    Kolejne kroki, znane z poprzedniego ćwiczenia:
                    wyznaczamy współrzędne satelity xs (oraz błąd zegara satelity dts) na czas ts (UWAGA, w kolejnych iteracjach
                    czas ts będzie się zmieniał i aktualizował, niezależnie dla każdego satelity!!!)
                    
                    Odległosć geometryczna:
                        1. rotacja do układu chwilowego - otrzymujemy xs_rot
                        2. Na podstawie xs_rot obliczamy odległosć geometryczną rho
                        
                    Obliczamy elewację i azymut
                    Macierz Rneu definiujemy na podstawie x0, przeliczonego do współrzędnych
                    phi, lambda, algorytmem Hirvonena
                    
                    Odrzucamy satelity znajdujące się poniżej maski
                    
                        Obliczamy poprawki atmosferyczne - dopiero wówczas, kiedy działać będzie nam program bez uwzględniania poprawek:
                            trop oraz iono
                    
                    Wyznaczamy pseudoodległosć przybliżoną (obliczoną), jako:
                        Pcalc = rho - cdts + dtr + trop + iono
                        
                    Wyznaczamy kolejne elementy wektora wyrazów wolnych y, jako:
                        y = Pobs - Pcalc
                        
                    Budujemy kolejne wiersze macierzy A:
                
                Kończymy pętle dla kolejnych satelitów
                
                1. Łączymy ze sobą elementy wektora wyrazów wolych w jeden wektor
                2. Łączymy ze sobą kolejnę wiersze macierz współczynników A
                3. Rozwiązujemy układ równań, metodą najmniejszych kwadratów
                
                               
                Aktualizujemy wartosci przybliżone o odpowiednie elementy wektora x
                xr[0] = x0[0] + x[0]
                xr[1] = x0[1] + x[1]
                xr[2] = x0[2] + x[2]
                dtr = dtr + x[3]/c 
                
                Tak obliczone wartoci xr oraz dtr stanowią wartoci wejsciowe do kolejnej iteracji, itd 
                do skończenia piątej iteracji lub spełnienia warunku zbieżnoci współrzędncyh
            
            
            Po skończeniu 5. iteracji, zbieramy obliczone współrzędne xr - warto zebrać również
            liczby obserwowanych satelitów, obliczone wartoci współczynników DOP (przynajmniej PDOP)
            
'''








