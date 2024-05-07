from readrnx_studenci import readrnxnav, date2tow
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
print(week, tow, dow)

t_z_tygodniami = week * 7 * 86400 + tow
print(t_z_tygodniami)

toe = nav_satelity[:, 17]

print(toe)
# gps_week = nav_satelity[:, 27]
# toe_z_tygodniami = gps_week

roznica = np.abs(tow - toe)
print(roznica)
id_min_roznica = np.argmin(roznica)
print("id_min", id_min_roznica)

nav_satelity_wybrana_epoka = nav_satelity[id_min_roznica]
# print(nav_satelity_wybrana_epoka)
tst = toe[7]
print("tst", tst, "\n\n\n\n")
mi = 3.986005 * (10 ** 14)
omegaE = 7.2921151467 * (10 ** -5)
c = 299792458

def satpos(week, tow, nav_sat):
    t_full = tow + week * 7 * 24 * 60 * 60
    toe_full = nav_sat[17] + nav_sat[27] * 7 * 24 * 60 * 60
    tk = t_full - toe_full
    # if tk != -2700: print("tk", tk)

    a = nav_sat[16] ** 2
    # if a != 26560761.5919: print("a", a)
    n0 = np.sqrt(mi/(a ** 3))
    # if n0 != 0.00014585057113006928: print("n0", n0)
    n = n0 + nav_sat[11]
    # if n != 0.00014585493166884608: print("n", n0)
    Mk = nav_sat[12] + n * tk
    # if Mk != -1.6609160957678846: print("Mk", Mk)

    E1 = Mk
    E2 = Mk + nav_sat[14] * np.sin(Mk)
    while abs(E2 - E1) > pow(10, -12):
        E1 = E2
        E2 = Mk + nav_sat[14] * m.sin(E2)
    Ek = E2

    vk = np.arctan2((np.sqrt(1 - nav_sat[14] ** 2) * np.sin(Ek)), np.cos(Ek) - nav_sat[14])
    # if vk != -1.6930032989900015: print("vk", vk)

    fik = vk + nav_sat[23]
    # if fik != -2.9199810375720014: print("fik", fik)

    delta_uk = nav_sat[15] * np.sin(2 * fik) + nav_sat[13] * np.cos(2 * fik)
    delta_rk = nav_sat[10] * np.sin(2 * fik) + nav_sat[22] * np.cos(2 * fik)
    delta_ik = nav_sat[20] * np.sin(2 * fik) + nav_sat[18] * np.cos(2 * fik)

    # if delta_uk == 4.908778137712827 * 10 ** -6: print("delta_uk", delta_uk)
    # if delta_rk != 311.26777: print("delta_rk", delta_rk)
    # if delta_ik != -2.3737800647519315 * 10 ** -7: print("delta_ik", delta_ik)

    uk = fik + delta_uk
    rk = a * (1 - nav_sat[14] * np.cos(Ek)) + delta_rk
    ik = nav_sat[21] + nav_sat[25] * tk + delta_ik

    # if uk != -2.9199761287938637: print("uk", uk)
    # if rk != 26606504.23982: print("rk", rk)
    # if ik != 0.9677757076332147: print("ik", ik)

    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)

    # if xk != -25955799.62960: print("xk", xk)
    # if yk != -5848293.20842: print("yk", yk)

    if np.abs(rk - np.sqrt(xk ** 2 + yk ** 2)) >= 0.01: print("bład")

    Omega_k = nav_sat[19] + (nav_sat[24] - omegaE) * tk - omegaE * nav_sat[17]
    # if Omega_k != -18.39607694561746: print("Omeka_k", Omega_k)

    Xk = xk * np.cos(Omega_k) - yk * np.cos(ik) * np.sin(Omega_k)
    Yk = xk * np.sin(Omega_k) + yk * np.cos(ik) * np.cos(Omega_k)
    Zk = yk * np.sin(ik)

    # if Xk != -21879348.443: print("Xk", Xk)
    # if Yk != -14352649.245: print("Yk", Yk)
    # if Zk != -4816807.994: print("Zk", Zk)

    #if np.abs(rk - np.sqrt(Xk ** 2 + Yk ** 2 + Zk ** 2)) >= 0.01:
    #    print("bład")

    delta_ts = nav_sat[6] + nav_sat[7] * (tow - nav_sat[17]) + nav_sat[8] * ((tow - nav_sat[17]) ** 2)
    delta_trel = (-2 * np.sqrt(mi) * nav_sat[14] * np.sqrt(a) * np.sin(Ek))/(c ** 2)
    delta_trels = delta_ts + delta_trel

    # if delta_ts != -0.0004752523245770398: print("delta_ts", delta_ts)
    # if delta_trels != 3.674986588631067 * 10 ** -8: print("delta_trel", delta_trel)
    # if delta_trels != -0.0004752155747111535: print("delta_trels", delta_trels)
    return ((Xk, Yk, Zk), delta_trels)

satpos(week, tow, nav_satelity_wybrana_epoka)
