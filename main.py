import math as m
import numpy as np
from datetime import date

''' RINEX NAWIGACYJNY'''


def readrnxnav(file):
    m = 1
    nav = np.zeros((2000, 37))
    inav = np.zeros(2000)
    n = -1
    with open(file, "r") as f:
        for s in f:
            answer = s.find('END OF HEADER')  # skip header
            if answer != -1:
                break
        for s in f:
            s = s.replace('D', 'E')
            if m == 1:
                prn = int(s2n(s, 1, 2))
                a = np.empty((1, 6))
                a[:] = np.NaN
                a[0, 0:6] = np.array(s2e(s, 4, 23))
            else:
                a = np.append(a, s2n(s, 4, 19))
            for x in range(3):
                p = 23 + x * 19
                a = np.append(a, s2n(s, p, 19))
            if m < 8:
                m += 1
            else:
                n += 1
                nav[n, :] = a
                inav[n] = prn
                m = 1
        nav = nav[0:n + 1, :]
        inav = inav[0:n + 1]
        inav = inav.astype(int)
    f.close()
    return nav, inav


def readrnxobs(file, time_start, time_end, GNSS='G'):
    with open(file, "r") as f:
        for s in f:
            label = s[59:]
            if label.find('SYS / # / OBS TYPES') == 1:
                if s[0] == GNSS:
                    p = 7
                    types_header = []
                    for i in range(int(s[4:4 + 2])):
                        if p > 58:
                            p = 7
                            s = next(f)
                        types_header.append(s[p:p + 3])
                        p += 4

            elif label.find('APPROX POSITION XYZ') != -1:
                xr = np.array(([float(s[1:1 + 14]), float(s[15:15 + 14]), float(s[29:29 + 14])]))

            elif label.find('END OF HEADER') == 1:
                break
            types_of_obs = ['C1C']
        ind = np.zeros((len(types_header)))
        for n in range(len(types_of_obs)):
            i = (types_header.index(types_of_obs[n])) if types_of_obs[n] in types_header else -1  # np.empty((0))
            if i > -1:
                ind[i] = n + 1

        obs = np.zeros((150000, len(types_of_obs))) * np.nan
        iobs = np.zeros((150000, 3))
        n = 0
        for s in f:
            label = s[0]
            if label == '>':
                epoch = s2e(s, 2, 29)
                y = epoch[0]
                # tt = (date.toordinal(date(epoch[0],epoch[1],epoch[2]))+366-t0)*86400+np.dot((epoch[3:6]), ([3600,60,1])) + 6*86400
                tt = date2tow(epoch)[1] - date2tow(epoch)[2] * 86400
                if tt > (date2tow(time_end)[1] - date2tow(epoch)[2] * 86400):
                    break
                else:
                    flag = int(round(tt)) >= (date2tow(time_start)[1] - date2tow(time_start)[2] * 86400)
                if flag:
                    number_of_all_sats = int(s[33:33 + 2])
                    iobs[n + np.arange(0, number_of_all_sats), 1] = tt
                    iobs[n + np.arange(0, number_of_all_sats), 2] = date2tow(epoch)[1]
                    for sat in range(number_of_all_sats):
                        s = next(f)
                        p = 3
                        if s[0] == GNSS:
                            for i in range(len(types_header)):
                                if ind[i] != 0:
                                    obs[n + sat, int(ind[i] - 1)] = s2n(s, p, 16)
                                    iobs[n + sat, 0] = s2n(s, 1, 2)
                                p += 16
                    n += number_of_all_sats
        obs = obs[0:n, :]
        iobs = iobs[0:n, :]
        obs = np.delete(obs, iobs[:, 0] == 0, axis=0)
        iobs = np.delete(iobs, iobs[:, 0] == 0, axis=0)
        f.close()
        iobs = iobs.astype(int)
    return obs, iobs


def s2e(s, p, n):
    epoch = [int(s[p:p + 4]), int(s[p + 5:p + 5 + 2]), int(s[p + 8:p + 8 + 2]), int(s[p + 11:p + 11 + 2]),
             int(s[p + 14:p + 14 + 2]), float(s[p + 17:n])]
    return epoch


def date2tow(data):
    dday = date.toordinal(date(data[0], data[1], data[2])) - (date.toordinal(date(1980, 1, 6)))
    week = dday // 7
    dow = dday % 7
    tow = dow * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    return week, tow, dow

def date2tow2(data, rollover=False):
    dday = date.toordinal(date(data[0], data[1], data[2])) - (date.toordinal(date(1980, 1, 6)))
    week = dday // 7
    day = dday % 7
    tow = day * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    if rollover:
        if 0 < week <= 1023:
            week = week
        elif 1023 < week <= 2047:
            week = week - 2 ** 10
        elif 2047 < week <= 2 ** 12 - 1:
            week = week - 2 ** 11
    return week, tow

def s2n(s, p, n):
    a = s[p:p + n]
    if (not (a and not a.isspace())):
        a = np.nan
    else:
        a = float(a)
    return a

mi = 3.986005 * (10 ** 14)
omegaE = 7.2921151467 * (10 ** -5)
c = 299792458



nav_file = 'WROC00POL_R_20230750000_01D_GN.rnx'
nav, inav = readrnxnav(nav_file)
time_start = [2023, 3, 16, 0, 0, 0]
time_end = [2023, 3, 16, 23, 59, 59]
obs_file = 'WROC00POL_R_20230750000_01D_30S_MO.rnx'
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
#print(inav)
#print('--------------------')
#print(iobs)

sat = 1
ind_sat = inav==sat
nav_sat = nav[ind_sat]

week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

#week, t = date2tow2([2023, 3, 16, 0, 15, 0])
t_full = tow + week * 7 * 86400

xr, yr, zr = 3836000.0, 1177000.0, 4942000.0
xyzr = [xr,yr,zr]

toe_full = nav_sat[:, 17] + nav_sat[:, 27] * 7 * 86400
dt = t_full - toe_full
ind_dt = np.argmin(np.abs(dt))
#print(ind_dt)

nav_sat = nav_sat[ind_dt, :]
asqrt = nav_sat[16]

incl = nav_sat[21]

tk = dt[ind_dt]
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
    #print(xk, yk)

    if np.abs(rk - np.sqrt(xk ** 2 + yk ** 2)) >= 0.01:
        print("bład")

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



def satpos_received2sent_time(sat_pos : 'list[float]', tau):
    omegaE = 7.2921151467e-5
    tau = tau
    R = np.column_stack([[np.cos(omegaE*tau), -np.sin(omegaE*tau), 0],
                        [np.sin(omegaE*tau), np.cos(omegaE*tau), 0],
                        [0,0,1]])
    new_pos = R * np.array(sat_pos)
    return new_pos

a = 6378137
e2 = 0.00669438002290

def Np(B):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N


def hirvonen(X, Y, Z):
    r = (X ** 2 + Y ** 2) ** 0.5
    B = m.atan(Z / (r * (1 - e2)))

    while 1:
        N = Np(B)
        H = r / np.cos(B) - N
        Bst = B
        B = m.atan(Z / (r * (1 - (e2 * (N / (N + H))))))
        if abs(Bst - B) < (0.00001 / 206265):
            break
    L = m.atan(Y / X)
    N = Np(B)
    H = r / np.cos(B) - N
    return B, L, H

def wektor_odleglosci (xyz_s : 'np.array', xyz_r : 'np.array'):
    wektor_odleglosci = xyz_s - xyz_r
    return wektor_odleglosci

def RT(fi, lam):
    R = np.matrix([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                   [-np.sin(fi) * np.sin(lam), np.cos(lam), np.cos(fi) * np.sin(lam)],
                   [np.cos(fi), 0, np.sin(fi)]])
    RT = R.transpose()
    return RT


def xyz_chwilowe(wektor_odleglosci, RT):
    wektor = RT * wektor_odleglosci
    return wektor


def azymut(e, n):
    Az = np.arctan2(e, n)
    if(Az < 0):
        Az += 2 * np.pi
    return Az

def elewacja(xyz_chwilowe):
    el = np.arcsin(xyz_chwilowe[2] /
                   (np.sqrt(xyz_chwilowe[0] ** 2 + xyz_chwilowe[1] ** 2 + xyz_chwilowe[2] ** 2))) * 180 / np.pi
    return el

def odleglosc(wektor_odleglosci):
    s = np.sqrt(np.square(wektor_odleglosci[0]) + np.square(wektor_odleglosci[1]) + np.square(wektor_odleglosci[2]))
    return s


interval = 30
for tr in [345600]: # range(start_tow, end_tow + 1, interval):#gdzie interval równe 30
    ρs0 = []
    # indeksy, gdzie seskunda tygodnia jest równa naszej aktualnej sekundzie
    indt = iobs[:,2]==tr
    # dla tych indeksów satelity indt, zerowa kolumna
    sats = iobs[indt,0]

    tau = [0.07] * len(sats)
    # błąd zegara odbiornika
    dtr = 0

    for i in range(5):
        for sat_i, sat in enumerate(sats):
            ind_sat = inav==sat
            nav_sat = nav[ind_sat,:]
            # liczba sekund efemerydy dla wszystkich godzin, które są w efemerydzie (to jest lista długości 8)
            toe_full = nav_sat[:,17] + nav_sat[:,27] * 7 * 24 * 60 * 60
            # ile czsu upłynęło od toe_full do rozważanej epoki (to też jest lista długości 8)
            dt = tr - toe_full
            # wybór indeksu, który jest najbliżej godziny rozważanej, czyli ma najmniejszą różnicę dt. Innymi słowy wybór indeksu najmniejszej wartości z tablicy dt.
            ind_dt = np.argmin(abs(dt))
            # miejsce w efemerydzie dla tego satelity, które jest najbliżej naszego czasu.
            nav_sat_for_this_t = nav_sat[ind_dt,:]


            ts = tr + dtr - tau[sat_i]
            pos, dts = satpos(week, ts, nav_sat_for_this_t)
            new_pos = satpos_received2sent_time(pos, tau[sat_i])
            ρs0.append(np.linalg.norm(np.array(new_pos) - np.array(xyzr)))
            tau[sat_i] = ρs0[-1] / c
        #print(ρs0)
        #exit()


''''
delta_tr = 0
tau = 0.07 #zrób z tego tablicę [0.07, 0.07 ...] o długości liczby satelit
tr = 345600

#def receiver_position(wsp_przyblizone, nav_sat, ind_dt):
ts = tr + delta_tr - tau
x0s, y0s, z0s, delta_ts = satpos(week, tr, nav_sat)
x0, y0, z0 = 3836000.0, 1177000.0, 4942000.0
xyz = np.array(([np.cos(omegaE * tau), np.sin(omegaE * tau), 0],
                [-np.sin(omegaE * tau), np.cos(omegaE * tau), 0],
                [0, 0, 1]))
xyz_s = np.matmul(xyz, np.array([[x0s], [y0s], [z0s]]))
ro = np.sqrt((xyz_s[0] - x0)**2 + (xyz_s[1] - y0)**2 + (xyz_s[2] - z0)**2)
print("xyz_s: ", xyz_s, "ro: ", ro)

tau = ro/c
print("tau: ", tau)
fi = 51.114337987471046 * np.pi/180
lam = 17.057586085564054 * np.pi/180
h = 591.3087742133066

wektor_sat_odb = np.matrix([[x0s - x0], [y0s - y0], [z0s - z0]])
R = np.matrix([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                [-np.sin(fi) * np.sin(lam), np.cos(lam), np.cos(fi) * np.sin(lam)],
                [np.cos(fi), 0, np.sin(fi)]])
R = R.transpose()
Xsr_neu = R * wektor_sat_odb
print("Xsr_neu: ", Xsr_neu)
Az = np.arctan2(Xsr_neu[1], Xsr_neu[0]) * 180/np.pi
el = np.arcsin(Xsr_neu[2] / (np.sqrt(Xsr_neu[0] ** 2 + Xsr_neu[1] ** 2 + Xsr_neu[2] ** 2))) * 180/np.pi
print("Az: ", Az)
print("El: ", el)



for tr in range(coś):
    ind_t = iobs:[:,2] == tr
    pseudo.obs = obs[ind.t, 0]
    satelity = iobs[ind_t], 0]
        
for i in range(1): #isat enumerate ???
    for sat in satelity:
        i(nie mogę rozczytać)satnav = inav = sat
        efemerydy.1.sat = hov[i_satnav, :]
            
        ts = tr + dtr - tau[isat]
        [Xs, Ys, Zs] = f.satpos(efemerydy_1_sat, ts)
        tau[isat] = 

'''


def dTroposfera(Hr, el):
    p0 = 1013.25
    t0 = 291.15
    Rh0 = 0.5
    #Hel = 180.818
    #N = 40.231
    #Hort = Hel - N
    Hort = Hr
    p = p0 * (1 - (0.0000226 * Hort)) ** 5.225
    temp = t0 - (0.0065 * Hort)
    Rh = Rh0 * (np.e ** (-0.0006396 * Hort))
    e = 6.11 * Rh * 10 ** ((7.5 * (temp - 273.15)) /
                           temp - 35.85)
    c1 = 77.64
    c2 = -12.96
    c3 = 3.718 * (10 ** 5)
    Nd0 = c1 * p / temp
    Nw0 = (c2 * e / temp) + (c3 * e / (temp ** 2))
    hd = 40136 + 148.72 * (temp - 273.15)
    hw = 11000

    dTd0 = (10 ** -6) * Nd0 * hd / 5
    dTw0 = (10 ** -6) * Nw0 * hw / 5
    #mw = 1 / np.sin((el ** 2) + 6.25)
    #md = 1 / np.sin((el ** 2) + 2.25)
    mw = 1 / np.sin(np.radians((el ** 2) + 6.25))
    md = 1 / np.sin(np.radians((el ** 2) + 2.25))

    dTd = md * dTd0
    dTw = mw * dTw0
    dT = dTd + dTw

    return dT
for i in range(0, 100, 10):
    print("stopnie: " + str(i) + " || " + str(dTroposfera(140.587, i)))



'''
def dJonosfera(fi_r, lam_r)
    geocentr_angle = 0.0137/(elev + 0.11) - 0.22
    fiIPP = fi_r + geocentr_angle * np.cos(az)
    if fiIPP > 0.416:
        fiIPP = 0.4126
    if fiIPP < -0.416
        fiIPP = -0.4126
    lamIPP = lam_r + (geocentr_amgle * np.sin(az))/np.cos(lamIPP)
    fim = fiIPP + 0.064 * np.cos(lamIPP * np.pi - 1.617)
    t_local = 43200 * lamIPP + tGPS
    sod = fmod(sow, 86400)
'''


