import numpy as np
#SNS_iono na teams

alfa = [2.4214e-08, 0.0000e0, -1.19]
beta = [1.2902e5, ]

el = 30
az = 120
tgps = 43200
fi = 52
lam = 21

els = el / 180
azs = az / 180
fis = fi / 180
lams = lam / 180

psi = 0.0137 / (els + 0.11) - 0.022
fi_ipp = fis + psi * np.cos(np.deg2rad(az)

if fi_ipp > 0.416:
    fi_ipp = 0.416
elif fi_ipp < -0.416:
    fi_ipp = -0.416

lam_ipp = lams + psi * np.sin(np.deg2rad(az))/np.cos(fi_ipp * np.pi)

def klobuchar(tgps, el, az, fi, lam, alfa, beta):
