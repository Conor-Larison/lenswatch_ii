import sncosmo
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from matplotlib.ticker import MaxNLocator
from astropy.coordinates import SkyCoord
import corner
import pickle


### Function to convert the SN redshift to the CMB frame

def convert_redshift(ra,dec,redshift):
    sky = SkyCoord(ra,dec,unit=(u.hourangle, u.deg))

    l,b = sky.galactic.l.deg,sky.galactic.b.deg

    redshift_cmb = (1+redshift)/(1-369.82/3e5*(np.sin(b*np.pi/180)*np.sin(48.253*np.pi/180)+np.cos(b*np.pi/180)*np.cos(48.253*np.pi/180)*np.cos((l-264.021)*np.pi/180)))-1
    return redshift_cmb


### Set up the cosmology

cosmo = FlatLambdaCDM(H0 = 70 * u.km/u.s/u.Mpc,Tcmb0 = 2.725*u.K,Om0=0.3)


nra,ndec = 263.934701511,4.83245985403


### These are our standardization parameters

alpha = 0.14
beta = 3.1
x1 = 1.16
c = 0.11
MB = -19.36

z = 0.3544

z_cmb = convert_redshift(nra,ndec,z)


fig,axs = plt.subplots(2,2,figsize=(6,4),sharey=True,sharex=True)
fig.tight_layout(pad=2.5)

plt.subplots_adjust(hspace=0,wspace=0)


### The old measured magnifications and joint lens model predictions from P23

p23_mus = [8.31, 3.24, 6.73, 3.39]
p23_mus_pred = [1.81, 3.72, 2.87, 4.12]

images = 'ABCD'

mag_dict = {}

for i in range(4):

    ### Load in x0 values from our SALT2 fits

    sample_file = open('fig4_data/' + images[i] + '.pkl', 'rb')
    samples= np.array(pickle.load(sample_file))
    sample_file.close()

    x0 = np.array([])
    mB = np.array([])
    for l in range(len(samples)):
        x0 = np.append(x0,samples[l][0])
        mB = np.append(mB,-2.5*math.log10(x0[-1]) + 10.5)

    x0_errs = corner.quantile(x0,np.array([0.16,0.5,0.84]))


    ### We take our x0 posteriors and standardize using the above parameters, then compare
    ### to the cosmological prediction to get absolute magnitudes for the images

    mB_upper = -2.5*math.log10(x0_errs[2]) + 10.5
    mB_lower = -2.5*math.log10(x0_errs[0]) + 10.5

    mu_upper = mB_upper - MB + alpha*x1 - beta*c
    mu_lower = mB_lower - MB + alpha*x1 - beta*c

    mu_upper = np.sqrt(mu_upper**2 + 0.1**2)
    mu_lower = np.sqrt(mu_lower**2 - 0.1**2)

    mu = mB - MB + alpha*x1 - beta*c

    dist = cosmo.luminosity_distance(z_cmb)

    mu_pred = 5*np.log10(dist / u.Mpc) + 25

    mag_upper = 2.512**(np.abs(mu_upper - float(mu_pred)))
    mag_lower = 2.512**(np.abs(mu_lower - float(mu_pred)))


    magnifications = np.array(2.512**(np.abs(np.array(mu) - float(mu_pred))))
    
    bins=np.arange(-0,15, 2)
    image = images[i]

    mag_med = np.median(magnifications)

    mag_dict[image] = [mag_med - mag_lower,mag_med,mag_upper - mag_med]

    ### Now we can plot our absolute magnification posteriors with the predictions and old
    ### measurements overlaid

    if i > 1:
        l = i - 2
        k = 1
    else:
        l = i
        k = 0
    axs[k][l].hist(magnifications,density=True,bins=bins,histtype='step',color='k',linewidth=2)
    axs[k][l].set_xlim(0,15)
    axs[k][l].xaxis.set_major_locator(MaxNLocator(integer=True))
    axs[k][l].axvline(p23_mus[i],color='slateblue',linewidth=2,linestyle='dashed',label=r'P23 $\mu_{meas}$')
    axs[k][l].axvline(p23_mus_pred[i],color='orangered',linewidth=2,linestyle='dashdot',label=r'P23 $\mu_{pred}$')
    axs[k][l].text(0.3,0.32,image,fontsize=15)
    axs[k][l].xaxis.set_ticks(np.arange(1, 14, 2))

    axs[k][l].tick_params(axis='both', which='major', labelsize=12)
    axs[k][l].tick_params(axis='both', which='minor', labelsize=10)


fig.text(0.5, 0.025, r'absolute magnification (|$\mu$|)',fontsize=14, ha='center')
fig.text(0.01, 0.5, 'posterior (density)',fontsize=14, va='center',rotation='vertical')
plt.legend(fontsize=11)

plt.show()
plt.close()