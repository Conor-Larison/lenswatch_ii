from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import space_phot
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from photutils.aperture import ApertureStats
from photutils.background import LocalBackground,MMMBackground


### Read in the DOLPHOT photometry

phot = Table.read('fig6_data/22qmx_F160W.phot', format='ascii.no_header')


qmet_idx = ((phot['col16']<99)& # VEGMAG in F475W
            (phot['col23']<0.4)& # Crowding
            (phot['col21']**2<0.1)& # Sharpness^2
            (phot['col11']<=2) # Objtype
           )

qmet_spatial = (qmet_idx) & (phot['col3'] > 255) & (phot['col3'] < 270) \
               & (phot['col4'] > 260) & (phot['col4'] < 275)

og_data1 = fits.open('fig6_data/iebc16ucq_flt_astro_replaced.fits')['SCI',1].data
og_data2 = fits.open('fig6_data/iebc16udq_flt_astro_replaced.fits')['SCI',1].data
og_data3 = fits.open('fig6_data/iebc16ueq_flt_astro_replaced.fits')['SCI',1].data

files = ['fig6_data/iebc16ucq_flt_astro_replaced.fits','fig6_data/iebc16udq_flt_astro_replaced.fits','fig6_data/iebc16ueq_flt_astro_replaced.fits']


positions = [[264,262],[267,266],[270,269]]
apertures = [CircularAperture(positions[i], r=6.0) for i in range(len(positions))]


localbkg_estimator = LocalBackground(7, 15, bkg_estimator=MMMBackground())
bkgs = [localbkg_estimator(og_data1,positions[i][0],positions[i][1]) for i in range(len(positions))]


### Get aperture photometry of the SN images and the associated errors

sky = SkyCoord(263.934701511,4.83245985403,unit=u.deg)

hst_obs = space_phot.observation2(files)
hst_obs.aperture_photometry(sky,radius=6,
                    skyan_in=7,skyan_out=15)

aperture_errs = hst_obs.aperture_result.phot_cal_table['flux_cal_err']


phot_table1 = aperture_photometry(og_data1 - bkgs[0], apertures[0])
phot_table2 = aperture_photometry(og_data2 - bkgs[1], apertures[1])
phot_table3 = aperture_photometry(og_data3 - bkgs[2], apertures[2])

test = ApertureStats(og_data1 - bkgs[0], apertures[0])


phot_table1['aperture_sum'].info.format = '%.8g'  # for consistent table output
phot_table2['aperture_sum'].info.format = '%.8g'  # for consistent table output
phot_table3['aperture_sum'].info.format = '%.8g'  # for consistent table output


ftot = np.mean([phot_table1['aperture_sum'][0],phot_table2['aperture_sum'][0],phot_table3['aperture_sum'][0]])
ftot_err = np.mean(aperture_errs)
fluxes = 10**((phot[qmet_spatial]['col29'] - 24.662)/-2.5)
fluxerrs = phot[qmet_spatial]['col31']/(2.5/np.log(10))*fluxes

### Image B had index 2, so we switch it to be in order


fluxes[[1,2]] = fluxes[[2,1]]
fluxerrs[[1,2]] = fluxerrs[[2,1]]


mu_pred = np.array([1.81,3.72,2.87,4.12])
mu_pred_errs_upper = np.array([0.9,1.04,1.51,1.19])
mu_pred_errs_lower = np.array([0.89,1.24,1.50,1.36])


### Next, we make our predictions for the image fluxes in F160W and propagate the errors

fluxes_pred = [ftot/np.sum(mu_pred)*mu_pred[i] for i in range(len(mu_pred))]


ind_flux = ftot/np.sum(mu_pred)

ind_flux_err_upper = np.sqrt((1/np.sum(mu_pred))**2*(ftot_err)**2 + np.sum((ftot/(np.sum(mu_pred))**2)**2*(mu_pred_errs_upper)**2))
ind_flux_err_lower = np.sqrt((1/np.sum(mu_pred))**2*(ftot_err)**2 + np.sum((ftot/(np.sum(mu_pred))**2)**2*(mu_pred_errs_lower)**2))


fluxes_pred_upper = [fluxes_pred[i] + np.sqrt((mu_pred[i])**2*(ind_flux_err_upper)**2 + (ind_flux)**2*(mu_pred_errs_upper[i])**2) for i in range(len(mu_pred))]
fluxes_pred_lower = [fluxes_pred[i] - np.sqrt((mu_pred[i])**2*(ind_flux_err_lower)**2 + (ind_flux)**2*(mu_pred_errs_lower[i])**2) for i in range(len(mu_pred))]

formats = ['s','8','d','o']
colors = ['orangered','slateblue','darkseagreen','yellowgreen']
labels = 'ABCD'


fluxes_pred = np.array(fluxes_pred)

fluxerrstest = np.sqrt(fluxerrs**2)

fig,axs = plt.subplots(1)

hatches = ['//', '\\\\', '||', '--']
linestyles = ['dashed','dashdot','solid','dotted']
formats = ['s','8','d','o']

pred_errs_upper = fluxes_pred_upper - fluxes_pred
pred_errs_lower = np.abs(fluxes_pred_lower - fluxes_pred)


[axs.vlines(fluxes[i],ymin=-5,ymax=5,color=colors[i],label=labels[i],linestyle=linestyles[i],linewidth=2) for i in range(len(fluxes_pred))]


[axs.plot(np.linspace(fluxes_pred[i],fluxes_pred[i] + 8*pred_errs_upper[i],1000), 3 *
        np.exp( - (np.linspace(fluxes_pred[i],fluxes_pred[i] + 8*pred_errs_upper[i],1000) - fluxes_pred[i])**2 / (2 * pred_errs_upper[i]**2) ),
         linewidth=2, color=colors[i],linestyle=linestyles[i]) for i in range(len(fluxes_pred))]

[axs.plot(np.linspace(fluxes_pred[i] - 8*pred_errs_lower[i],fluxes_pred[i],1000), 3 *
        np.exp( - (np.linspace(fluxes_pred[i] - 8*pred_errs_lower[i],fluxes_pred[i],1000) - fluxes_pred[i])**2 / (2 * pred_errs_lower[i]**2) ),
         linewidth=2, color=colors[i],linestyle=linestyles[i]) for i in range(len(fluxes_pred))]


axs.set_ylim(0.0,4)
[axs.axvspan(fluxes[i] - fluxerrstest[i], fluxes[i] + fluxerrstest[i], alpha=0.4, color=colors[i],hatch=hatches[i]) for i in range(len(colors))]

axs.set_xlabel('flux (counts / s)',fontsize=14)
axs.set_ylabel('Prediction PDF',fontsize=14)
axs.set_yticks([])
axs.set_xlim(8,50)
axs.tick_params(labelsize=12)
plt.title('F160W Predicted Fluxes vs. Measured')
plt.legend(loc=2,fontsize=11)

plt.show()
plt.close()