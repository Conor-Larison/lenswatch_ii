import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from photutils.background import LocalBackground,MMMBackground
import space_phot
from scipy.odr import Model,RealData,ODR


### Read in our photometry values and drizzled images


old_phot_file = 'fig2_data/lenswatch_phot_mags.csv'
new_phot_file = 'fig2_data/sub_phot_final.csv'


temp_images_drz = glob.glob('fig2_data/temp_files/*drz.fits')


sn_location = SkyCoord(263.9347929,4.8324551,unit=u.degree)

lens_flux = {'F814W':[0,0],'F625W':[0,0],'F475W':[0,0]}

temp_filts = ['F475W','F625W','F814W']


for i,file in enumerate(temp_images_drz):
    ### Here, we calculate the aperture photometry for the lens galaxy in our UVIS filters,
    ### along with their errors.

    data = fits.open(file)['SCI',1].data

    wcs = WCS(fits.open(file)['SCI',1].header)


    x,y = wcs.world_to_pixel(sn_location)

    y,x = int(y),int(x)+2

    sky = SkyCoord(263.934701511,4.83245985403,unit=u.deg)

    hst_obs = space_phot.observation3(file)
    hst_obs.aperture_photometry(sky,xy_positions=[x,y],radius=14,
                        skyan_in=14,skyan_out=20)


    lflux= hst_obs.aperture_result.phot_cal_table['flux']

    lens_flux[temp_filts[i]][0] = lflux

    aperture_errs = hst_obs.aperture_result.phot_cal_table['fluxerr']

    lens_flux[temp_filts[i]][1] = aperture_errs


old_data = ascii.read(old_phot_file)
new_data = ascii.read(new_phot_file)

zps = np.unique(new_data['zp'])

### I ended up hard-coding in the old photometry

old_phot_dict = {'F814W_1': [20.67,0.01,zps[0]],'F625W_1': [21.67,0.02,zps[1]],
                  'F475W_1': [23.22,0.04,zps[2]],
                  'F814W_2': [21.71,0.02,zps[0]],
                  'F625W_2': [22.65,0.03,zps[1]],'F475W_2': [24.31,0.07,zps[2]],
                  'F814W_3': [20.88,0.02,zps[0]],
                  'F625W_3': [21.90,0.02,zps[1]],'F475W_3': [23.35,0.04,zps[2]],
                  'F814W_4': [21.60,0.02,zps[0]],
                  'F625W_4': [22.72,0.04,zps[1]],'F475W_4': [24.26,0.07,zps[2]]}

filters = ['UVF814W','UVF625W','UVF475W']
images = ['image_1','image_2','image_3','image_4']

flux_differences = np.array([])
flux_diff_errs = np.array([])

norm_lens_fluxes = np.array([])
norm_lens_errs = np.array([])
for i,image in enumerate(images):
    for j,filt in enumerate(filters):

        ### Here we calculate the flux ratios and propagate their uncertainties

        new_flux = new_data[np.where((new_data['filter'] == filt) & (new_data['image'] == image))]['flux']
        new_flux_err = new_data[np.where((new_data['filter'] == filt) & (new_data['image'] == image))]['fluxerr']

        key_start = filt[2:]

        lflux = lens_flux[key_start][0]
        lflux_err = lens_flux[key_start][1]

        key = key_start + image[-2:]

        old_flux = 10**((old_phot_dict[key][0] - old_phot_dict[key][2])/-2.5)
        old_flux_err = old_flux - 10**((old_phot_dict[key][0] + old_phot_dict[key][1] - old_phot_dict[key][2])/-2.5)

        flux_differences = np.append(flux_differences,new_flux/old_flux)
        flux_diff_err = np.sqrt((1/old_flux)**2*(new_flux_err)**2 + (new_flux/old_flux**2)**2*(old_flux_err)**2)

        flux_diff_errs = np.append(flux_diff_errs,flux_diff_err)

        norm_lens_fluxes = np.append(norm_lens_fluxes,old_flux/lflux)
        norm_lens_err = np.sqrt((1/lflux)**2*(new_flux_err)**2 + (old_flux/lflux**2)**2*(lflux_err)**2)

        norm_lens_errs = np.append(norm_lens_errs,norm_lens_err)


def linear(B, x):
    return B[0]*x + B[1]

# Fit a linear model using scipy ODR
lin_model = Model(linear)

# Create a RealData object using our initiated data from above
data = RealData(norm_lens_fluxes, flux_differences, sx=norm_lens_errs, sy=flux_diff_errs)

# Set up ODR with the model and data
odr = ODR(data, lin_model, beta0=[-1.,1.5])

# Run the regression.
out = odr.run()

x_imp = (1.0-out.beta[1])/out.beta[0]
x_imp_err = np.sqrt((1/out.beta[0])**2*(out.sd_beta[1])**2 + ((1.0-out.beta[1])/out.beta[0]**2)**2*(out.sd_beta[0])**2)


x_fit = np.linspace(0.2, 1.3, 1000)
y_fit = linear(out.beta, x_fit)

fig = plt.figure()

plt.plot(x_fit, y_fit,linewidth=2,color='k')

plt.hlines(1.0,0.2,1.3,linestyle='dashed',linewidth=2,color='red')

colors = ['orangered','slateblue','darkseagreen','yellowgreen']
formats = ['s','8','d','o']


plt.errorbar(norm_lens_fluxes[0:3],flux_differences[0:3],xerr=norm_lens_errs[0:3],yerr=flux_diff_errs[0:3],label='Image A',fmt=formats[0],color=colors[0],markersize=8)
plt.errorbar(norm_lens_fluxes[3:6],flux_differences[3:6],xerr=norm_lens_errs[3:6],yerr=flux_diff_errs[3:6],label='Image B',fmt=formats[1],color=colors[1],markersize=8)
plt.errorbar(norm_lens_fluxes[6:9],flux_differences[6:9],xerr=norm_lens_errs[6:9],yerr=flux_diff_errs[6:9],label='Image C',fmt=formats[2],color=colors[2],markersize=8)
plt.errorbar(norm_lens_fluxes[9:],flux_differences[9:],xerr=norm_lens_errs[9:],yerr=flux_diff_errs[9:],label='Image D',fmt=formats[3],color=colors[3],markersize=8) 

plt.ylabel(r'$\mathrm{F_{new}/F_{old}}$',fontsize=15)
plt.xlabel(r'$\mathrm{F_{old}/F_{contam}}$',fontsize=15)
plt.xlim(0.2,1.3)
plt.legend(fontsize=15)
plt.tick_params(labelsize=13)
plt.show()
plt.close()
