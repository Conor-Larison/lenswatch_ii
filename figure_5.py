from ned_extinction_calc import request_extinctions
import sncosmo
import numpy as np
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

### Set up our SALT2-extended model

dust = sncosmo.F99Dust()
mod = sncosmo.Model('salt2-extended',effects=[dust],effect_names=['mw'],effect_frames=['obs'])
filts = ['Landolt B','Landolt V']

nra,ndec = 263.934701511,4.83245985403
A_lams = request_extinctions(nra, ndec, filters=filts)
ebv = A_lams[0] - A_lams[1]
mod.set(mwebv=ebv)

zwicky_data = ascii.read('fig5_data/sub_phot_ir.csv')

zwicky_data_mags = ascii.read('fig5_data/sub_phot_ir.csv')


zwicky_data['time'] = zwicky_data['mjd']
zwicky_data['band'] = zwicky_data['filter']

fluxerr1 = 10**(-0.4*(zwicky_data['mag']+zwicky_data['magerr']-zwicky_data['zp'])) - 10**(-0.4*(zwicky_data['mag']-zwicky_data['zp']))
fluxerr2 = 10**(-0.4*(zwicky_data['mag']-zwicky_data['magerr']-zwicky_data['zp'])) - 10**(-0.4*(zwicky_data['mag']-zwicky_data['zp']))
zwicky_data['fluxerr'] = (np.abs(fluxerr1) + np.abs(fluxerr2))/2.0
zwicky_data.remove_columns(['ra','dec','exp','mjd','filter','mag','magerr'])

zwicky_data['time'] = np.median(zwicky_data['time'])


zwicky_data = zwicky_data[zwicky_data.argsort('image')]

zwicky_data_a = zwicky_data[(zwicky_data['image'] == 'image_1')]
zwicky_data_b = zwicky_data[(zwicky_data['image'] == 'image_2')]
zwicky_data_c = zwicky_data[(zwicky_data['image'] == 'image_3')]
zwicky_data_d = zwicky_data[(zwicky_data['image'] == 'image_4')]

image_data = [zwicky_data_a,zwicky_data_b,zwicky_data_c,zwicky_data_d]


### We offset our t0 based on the measured time delays for each image

t0 = 59808.6
time_delays = [t0,t0+(0.52),t0+(1.97),t0+(0.78)]

z = 0.3544

images = 'ABCD'

fig,axs = plt.subplots(2,2,sharex=True)
fig.tight_layout(pad=2.5)

plt.subplots_adjust(hspace=0,wspace=0.15)


for i in range(len(image_data)):

    mod.set(t0=time_delays[i])
    mod.set(x1=1.16)
    mod.set(c=0.11)
    mod.set(z=z)


    ### Fit only the UVIS data

    res, fitted_model = sncosmo.nest_lc(image_data[i][np.where(image_data[i]['band']!='F160W')],mod,['x0'],bounds={},guess_amplitude_bound=True,modelcov=True)

    if i > 1:
        l = i - 2
        k = 1
    else:
        l = i
        k = 0

    formats = ['s','8','d','o']
    colors = ['orangered','slateblue','darkseagreen','yellowgreen']


    filter_ind = (image_data[i]['band'] == 'F160W')
    axs[k][l].errorbar((image_data[i]['time'][filter_ind]-t0)/(1+z),image_data[i]['flux'][filter_ind], \
            yerr=image_data[i]['fluxerr'][filter_ind],marker=formats[i],markersize=5,capsize=5, \
            elinewidth=3,zorder=1,color=colors[i])
    time_min = t0 - 8
    time_max = t0 + 40

    axs[k][l].text(0.9, 0.875, images[i],fontsize=12,fontweight='bold',
     horizontalalignment='center',
     verticalalignment='center',
     transform = axs[k][l].transAxes)

    times= np.linspace(time_min,time_max,1000)

    ### Get model flux errors

    f,mcov = fitted_model.bandfluxcov('f160w',times,24.663,'vega')
    mcov = 1.0857*np.sqrt(np.diagonal(mcov))/f

    mfluxerrs = mcov/(2.5/np.log(10))*fitted_model.bandflux('f160w',times,24.663,'vega')

    
    axs[k][l].fill_between((times-t0)/(1+z),fitted_model.bandflux('f160w',times,24.663,'vega')-mfluxerrs,fitted_model.bandflux('f160w',times,24.663,'vega')+mfluxerrs,alpha=0.4,color='gray')

    model_flux = fitted_model.bandflux('f160w',times,24.663,'vega')
    axs[k][l].plot((times-t0)/(1+z),model_flux,color='orangered')

    axs[k][l].tick_params(axis='both', which='major', labelsize=12)
    axs[k][l].tick_params(axis='both', which='minor', labelsize=10)

    mod_flux_ind = fitted_model.bandflux('f160w',image_data[i]['time'][0] + (time_delays[i]-t0),24.663,'vega')

    ftest,mcovtest = fitted_model.bandfluxcov('f160w',image_data[i]['time'][0]+(time_delays[i]-t0),24.663,'vega')

    mcovtest = 1.0857*np.sqrt(mcovtest)/ftest

    mod_flux_ind_err = mcovtest/(2.5/np.log(10))*fitted_model.bandflux('f160w',image_data[i]['time'][0]+(time_delays[i]-t0),24.663,'vega')

fig.text(0.0025, 0.5, r'flux ($ZP_{Vega}=24.663$)', fontsize=14, va='center',rotation='vertical')
fig.text(0.5, 0.025, 'phase (days)',fontsize=14, ha='center')


plt.show()
plt.close()