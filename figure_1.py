from astropy.io import fits
import space_phot
import matplotlib.pyplot as plt
import astropy
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
import montage_wrapper as montage


### HST files used in Figure 1

orig_files = ['fig1_data/orig_files/iebc16u6q_flc.fits','fig1_data/orig_files/iebc16uaq_flc.fits','fig1_data/orig_files/iebc16ubq_flc.fits']

temp_files = ['fig1_data/temp_files/iebc17eiq_flc.fits','fig1_data/temp_files/iebc17emq_flc.fits','fig1_data/temp_files/iebc17enq_flc.fits']

sub_files = ['fig1_data/sub_files/iebc16u6q_flc.fits','fig1_data/sub_files/iebc16uaq_flc.fits','fig1_data/sub_files/iebc16ubq_flc.fits']

psf_sub_files = ['fig1_data/psf_sub_files/F475W_residuals.fits','fig1_data/psf_sub_files/F625W_residuals.fits','fig1_data/psf_sub_files/F814W_residuals.fits']


## We separate NIR data because the scaling tends to be different

f160_orig_file = 'fig1_data/orig_files/iebc16udq_flt.fits'
f160_temp_file = 'fig1_data/temp_files/iebc17epq_flt.fits'
f160_sub_file = 'fig1_data/sub_files/iebc16udq_flt_astro_replaced.fits'
f160_psf_file = 'fig1_data/psf_sub_files/F160W_residuals.fits'


filters = ['F475W','F625W','F814W']


fig,axs = plt.subplots(4,4,layout='tight',sharex=False,sharey=False,figsize=(7.5,7.5),subplot_kw=dict(projection=WCS(fits.open(orig_files[0])[0].header)))


### These centers were hard-coded from examining the files, thus,
### the alignments between filters may not be perfect for this figure


centers = {'F814W':[(345,370)],'F625W':[(345,369)],
           'F475W':[(345,366)]}
centers_o = {'F814W':[(centers['F814W'][0][0]+2,centers['F814W'][0][1]+3)],
            'F625W':[(centers['F625W'][0][0]+2,centers['F625W'][0][1]+3)],
            'F475W':[(centers['F475W'][0][0]+2,centers['F475W'][0][1]+3)]}

centers_temp = {'F814W':[(centers['F814W'][0][0]+3,centers['F814W'][0][1])],
            'F625W':[(centers['F625W'][0][0]+3,centers['F625W'][0][1])],
            'F475W':[(centers['F475W'][0][0]+3,centers['F475W'][0][1])]}

ir_centers = {'F160W':[(381,401)]}
ir_centers_o = {'F160W':[(ir_centers['F160W'][0][0],ir_centers['F160W'][0][1]+1)]}


for i, filt in enumerate(filters):


    ### We read in the files and reproject them, so North is up and East is left
    
    orig_hdu = montage.reproject_hdu(fits.open(orig_files[i])[2], north_aligned=True)
    orig_data = orig_hdu.data

    temp_hdu = fits.open(orig_files[i])[2]
    temp_hdu.data = fits.open(temp_files[i])['SCI',1].data

    temp_data = montage.reproject_hdu(temp_hdu, north_aligned=True).data


    sub_hdu = fits.open(orig_files[i])[2]
    sub_hdu.data = fits.open(sub_files[i])['SCI',1].data

    sub_data = montage.reproject_hdu(sub_hdu, north_aligned=True).data


    psf_hdu = fits.open(orig_files[i])[2]
    psf_hdu.data = fits.open(psf_sub_files[i])['SCI',1].data

    psf_data = montage.reproject_hdu(psf_hdu, north_aligned=True).data


    for l in range(4):
        axs[i][l].tick_params(
        axis='both',
        length=0, 
        width=0,     
        labelsize=0
        )


    ### Norms may vary between UVIS filters and F160W filter panels

    norm1 = astropy.visualization.simple_norm(orig_data[centers_o[filt][0][0]-10:centers_o[filt][0][0]+11,centers_o[filt][0][1]-10:centers_o[filt][0][1]+11],stretch='log',invalid=0)
    norm2 = astropy.visualization.simple_norm(sub_data[centers_o[filt][0][0]-5:centers_o[filt][0][0]+6,centers_o[filt][0][1]-5:centers_o[filt][0][1]+6],stretch='log',invalid=0)


    axs[i][0].imshow(orig_data[centers_o[filt][0][0]-13:centers_o[filt][0][0]+14,centers_o[filt][0][1]-13:centers_o[filt][0][1]+14],norm=norm1,cmap='Greys',origin='lower')
    axs[i][1].imshow(temp_data[centers_temp[filt][0][0]-13:centers_temp[filt][0][0]+14,centers_temp[filt][0][1]-13:centers_temp[filt][0][1]+14],norm=norm2,cmap='Greys')
    axs[i][2].imshow(sub_data[centers_o[filt][0][0]-13:centers_o[filt][0][0]+14,centers_o[filt][0][1]-13:centers_o[filt][0][1]+14],norm=norm2,cmap='Greys')
    axs[i][3].imshow(psf_data[centers_o[filt][0][0]-13:centers_o[filt][0][0]+14,centers_o[filt][0][1]-13:centers_o[filt][0][1]+14],norm=norm2,cmap='Greys')
    
    axs[i][0].set_ylabel(filt,fontsize=14)
    
    axs[0][0].set_title('Original',fontsize=14)
    axs[0][1].set_title('Template',fontsize=14)
    axs[0][2].set_title('Template Sub',fontsize=14)
    axs[0][3].set_title('PSF Sub',fontsize=14)


### Now we do the same for the F160W images

f160_orig_hdu = montage.reproject_hdu(fits.open(f160_orig_file)[2], north_aligned=True)
f160_orig_data = f160_orig_hdu.data

f160_temp_hdu = fits.open(f160_orig_file)[2]
f160_temp_hdu.data = fits.open(f160_temp_file)['SCI',1].data

f160_temp_data = montage.reproject_hdu(f160_temp_hdu, north_aligned=True).data

f160_sub_hdu = fits.open(f160_orig_file)[2]
f160_sub_hdu.data = fits.open(f160_sub_file)['SCI',1].data

f160_sub_data = montage.reproject_hdu(f160_sub_hdu, north_aligned=True).data

f160_psf_hdu = fits.open(f160_orig_file)[2]
f160_psf_hdu.data = fits.open(f160_psf_file)['SCI',1].data

f160_psf_data = montage.reproject_hdu(f160_psf_hdu, north_aligned=True).data

x_width,y_width = 20,20

norm1 = astropy.visualization.simple_norm(f160_orig_data[ir_centers_o['F160W'][0][0]-x_width:ir_centers_o['F160W'][0][0]+x_width+1,ir_centers_o['F160W'][0][1]-y_width:ir_centers_o['F160W'][0][1]+y_width+1],stretch='log',invalid=0)
z1, z2 = ZScaleInterval().get_limits(f160_sub_data[ir_centers['F160W'][0][0]-x_width:ir_centers['F160W'][0][0]+x_width + 1,ir_centers['F160W'][0][1]-y_width:ir_centers['F160W'][0][1]+y_width+1])

lower_width,upper_width = [4,4]


axs[3][0].imshow(f160_orig_data[ir_centers_o['F160W'][0][0]-lower_width:ir_centers_o['F160W'][0][0]+upper_width + 1,ir_centers_o['F160W'][0][1]-lower_width:ir_centers_o['F160W'][0][1]+upper_width+1],vmin=-1.5,cmap='Greys')
axs[3][1].imshow(f160_temp_data[ir_centers['F160W'][0][0]-lower_width:ir_centers['F160W'][0][0]+upper_width + 1,ir_centers['F160W'][0][1]-lower_width:ir_centers['F160W'][0][1]+upper_width+1],norm=norm1,cmap='Greys')
axs[3][2].imshow(f160_sub_data[ir_centers['F160W'][0][0]-lower_width:ir_centers['F160W'][0][0]+upper_width + 1,ir_centers['F160W'][0][1]-lower_width:ir_centers['F160W'][0][1]+upper_width+1],vmin=z1-100,vmax=z2+1000,cmap='Greys')
axs[3][3].imshow(f160_psf_data[ir_centers['F160W'][0][0]-lower_width:ir_centers['F160W'][0][0]+upper_width + 1,ir_centers['F160W'][0][1]-lower_width:ir_centers['F160W'][0][1]+upper_width+1],vmin=z1-100,vmax=z2+1000,cmap='Greys')


for l in range(4):
    axs[3][l].tick_params(
        axis='both',
        length=0, 
        width=0,     
        labelsize=0
        )

axs[3][0].set_ylabel('F160W',fontsize=14)

axs[0][2].text(16,19,r'A$^{\mathbf{-}}$',color='k',fontsize=15,fontweight='heavy')
axs[0][2].text(20,7,r'B$^{\mathbf{+}}$',color='k',fontsize=15,fontweight='heavy')
axs[0][2].text(8.0,4.5,r'C$^{\mathbf{-}}$',color='k',fontsize=15,fontweight='heavy')
axs[0][2].text(5.0,15,r'D$^{\mathbf{+}}$',color='k',fontsize=15,fontweight='heavy')


axs[0][2].arrow(13,1,13, 0, 
         head_width=0, head_length=0, 
         fc='k', ec='k', width=0.5)
axs[0][2].text(17,2, '0.5"', fontweight='heavy',
         color='k')

    
plt.show()
plt.close()