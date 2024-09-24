import matplotlib.pyplot as plt
import numpy as np

### Load in the posteriors for the different f_subs

name = 'ls1'

def load_mags(f_sub):
    magsA = np.loadtxt('fig7_data/'+name+'_magA_fsub'+str(f_sub)+'.txt')
    magsB = np.loadtxt('fig7_data/'+name+'_magB_fsub'+str(f_sub)+'.txt')
    magsC = np.loadtxt('fig7_data/'+name+'_magC_fsub'+str(f_sub)+'.txt')
    magsD = np.loadtxt('fig7_data/'+name+'_magD_fsub'+str(f_sub)+'.txt')
    return magsA, magsB, magsC, magsD

f_sub = 1
magsA1, magsB1, magsC1, magsD1 = load_mags(f_sub)
f_sub = 5
magsA2, magsB2, magsC2, magsD2 = load_mags(f_sub)
f_sub = 10
magsA3, magsB3, magsC3, magsD3 = load_mags(f_sub)

### Smooth model values are hard coded in, but can be calculated from the convergences and shears
### from P23

mag_image_A,mag_image_B,mag_image_C,mag_image_D = 2.4630541871921188,4.424778761061948,4.0,4.854368932038836

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.width'] = 3.5
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.minor.size'] = 6
plt.rcParams['ytick.minor.size'] = 6
plt.rcParams['ytick.major.width'] = 3.5
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 18

bins = 50
fig,axs = plt.subplots(2,2,layout='tight')
fig.set_size_inches(13.5, 9)


bins = np.arange(mag_image_A-0.1,mag_image_A+0.1,0.01)

ran=(mag_image_A*0.9, mag_image_A*1.1)
axs[0][0].hist(magsA1, bins=bins, range=ran,alpha=0.5,color='k',density=True, label=r'$f_{\rm{sub}} = 1 \%$')
axs[0][0].hist(magsA2, bins=bins, range=ran,alpha=0.5,color='r',density=True, label=r'$f_{\rm{sub}} = 5 \%$')
axs[0][0].hist(magsA3, bins=bins, range=ran,alpha=0.5,color='b',density=True, label=r'$f_{\rm{sub}} = 10 \%$')
axs[0][0].axvline(np.round(mag_image_A,2), color='k', label=r'$\rm{smooth \ model}$', lw=5, linestyle='--')
axs[0][0].legend(fontsize=15, frameon=False,loc=2)
axs[0][0].annotate('A', xy=(0.92, 0.88), xycoords='axes fraction', fontsize=24, fontweight='bold')


axs[0][0].set_yticks(np.round(np.linspace(axs[0][0].get_ylim()[0],axs[0][0].get_ylim()[1],5),1))
axs[0][0].set_xlim(mag_image_A-0.1,mag_image_A+0.1)
axs[0][0].set_xticks(np.round(np.linspace(mag_image_A-0.12,mag_image_A+0.12,5),2))


bins = np.arange(mag_image_B-0.2,mag_image_B+0.2,0.02)

ran=(mag_image_B*0.9, mag_image_B*1.1)
axs[0][1].hist(magsB1, bins=bins, range=ran,alpha=0.5,color='k',density=True, label=r'$f_{\rm{sub}} = 1 \%$')
axs[0][1].hist(magsB2, bins=bins, range=ran,alpha=0.5,color='r',density=True, label=r'$f_{\rm{sub}} = 5 \%$')
axs[0][1].hist(magsB3, bins=bins, range=ran,alpha=0.5,color='b',density=True, label=r'$f_{\rm{sub}} = 10 \%$')
axs[0][1].axvline(np.round(mag_image_B,2), color='k', label=r'$\rm{smooth \ model}$', lw=5, linestyle='--')
axs[0][1].annotate('B', xy=(0.92, 0.88), xycoords='axes fraction', fontsize=24, fontweight='bold')


axs[0][1].set_yticks(np.round(np.linspace(axs[0][1].get_ylim()[0],axs[0][1].get_ylim()[1],5),1))
axs[0][1].set_xlim(mag_image_B-0.2,mag_image_B+0.2)
axs[0][1].set_xticks(np.round(np.linspace(mag_image_B-0.2,mag_image_B+0.2,5),2))


bins = np.arange(mag_image_C-0.15,mag_image_C+0.15,0.015)

ran=(mag_image_C*0.9, mag_image_C*1.1)
axs[1][0].hist(magsC1, bins=bins, range=ran,alpha=0.5,color='k',density=True, label=r'$f_{\rm{sub}} = 1 \%$')
axs[1][0].hist(magsC2, bins=bins, range=ran,alpha=0.5,color='r',density=True, label=r'$f_{\rm{sub}} = 5 \%$')
axs[1][0].hist(magsC3, bins=bins, range=ran,alpha=0.5,color='b',density=True, label=r'$f_{\rm{sub}} = 10 \%$')
axs[1][0].axvline(mag_image_C, color='k', label=r'$\rm{smooth \ model}$', lw=5, linestyle='--')
axs[1][0].annotate('C', xy=(0.92, 0.88), xycoords='axes fraction', fontsize=24, fontweight='bold')


axs[1][0].set_yticks(np.round(np.linspace(axs[1][0].get_ylim()[0],axs[1][0].get_ylim()[1],5),1))
axs[1][0].set_xlim(mag_image_C-0.15,mag_image_C+0.15)
axs[1][0].set_xticks(np.round(np.linspace(mag_image_C-0.15,mag_image_C+0.15,5),2))


bins = np.arange(mag_image_D-0.2,mag_image_D+0.2,0.02)


ran=(mag_image_D*0.9, mag_image_D*1.1)
axs[1][1].hist(magsD1, bins=bins, range=ran,alpha=0.5,color='k',density=True, label=r'$f_{\rm{sub}} = 1 \%$')
axs[1][1].hist(magsD2, bins=bins, range=ran,alpha=0.5,color='r',density=True, label=r'$f_{\rm{sub}} = 5 \%$')
axs[1][1].hist(magsD3, bins=bins, range=ran,alpha=0.5,color='b',density=True, label=r'$f_{\rm{sub}} = 10 \%$')
axs[1][1].axvline(np.round(mag_image_D,2), color='k', label=r'$\rm{smooth \ model}$', lw=5, linestyle='--')
axs[1][1].annotate('D', xy=(0.92, 0.88), xycoords='axes fraction', fontsize=24, fontweight='bold')


axs[1][1].set_xlim(mag_image_D-0.2,mag_image_D+0.2)
axs[1][1].set_yticks(np.round(np.linspace(axs[1][1].get_ylim()[0],axs[1][1].get_ylim()[1],5),1),labelsize=20)
axs[1][1].set_xticks(np.round(np.linspace(mag_image_D-0.2,mag_image_D+0.2,5),2),labelsize=20)


fig.text(0.5, 0.005, r'absolute magnification (|$\mu$|)',fontsize=23, ha='center')
fig.text(0.005, 0.5, 'pyHalo posterior',fontsize=23, va='center',rotation='vertical')

fig.tight_layout(pad=3)
plt.subplots_adjust(hspace=0.2,wspace=0.2)


plt.show()
plt.close()