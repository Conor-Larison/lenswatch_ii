import matplotlib.pyplot as plt
import numpy as np
import pickle
import corner

### Read in our SNTD results

sample_file = open('fig3_data/sntd_samples.pkl', 'rb')
samples= np.array(pickle.load(sample_file))
sample_file.close()

weight_file = open('fig3_data/sntd_weights.pkl', 'rb')
weights= np.array(pickle.load(weight_file))
weight_file.close()


# Corner plot!! :)
    
fig = corner.corner(
    samples,
    weights=weights,
    labels=[r'$t_{pk}$','SALT2 c',r'$\Delta t_{BA}$',r'$\Delta t_{CA}$',r'$\Delta t_{DA}$'],
    quantiles=(0.16, .5, 0.84),
    bins=30,
    color='k',
    show_titles=True,
    title_fmt='.2f',
    smooth1d=False,
    smooth=True,
    fill_contours=False,
    plot_contours=True,
    plot_density=True,
    use_mathtext=True,
    title_kwargs={"fontsize": 15},
    label_kwargs={'fontsize': 16})
for ax in fig.get_axes():
    ax.tick_params(axis='both', labelsize=14)


plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plt.show()
plt.close()