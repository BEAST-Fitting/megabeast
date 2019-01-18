# system
from __future__ import (absolute_import, division, print_function)

# other packages
from tqdm import (tqdm, trange)
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

# beast
from beast.physicsmodel.prior_weights_dust import _lognorm, _two_lognorm

# megabeast
from ..read_megabeast_input import read_megabeast_input
from ..beast_data import (read_beast_data, read_lnp_data)



def simulate_av_plots(megabeast_input_file, log_scale=False,
                          input_lognormal=None, input_lognormal2=None):
    """
    Plot distributions of simulated AVs, and overplot the best fit lognormals

    Parameters
    ----------
    megabeast_input_file : string
        Name of the file that contains settings, filenames, etc

    log_scale : boolean (default=False)
        If True, make the histogram x-axis a log scale (to visualize log-normal
        A_V distribution)

    input_lognormal, input_lognormal2 : dict (default=None)
        Set these to the original values used to create the fake data, and they
        will also be plotted

    """

    # read in the settings from the file
    mb_settings = read_megabeast_input(megabeast_input_file)

    # get the project name
    projectname = mb_settings['projectname']

    # read in the beast data that is needed by all the pixels
    beast_data = read_beast_data(mb_settings['beast_seds_filename'],
                                 mb_settings['beast_noise_filename'],
                                 beast_params=['completeness',
                                               'Av'])#,'Rv','f_A'])
    av_grid = np.unique(beast_data['Av'])

    # also make a more finely sampled A_V grid
    if not log_scale:
        av_grid_big = np.linspace(np.min(av_grid), np.max(av_grid), 500)
    else:
        av_grid_big = np.geomspace(np.min(av_grid), np.max(av_grid), 500)

    # read in the nstars image
    nstars_image, nstars_header = fits.getdata(mb_settings['nstars_filename'], header=True)    
    # dimensions of images/plotting
    y_dimen = nstars_image.shape[0]
    x_dimen = nstars_image.shape[1]

    # read in the best fits
    label_list = mb_settings['fit_param_names']
    best_fits = {}
    for label in label_list:
        with fits.open('./'+projectname+'_megabeast/'+projectname+'_'+label+'_bestfit.fits') as hdu:
            best_fits[label] = hdu[0].data


    # set colors for plots
    cmap = matplotlib.cm.get_cmap('inferno')
    color_data = cmap(0.0)
    color_fit = cmap(0.5)
    if input_lognormal is not None:
        color_input = cmap(0.85)
    
    # -----------------
    # plotting
    # -----------------


    # set up figure
    fig = plt.figure(figsize=(x_dimen*2,y_dimen*2))


    
    for i in tqdm(range(y_dimen), desc='y pixels'):
        for j in tqdm(range(x_dimen), desc='x pixels'):
    #for i in [0]:
    #    for j in [12]:

            if nstars_image[i,j] > 20:


                # -------- data

                # read in the original lnp data
                lnp_filename = mb_settings['lnp_file_prefix']+"_%i_%i_lnp.hd5"%(j, i)
                lnp_data = read_lnp_data(lnp_filename, nstars_image[i,j])
                lnp_vals = np.array(lnp_data['vals'])

                # completeness for each of the values
                lnp_comp = beast_data['completeness'][lnp_data['indxs']]

                # best A_V for each star
                best_av = []
                for k in range(lnp_vals.shape[1]):
                    vals = lnp_vals[:,k]
                    lnp_vals[:,k] = np.log(np.exp(vals) / np.sum(np.exp(vals)))
                    inds = lnp_data['indxs'][:,k]
                    best_val_ind = np.where(vals == np.max(vals))[0][0]
                    best_av.append(beast_data['Av'][inds[best_val_ind]])
                best_av = np.array(best_av)
                

                # stack up some representation of what's being maximized in ensemble_model.py
                prob_stack = np.sum(lnp_comp * np.exp(lnp_vals), axis=1)

                # normalize it (since it's not clear what the numbers mean anyway)
                #prob_stack = prob_stack / np.sum(prob_stack)
                prob_stack = prob_stack / np.trapz(prob_stack, av_grid)

                ## stack up the probabilities at each A_V
                #prob_stack = np.sum(np.exp(lnp_vals), axis=1)

                # set up the subplot
                plt.subplot(y_dimen, x_dimen, (y_dimen-i-1)*(x_dimen) + j + 1)

                # make a histogram
                if not log_scale:
                    plt.plot(av_grid, prob_stack,
                                 marker='.', ms=0, mew=0,
                                 linestyle='-', color=color_data, linewidth=4)                   
                if log_scale:
                    plt.plot(np.log10(av_grid), prob_stack,
                                 marker='.', ms=0, mew=0,
                                 linestyle='-', color=color_data, linewidth=4)
                    
                ax = plt.gca()



                # -------- input lognormal(s)

                if input_lognormal is not None:

                    # create lognormal
                    lognorm = _lognorm(av_grid_big, input_lognormal['max_pos'],
                                           input_lognormal['sigma'], input_lognormal['N'])

                    # if there's a second lognormal
                    if input_lognormal2 is not None:
                        lognorm += _lognorm(av_grid_big, input_lognormal2['max_pos'],
                                           input_lognormal2['sigma'], input_lognormal2['N'])


                    # normalize it
                    #lognorm = lognorm / np.sum(lognorm)
                    lognorm = lognorm / np.trapz(lognorm, av_grid_big)

                    # plot it
                    #yrange_before = ax.get_ylim()
                    if not log_scale:
                        plt.plot(av_grid_big, lognorm,
                                    marker='.', ms=0, mew=0, linestyle='-', color=color_input, linewidth=2, alpha=0.85)
                    if log_scale:
                        plt.plot(np.log10(av_grid_big), lognorm,
                                    marker='.', ms=0, mew=0, linestyle='-', color=color_input, linewidth=2, alpha=0.85)
                    #ax.set_ylim(yrange_before)
                    

                # -------- best fit

                # generate best fit
                lognorm = _two_lognorm(av_grid_big,
                                           best_fits['Av1'][i,j], best_fits['Av2'][i,j],
                                           sigma1=best_fits['sigma1'][i,j],
                                           sigma2=best_fits['sigma2'][i,j],
                                           N1=nstars_image[i,j]*(1 - 1/(best_fits['N12_ratio'][i,j]+1)),
                                           N2=nstars_image[i,j]/(best_fits['N12_ratio'][i,j]+1) )

                # normalize it
                #lognorm = lognorm / nstars_image[i,j]
                #lognorm = lognorm / np.sum(lognorm)
                lognorm = lognorm / np.trapz(lognorm, av_grid_big)
                
                # plot it
                yrange_before = ax.get_ylim()
                if not log_scale:
                    plt.plot(av_grid_big, lognorm,
                                marker='.', ms=0, mew=0, dashes=[3,1.5], color=color_fit, linewidth=2)
                if log_scale:
                    plt.plot(np.log10(av_grid_big), lognorm,
                                marker='.', ms=0, mew=0, dashes=[3,1.5], color=color_fit, linewidth=2)
                ax.set_ylim(yrange_before)


                
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    if not log_scale:
        plt.xlabel(r'$A_V$', size=15)
    else:
        plt.xlabel(r'Log $A_V$', size=15)
    plt.ylabel('PDF', size=15)
    plt.tight_layout()
    
    # save figure
    plt.savefig('./'+projectname+'_megabeast/'+projectname+'_bestfit_plot.pdf')
    plt.close()
