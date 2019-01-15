# system
from __future__ import (absolute_import, division, print_function)

# other packages
import numpy as np
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits

# megabeast
from .read_megabeast_input import read_megabeast_input



def parameter_maps(megabeast_input_file, n_col=2):
    """
    Create maps of the best-fit parameters from the MegaBEAST

    Parameters
    ----------
    megabeast_input_file : string
        Name of the file that contains settings, filenames, etc

    n_col : int (default=2)
        number of columns of plots in the output file

    """

    # read in the settings from the file
    mb_settings = read_megabeast_input(megabeast_input_file)

    # get the project name
    projectname = mb_settings['projectname']

    # list of parameters to plot
    plot_params = mb_settings['fit_param_names']

    # aspect ratio of the field
    with fits.open('./'+projectname+'_megabeast/'+projectname+'_'+plot_params[0]+'_bestfit.fits') as hdu:
        im_size = hdu[0].data.shape
    
    
    # initialize figure
    size = 4
    n_row = len(plot_params)//n_col + np.count_nonzero(len(plot_params) % n_col)
    
    fig = plt.figure(figsize=(size*n_col, size*im_size[0]/im_size[1]*n_row ))


    
    for p,param in enumerate(plot_params):

        print(param)

        # image file name
        im_file = './'+projectname+'_megabeast/'+projectname+'_'+param+'_bestfit.fits'

        # subplot info
        # - current column and row
        #col_num = p % n_col
        #row_num = p // n_col
        # - corresponding dimensions
        #subplot_dimen = [1/n_col * col_num, 1-(1/n_row * (row_num+1)),
        #                     1/n_col * (col_num+1), 1-(1/n_row * row_num)]
        #print('p='+str(p)+' col_num='+str(col_num)+' row_num='+str(row_num))
        #print(subplot_dimen)

        #pdb.set_trace()


        f1 = aplpy.FITSFigure(im_file, figure=fig, subplot=(n_row, n_col, p+1))
        #f1 = aplpy.FITSFigure(im_file, figure=fig, subplot=subplot_dimen)
        f1.show_colorscale(cmap='magma')
        f1.add_colorbar()
        #f1.colorbar.set_box([0.12, 0.04, 0.33, 0.02], box_orientation='horizontal') # [xmin, ymin, dx, dy]
        f1.colorbar.set_font(size=10)
        f1.axis_labels.set_xtext(param)
        f1.axis_labels.set_font(size=15)
        #f1.axis_labels.hide_x()
        f1.axis_labels.hide_y()
        f1.tick_labels.hide_x()
        f1.tick_labels.hide_y()
        f1.frame.set_linewidth(0)
        plt.tight_layout()


        
    plt.savefig('./'+projectname+'_megabeast/'+projectname+'_bestfit_maps.pdf')
    plt.close()
