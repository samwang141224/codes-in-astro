#!/usr/bin/env python

import time

from wildhunt import image, catalog, plotting,plottingproposal

from IPython import embed

if __name__ == "__main__":

    t0 = time.time()

    cat = catalog.Catalog('example', 'RA', 'DEC', 'Source_Name', datapath='/Users/samwang/Downloads/Gdrop/test/LOFAR/LOFAR_sel/sel2/nols_candi_sel2_clean.csv')

    
    FOV = 60



    survey_dict = [

        
       
        {'survey': 'DELSDR9', 'bands': ['g','r','z'],
         'fov': FOV},
         {'survey': 'PS1', 'bands': ['y'], 'fov':FOV},
         
         
         
         
         {'survey': 'UNWISE-neo2', 'bands': ['1','2'], 'fov': FOV},
        
        



         
        
        {'survey': 'allWISE', 'bands': ['3', '4'], 'fov': FOV}
    ]

    #{'survey': 'UKIDSSDR11PLUSLAS', 'bands': ['Y','J','H','K'],
          #'fov':FOV}
    #{'survey': 'PS1', 'bands': ['y'], 'fov':FOV}
    #{'survey': 'PS1', 'bands': ['g', 'r', 'i', 'z', 'y'], 'fov':FOV},
    #{'survey': 'UNWISE-neo2', 'bands': ['1','2'], 'fov': FOV},


    ra = cat.df.compute()[cat.ra_colname][:]
    dec = cat.df.compute()[cat.dec_colname][:]
    plotting.generate_cutout_images(ra, dec, survey_dict, imgsize = 60,aperture = 5,n_col = 4,
                                    image_folder_path = '/Users/samwang/Downloads/Gdrop/test/LOFAR/LOFAR_sel/sel2/cutouts',
                                    n_jobs=4)

    #plot spitzer for JWST proposal
    #ra = [11.87188]
    #dec = [0.066349]
    #survey_dict = [{'survey': 'SPITZER', 'bands': ['3.6um'],
    #     'fov': FOV},]
    #plottingproposal.generate_cutout_images_proposal(ra, dec, survey_dict,download_images = False, imgsize = 60,aperture = 3,n_col = 4,
    #                                image_folder_path = '/Users/samwang/Downloads/Gdrop/test/cutouts',
    #                                n_jobs=4)

    print("{:.1f} s: ".format( time.time() - t0))

    '''
    {'survey': 'UKIDSSDR11PLUSLAS', 'bands': ['J'],
          'fov':120}
    '''