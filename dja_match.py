import os
import numpy as np

# astropy functions for joining tables
from astropy.table import Table, hstack, vstack, join

# grizli and eazy
import grizli
import grizli.catalog
from grizli import utils
from grizli.aws import db # grizli.aws for SQL queries
import eazy

import warnings
warnings.filterwarnings('ignore')


# --------------------------------------------------------------
# ------------------------------- READ SPECTRSCOPY FROM DATABASE
# --------------------------------------------------------------

def spec_tab_from_db():
    """
    Downloads basic info for spectra on the AWS database
    
    Parameters
    ----------
    None

    Returns
    -------
    spec_tab : table containing spectra info, format=astropy Table
    """
    spec_tab = db.SQL("""
                SELECT ne.root, ne.file, ne.ra, ne.dec
                FROM nirspec_extractions ne
                """)
    return(spec_tab)


# --------------------------------------------------------------
# --------------------- READ UMNATCHED SPECTRSCOPY FROM DATABASE
# --------------------------------------------------------------


def unmatched_spec_from_db():
    """
    Downloads basic info for spectra on the AWS database
    that have NOT already been matched to photometry
    
    Parameters
    ----------
    None

    Returns
    -------
    spec_tab : table containing spectra from nirspec_extractions that are NOT in nirspec_phot_match, format=astropy Table
    """
    spec_tab = db.SQL("""
                SELECT ne.root, ne.file, ne.ra, ne.dec
                FROM nirspec_extractions ne
                WHERE ne.file NOT IN (SELECT DISTINCT nrp.file_spec
                                        FROM nirspec_phot_match nrp
                                        )
                """)
    return(spec_tab)


# --------------------------------------------------------------
# ------------------------- READ GRIZLI PHOTOMETRY FROM DATABASE
# --------------------------------------------------------------


def grizli_phot_from_db():
    """
    Downloads basic photometry info on the AWS database
    
    Parameters
    ----------
    None

    Returns
    -------
    phot_tab : table containing phot info from grizli_photometry, format=astropy Table
    """
    phot_tab = db.SQL("""
                SELECT gr.id_phot as id, gr.ra_phot as ra, gr.dec_phot as dec, gr.file_phot
                FROM grizli_photometry gr
                """)
    return(phot_tab)



# --------------------------------------------------------------
# -------------------------------- SPECTROSCOPY-PHOTOMETRY MATCH
# --------------------------------------------------------------

def match_spec_phot_cats(
photcat_concat=None,
spec_tab=None,
sep=30,
outfile_name='spec_phot_match_temp',
path_to_outfile='cats/'
):
    """
    Matches spectroscopic and photometric data via an [ra,dec] match,
    writes output to a 'fits' file
    
    Parameters
    ----------
    photcat_concat : concatenated photcats, format=astropy Table
    spec_tab : spectroscopy table, format=astropy Table
    sep : matching radius in arcsec, format=float
    outfile_name : name of output file, format=string
    path_to_outfile : path to output files, format=string

    Returns
    -------
    empty
    """

    photcat_list = np.unique(photcat_concat['file_phot'])

    # looping over photcat files in order to marry the cats
    for i in range(len(photcat_list)):
    
        # name of photcat file
        photcat_file = photcat_list[i]
        # subset of photcat
        photcat_sub = photcat_concat[photcat_concat['file_phot']==photcat_file]
    
        # match to spec table
        # find nearest match to spec source
        idx, dr, dra, ddec = photcat_sub.match_to_catalog_sky(spec_tab, get_2d_offset=True)
        # mask for sources with < 30 arcsec separation
        _dr_mask_ = dr.value<sep
    
        # marrying the catalogues
        married_tab = hstack([spec_tab[_dr_mask_]['root','file'], photcat_sub[idx][_dr_mask_]['file_phot','id']])
        married_tab.add_column(dra[_dr_mask_], name='dra')
        married_tab.add_column(ddec[_dr_mask_], name='ddec')
        married_tab.add_column(dr[_dr_mask_], name='dr')
        married_tab.rename_column('id', 'id_phot')
        married_tab.rename_column('file', 'file_spec')
        married_tab.rename_column('root', 'root_spec')
        
        if i==0:
            marriedcat_concat = married_tab
        elif i>0:
            marriedcat_concat = vstack([marriedcat_concat,married_tab])

    marriedcat_concat.write(path_to_outfile+outfile_name+'.fits', format='fits', overwrite=True)

    return()


    

# --------------------------------------------------------------
# --------------------------------- CONCATENATING ALL PHOTOMETRY
# --------------------------------------------------------------

def concat_photometry(
photcat_list=None,
outfile_name='photcat_concat_temp',
path_to_outfile='cats/'
):
    """
    Concatenates ['photcat_name','id','ra','dec'] from each photometric catalogue in 'photcat_list'
    by looping over 'photcat_list' and appending info to a master file;
    writes output to a 'fits' file
    
    Parameters
    ----------
    photcat_list : list of names of photometric catalogues on the DJA, format=astropy Table
    outfile_name : name of output file, format=string
    path_to_outfile : path to output files, format=string

    Returns
    -------
    empty
    """

    for i in range(len(photcat_list)):
    
        # read i'th photometric catalogue
        photcat = utils.read_catalog(photcat_list['download'][i])
    
        # extract relevant info from phot file
        photcat_sub = photcat['id','ra','dec']
        # make list of photcat filename
        photfield_list = [str(photcat_list['file'][i])] * len(photcat_sub)
        # add photcat filename to table
        photcat_sub.add_column(photfield_list, name='file_phot', index=0)
    
        if i==0:
            photcat_concat = photcat_sub
        elif i>0:
            photcat_concat = vstack([photcat_concat,photcat_sub])
    
    photcat_concat.write(path_to_outfile+outfile_name+'.fits',format='fits', overwrite=True)

    return()


# --------------------------------------------------------------
# ---------------------------------------- PHOTOMETRY-EAZY MATCH
# --------------------------------------------------------------

def match_grizli_phot_to_eazy(
photcat_list=None,
path_to_files = '../cats/zout_data',
outfile_name='phot_zout_match_temp',
path_to_outfile='cats/'
):
    """
    Matches photometric data to EAZY output from the DJA via filenames and source IDs,
    writes output to a 'fits' file

    Parameters
    ----------
    photcat_list : list of names of photometric catalogues on the DJA, format=astropy Table
    path_to_files: path to input files, format=string
    outfile_name : name of output file, format=string
    path_to_outfile : path to output files, format=string
    
    Returns
    -------
    empty    
    """

    for i in range(len(photcat_list)):
        
        # get field name from file name
        field = (photcat_list[i]['file']).split('-fix')[0]

        """ opening and reading zout files """
        zout_filename = field+'-fix.eazypy.zout.fits'
        zout_tab = Table.read(f'{path_to_files}/{zout_filename}', format='fits')
        # remove certain columns from zout file
        zout_cols2remove = ['mass_p','sfr_p','Lv_p','LIR_p','energy_abs_p','Lu_p','Lj_p',
                            'L1400_p','L2800_p','LHa_p','LOIII_p','LHb_p','LOII_p','Av_p','ssfr_p']
        zout_tab.remove_columns(zout_cols2remove)
    
        """ opening and reading phot files """
        phot_filename = field+'-fix_phot_apcorr.fits'
        phot_file = Table.read(f'{path_to_files}/{phot_filename}', format='fits')
        # colnames of cols to keep from photometric catalogue
        phot_colnames = ['id', 'thresh', 'npix', 'tnpix', 'xmin', 'xmax', 'ymin', 'ymax',
                         'x2_image', 'y2_image', 'xy_image', 'errx2', 'erry2', 'errxy',	
                         'a_image', 'b_image', 'theta_image', 'cxx_image', 'cyy_image', 'cxy_image',
                         'cflux', 'flux', 'cpeak', 'peak', 'xcpeak', 'ycpeak', 'xpeak', 'ypeak', 'x_image', 'y_image',
                         'ra', 'dec',
                         'flux_iso', 'fluxerr_iso', 'area_iso', 'mag_iso',
                         'kron_radius', 'kron_rcirc',
                         'flux_auto', 'fluxerr_auto', 'flux_radius_20', 'flux_radius', 'flux_radius_90',
                         'area_auto', 'mag_auto', 'magerr_auto']
        # filter colnames of cols to keep from photometric catalogue
        filt_names = ['f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w'] # filter names
        ext1 = ['_tot', '_etot', '_mask_aper', '_flux_aper', '_fluxerr_aper'] # colname extensions 1
        ext2 = ['_0', '_1', '_2'] # colname extensions 2
        # generating all filt colnames
        filt_colnames = [filt_names_i + ext_j + ext_k for filt_names_i in filt_names for ext_k in ext2 for ext_j in ext1]
        # concatenating all phot colnames
        all_phot_colnames = phot_colnames+filt_colnames
        phot_colnames_in_phot_file =  [phot_col for phot_col in all_phot_colnames if phot_col in phot_file.colnames]
    
        phot_file_sub = phot_file[phot_colnames_in_phot_file]
        phot_file_sub.rename_column('id', 'id_phot')
        phot_file_sub.rename_column('ra', 'ra_phot')
        phot_file_sub.rename_column('dec', 'dec_phot')
        phot_file_sub.add_column([phot_filename]*len(phot_file_sub), name='file_phot', index=0)
        phot_file_sub.add_column([zout_filename]*len(zout_tab), name='file_zout', index=1)
    
        # matching phot and zout tables based on 'ids'
        phot_zout_match_tab = join(phot_file_sub,zout_tab, keys_left=['id_phot'], keys_right=['id'], join_type='inner')
    
        if i==0:
            phot_zout_match_concat = phot_zout_match_tab
        elif i>0:
            phot_zout_match_concat = vstack([phot_zout_match_concat,phot_zout_match_tab])
        else:
            print('Check index in for loop!')

    phot_zout_match_concat.write(path_to_outfile+outfile_name+'.fits', format='fits', overwrite=True)
    
    return()

