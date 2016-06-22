import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm


def calc_density(diam):

    from input_file import *

    diam1=np.ones(npart)*diam1
    diam2=np.ones(npart)*diam2

    rho1 = np.ones(npart)*rho1
    rho2 = np.ones(npart)*rho2

    density = np.zeros(npart)

    for i in range(npart):

        if ( diam[i] < diam1[i] ):

            density[i] = rho1[i]

        elif ( diam[i] < diam2[i] ):

            diam1_phi = -log(diam1(i_part)*1000)/log(2.0)
            diam2_phi = -log(diam2(i_part)*1000)/log(2.0)

            density[i] = rho1[i] + ( diam_phi[i] - diam1_phi ) / ( diam2_phi - diam1_phi ) * ( rho2[i] - rho1[i] )
       
        else:

           density[i] = rho2[i]
       
    return density
