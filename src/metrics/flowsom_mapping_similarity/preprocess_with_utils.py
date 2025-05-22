import os
import anndata as ad
import numpy as np

import sys
sys.path.append(sys.argv[1])

from helper_functions import get_obs_var_for_integrated, subset_nocontrols, remove_unlabelled


#Load data
input_integrated = ad.read_h5ad(sys.argv[2])
input_unintegrated = ad.read_h5ad(sys.argv[3])
input_validation = ad.read_h5ad(sys.argv[4])
#Format integrated data
input_integrated = get_obs_var_for_integrated(input_integrated,input_validation,input_unintegrated)
input_integrated = subset_nocontrols(input_integrated)
input_integrated = remove_unlabelled(input_integrated)
#Format validation data
input_validation = remove_unlabelled(input_validation)

#print(input_integrated, input_unintegrated, input_validation, sep="\n"*2)

if not os.path.exists('tmp/'):
    os.mkdir('tmp/')
input_integrated.write_h5ad('tmp/integrated.h5ad')
input_validation.write_h5ad('tmp/validation.h5ad')