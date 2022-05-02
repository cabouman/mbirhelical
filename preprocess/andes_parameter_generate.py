import json
import argparse
import os
import subprocess
import time
import matplotlib.pyplot as plt
from scipy.io import savemat
import math
import numpy as np
import pydicom
from io import StringIO
from math import ceil



for data_idx in range(0,200):
    if data_idx < 67:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_%03d/" % (data_idx)
    elif data_idx>=67 and data_idx <134:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_lung_lesion/dcm_%03d/" % (data_idx)
    else:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_liver/dcm_%03d/" % (data_idx)

    print("DIR ",DIR)

    print("length ",len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))

    #150 has no tube current
    if data_idx!=150:
        NumViews=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])-1  
    else:
        NumViews=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])


    if NumViews<10000:
        proj_name=DIR+"proj_0001.dcm"
    else:
        proj_name=DIR+"proj_00001.dcm"
    proj_info = pydicom.dcmread(proj_name)


    NumRows=proj_info.Columns
    NumChannels=proj_info.Rows

    print("NumRows ",NumRows, " NumChannels ",NumChannels," NumViews ",NumViews)
    print(type(NumViews)," ",type(NumRows)," ",type(NumChannels),type(NumViews))


    if data_idx <134:
        z_spacing = 0.5476
    else:
        z_spacing = 1.5


    #create the main parent folder and the detector1 folder
    parent_folder="/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/detector1" % (data_idx)

    from pathlib import Path
    Path(parent_folder).mkdir(parents=True, exist_ok=True)



    recon_z_start = 0.0
    recon_z_end = 0.0
    z_file=open("/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_zPositionList.txt" % (data_idx))   
    lines_to_read=[500, NumViews-500-1]
    for position, line in enumerate(z_file):
        if position == lines_to_read[0]:
            recon_z_start = float(line)
        elif position == lines_to_read[1]:
            recon_z_end = float(line)

    print("recon_z_start ",recon_z_start," recon_z_end ",recon_z_end)


    outputname=StringIO("/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/detector1/geom_recon.txt" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('number of rows\n')
        f.write('%d\n\n' % NumRows)
        f.write('number of channels\n')
        f.write('%d\n\n' % NumChannels)
        f.write('number of views\n')
        f.write('%d\n\n' % NumViews)
        f.write('views per rotation\n')
        f.write('%d\n\n' % 1000)
        f.write('src to iso (mm)\n')
        f.write('%d\n\n' % 575)
        f.write('src to det (mm)\n')
        f.write('%f\n\n' % 1050.0015)
        f.write('un-normalized pitch (det rows/rot)\n')
        f.write('%f\n\n' % 35.0475)
        f.write('X-ray source initial z position (mm)\n')
        f.write('%f\n\n' % 0)
        f.write('view angle spacing (rad)\n')
        f.write('%.10f\n\n' % 0.006283185307)
        f.write('detector channel angle spacing (rad)\n')
        f.write('%.10f\n\n' % 0.00095237959)
        f.write('detector channel offset (rad)\n')
        f.write('%.10f\n\n' % 0.001666664282500)
        f.write('detector row width (mm)\n')
        f.write('%d\n\n' % 1)
        f.write('detector row offset (mm)\n')
        f.write('%d\n\n' % 0)
        f.write('diameter of field of view (mm)\n')
        f.write('%d\n\n' % 500)
        f.write('initial photon counts (0 for noiseless)\n')
        f.write('%f\n\n' % 2.25)
        f.write('electronic noise variance (0 for noiseless)\n')
        f.write('%f\n\n' % 0.3)
        f.write('sinogram file location\n')
        sino_DIR="/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_proj.sino" % (data_idx)        
        f.write(sino_DIR+'\n\n')
        f.write('weight file location\n')
        wght_DIR="/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_weight.wght" % (data_idx)   
        f.write(wght_DIR+'\n\n')
        f.write('dosage file location\n')
        f.write('NA\n\n')
        f.write('offset file location\n')
        f.write('NA\n\n')
        f.write('detector mask location\n')
        f.write('NA\n\n')
        f.write('view angles list\n')
        view_DIR="/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_viewAnglesList.txt" % (data_idx)   
        f.write(view_DIR+'\n\n')
        f.write('source z position list\n')
        z_DIR="/gpfs/alpine/gen006/scratch/xf9/aapm-preprocess/dcm%03d_zPositionList.txt" % (data_idx)   
        f.write(z_DIR+'\n')


    outputname=StringIO("/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/ce.txt" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('SigmaLambda\n')
        if data_idx <134:
            f.write('%f\n\n' % 0.10)
        else:
            f.write('%f\n\n' % 0.15)
        
        f.write('Rho Consensus (damping):\n')
        f.write('%f\n' % 0.8)




    outputname=StringIO("/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/forward_model_directory.txt" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('../data/aapm-parameters/dcm_%03d/detector1/geom_recon.txt\n' % data_idx)


    outputname=StringIO("/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/info_recon.txt" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('number of voxels in x\n')
        f.write('%d\n\n' % 512)
        f.write('number of voxels in y\n')
        f.write('%d\n\n' % 512)
        f.write('number of voxels in z (good slices)\n')
        useful_z_slices = ceil((recon_z_end - recon_z_start)/z_spacing)+1
        if data_idx<134:
            total_z_slices = useful_z_slices+ceil(35.0475*2/z_spacing)+10
        else:
            total_z_slices = useful_z_slices+ceil(35.0475*2/z_spacing)

        f.write('%d\n\n' % total_z_slices)
        f.write('number of voxels in z (total)\n')
        f.write('%d\n\n' % total_z_slices)
        f.write('x coordinate of the center voxel (mm)\n')
        f.write('%d\n\n' % 0)
        f.write('y coordinate of the center voxel (mm)\n')
        f.write('%d\n\n' % 0)
        f.write('z coordinate of the center voxel (mm)\n')

        print("total ",total_z_slices," useful ",useful_z_slices)
        if total_z_slices%2==1:
            z_center = recon_z_start+(useful_z_slices-1)//2*z_spacing
        else:
            z_center = z_spacing/2+recon_z_start+(useful_z_slices-1)//2*z_spacing
        f.write('%f\n\n' % z_center)
        f.write('voxel spacing in xy (mm)\n')
        f.write('%f\n\n' % 0.976)
        f.write('voxel spacing in z (mm)\n')
        f.write('%f\n\n' % z_spacing)
        f.write('radius of the reconstruction mask (mm)\n')
        if data_idx<134:
            f.write('%d\n\n' % 290)
        else:
            f.write('%d\n\n' % 250)

        f.write('initial reconstruction image location\n')

        recon_dir="/gpfs/alpine/gen006/scratch/xf9/recon/dcm%03d/recon" % (data_idx)   
        #recon_dir="NA"

        f.write(recon_dir+'\n\n')
        f.write('mask file location\n')
        f.write('NA\n')



    skipped_slices =0
    if total_z_slices%2==0:
        skipped_slices = (recon_z_start-(z_center-z_spacing/2-(total_z_slices/2-1)*z_spacing)+1e-4)//z_spacing
    else:
        skipped_slices = (recon_z_start-(z_center-(total_z_slices//2)*z_spacing)+1e-4)//z_spacing

    print("skipped_slices ",skipped_slices)

    outputname=StringIO("/gpfs/alpine/gen006/proj-shared/xf9/mbirhelical/data/aapm-parameters/dcm_%03d/prior_qggmrf.txt" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('Q-GGMRF Prior Parameter, q (recommended value is 2) :\n')
        f.write('%f\n\n' % 2)
        f.write('Q-GGMRF Prior Parameter, p (permitted value between 1 to parameter q) :\n')
        f.write('%f\n\n' % 1.2)
        f.write('Q-GGMRF Prior Parameter, T (soft threshold for edges):\n')
        f.write('%f\n\n' % 1)
        f.write('Prior Regularization parameter, sigmaX, (mm^-1) (increasing sigmaX decreases regularization) :\n')
        if data_idx <134:
            f.write('%f\n' % 0.5)
        else:
            f.write('%f\n' % 0.3)
