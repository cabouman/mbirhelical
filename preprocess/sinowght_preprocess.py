import json
import argparse
import os
import subprocess
import torch
from torch import nn, optim
from torch.utils.data import DataLoader
from torch.nn import functional as F
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
import time
import matplotlib.pyplot as plt
from scipy.io import savemat
import math
import numpy as np
import pydicom
from io import StringIO

get_master = "echo $(cat {} | sort | uniq | grep -v batch | grep -v login | head -1)".format(os.environ['LSB_DJOB_HOSTFILE'])
os.environ['MASTER_ADDR'] = str(subprocess.check_output(get_master, shell=True))[2:-3]
os.environ['MASTER_PORT'] = "23456"
os.environ['WORLD_SIZE'] = os.environ['OMPI_COMM_WORLD_SIZE']
os.environ['RANK'] = os.environ['OMPI_COMM_WORLD_RANK']

world_size = int(os.environ['OMPI_COMM_WORLD_SIZE'])
world_rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
local_rank = int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])


device = torch.device('cuda:{}'.format(local_rank))


dist.init_process_group('nccl', rank=world_rank, world_size=world_size)



print("reach here ",world_rank)

total_datasets = 200
my_range_start =  total_datasets//world_size * world_rank
if world_rank != (world_size-1):
    my_range_end = np.clip(total_datasets//world_size * (world_rank+1),a_min=None, a_max = total_datasets)
else:
    my_range_end = total_datasets


print("world_rank ",world_rank," my_range_start ",my_range_start," my_range_end ",my_range_end)

for data_idx in range(my_range_start,my_range_end):
    if data_idx < 67:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_%03d/" % (data_idx)
    elif data_idx>=67 and data_idx <134:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_lung_lesion/dcm_%03d/" % (data_idx)
    else:
        DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_liver/dcm_%03d/" % (data_idx)

    print("DIR ",DIR)

    print("length ",len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))


    NumViews=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])-1  

    if NumViews<10000:
        proj_name=DIR+"proj_0001.dcm"
    else:
        proj_name=DIR+"proj_00001.dcm"
    proj_info = pydicom.dcmread(proj_name)


    NumRows=proj_info.Columns
    NumChannels=proj_info.Rows
    RescaleSlope = proj_info.RescaleSlope
    RescaleIntercept = proj_info.RescaleIntercept

    print("NumRows ",NumRows, " NumChannels ",NumChannels," NumViews ",NumViews," RescaleSlope ",RescaleSlope," RescaleIntercept ",RescaleIntercept)
    print(type(NumViews)," ",type(NumRows)," ",type(NumChannels),type(NumViews))





    focal1_Proj=torch.zeros(NumChannels,NumRows,NumViews)
    for iv in range (0, NumViews):
        if NumViews<10000:
            proj="proj_%04d.dcm" % (iv+1)
        else:
            proj="proj_%05d.dcm" % (iv+1)
        projname=DIR+proj
        temp = pydicom.dcmread(projname).pixel_array
        focal1_Proj[:,:,iv] = torch.from_numpy(temp *RescaleSlope +RescaleIntercept)

    print("focal1_Proj max ",torch.max(focal1_Proj)," min ",torch.min(focal1_Proj), " dtype ",focal1_Proj.dtype, " shape ",focal1_Proj.shape)
    focal1_Proj=torch.permute(focal1_Proj, (2, 0, 1))

    print("focal1_Proj after permute ",focal1_Proj.shape)


