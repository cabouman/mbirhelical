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

device = torch.device('cuda:0')


data_idx=0
DIR="/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_%03d/" % (data_idx)

print("DIR ",DIR)
proj_name=DIR+"proj_0001.dcm"
proj_info = pydicom.dcmread(proj_name)


print("length ",len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]))


NumViews=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])-1   
NumRows=proj_info.Columns
NumChannels=proj_info.Rows
RescaleSlope = proj_info.RescaleSlope
RescaleIntercept = proj_info.RescaleIntercept

print("NumRows ",NumRows, " NumChannels ",NumChannels," NumViews ",NumViews," RescaleSlope ",RescaleSlope," RescaleIntercept ",RescaleIntercept)
print(type(NumViews)," ",type(NumRows)," ",type(NumChannels),type(NumViews))





focal1_Proj=torch.zeros(NumChannels,NumRows,NumViews)
for iv in range (0, NumViews):
    proj="proj_%04d.dcm" % (iv+1)
    projname=DIR+proj
    temp = pydicom.dcmread(projname).pixel_array
    focal1_Proj[:,:,iv] = torch.from_numpy(temp *RescaleSlope +RescaleIntercept)

print("focal1_Proj max ",torch.max(focal1_Proj)," min ",torch.min(focal1_Proj), " dtype ",focal1_Proj.dtype, " shape ",focal1_Proj.shape)
focal1_Proj=torch.permute(focal1_Proj, (2, 0, 1))

print("focal1_Proj after permute ",focal1_Proj.shape)

weight = torch.exp(-0.3*focal1_Proj)


outputname=StringIO("/gpfs/alpine/med106/world-shared/xf9/aapm-preprocess/dcm%03d_proj.sino" % data_idx)

with open(outputname.getvalue(),'w') as f:
    f.write('%d ' % NumRows)
    f.write('%d ' % NumChannels)
    f.write('%d\n' % NumViews)


with open(outputname.getvalue(),'ab') as f:
    np.array(focal1_Proj.cpu().numpy(),dtype=np.float32).tofile(f)


outputname=StringIO("/gpfs/alpine/med106/world-shared/xf9/aapm-preprocess/dcm%03d_weight.wght" % data_idx)

with open(outputname.getvalue(),'w') as f:
    f.write('%d ' % NumRows)
    f.write('%d ' % NumChannels)
    f.write('%d\n' % NumViews)


with open(outputname.getvalue(),'ab') as f:
    np.array(weight.cpu().numpy(),dtype=np.float32).tofile(f)



