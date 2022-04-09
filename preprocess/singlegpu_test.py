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



proj_info = pydicom.dcmread("/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_000/proj_0001.dcm")



NumViews=9000;   
NumRows=proj_info.Columns;
NumChannels=proj_info.Rows;
RescaleSlope = proj_info.RescaleSlope;
RescaleIntercept = proj_info.RescaleIntercept;
#proj_info = dicominfo('proj_0001.dcm',"dictionary",'C:/Users/xiaow/Downloads/AAPM-2022/DICOM-CT-PD-dict_v10.txt');

#print("NumRows ",NumRows, " NumChannels ",NumChannels," RescaleSlope ",RescaleSlope," RescaleIntercept ",RescaleIntercept)

focal1_Proj=torch.zeros(NumChannels,NumRows,NumViews).to(device)
for iv in range (0, NumViews):
    dataname=StringIO("/gpfs/alpine/med106/world-shared/irl1/aapm_data/dcmproj_copd/dcm_000/proj_%04d.dcm" % (iv+1))
    
    print(dataname.getvalue())
    temp = pydicom.dcmread(dataname.getvalue()).pixel_array
    focal1_Proj[:,:,iv] = torch.from_numpy(temp *RescaleSlope +RescaleIntercept)

print("focal1_Proj max ",torch.max(focal1_Proj)," min ",torch.min(focal1_Proj), " dtype ",focal1_Proj.dtype)
