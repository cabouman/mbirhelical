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
from io import StringIO
import socket
import psutil
import re


def ceil(a, b):
        return -1 * (-a // b)




"""
Setup communication routines: Credit to https://github.com/ORNL/HydraGNN
"""

def init_comm_size_and_rank():
    world_size = None
    world_rank = 0
    if os.getenv("OMPI_COMM_WORLD_SIZE") and os.getenv("OMPI_COMM_WORLD_RANK"):
        ## Summit
        world_size = int(os.environ["OMPI_COMM_WORLD_SIZE"])
        world_rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
    elif os.getenv("SLURM_NPROCS") and os.getenv("SLURM_PROCID"):
        ## CADES, Frontier, Perlmutter
        world_size = int(os.environ["SLURM_NPROCS"])
        world_rank = int(os.environ["SLURM_PROCID"])
    ## Fall back to default
    if world_size is None:
        world_size = 1
    return int(world_size), int(world_rank)


def find_ifname(myaddr):
    """
    Find socket ifname for a given ip adress. This is for "GLOO" ddp setup.
    Usage example:
        find_ifname("127.0.0.1") will return a network interface name, such as "lo". "lo0", etc.
    """
    ipaddr = socket.gethostbyname(myaddr)
    ifname = None
    for nic, addrs in psutil.net_if_addrs().items():
        for addr in addrs:
            if addr.address == ipaddr:
                ifname = nic
                break
        if ifname is not None:
            break
    return ifname


def parse_slurm_nodelist(nodelist):
    """
    Parse SLURM_NODELIST env string to get list of nodes.
    Usage example:
        parse_slurm_nodelist(os.environ["SLURM_NODELIST"])
    Input examples:
        "or-condo-g04"
        "or-condo-g[05,07-08,13]"
        "or-condo-g[05,07-08,13],or-condo-h[01,12]"
    """
    nlist = list()
    for block, _ in re.findall(r"([\w-]+(\[[\d\-,]+\])*)", nodelist):
        m = re.match(r"^(?P<prefix>[\w\-]+)\[(?P<group>.*)\]", block)
        if m is None:
            ## single node
            nlist.append(block)
        else:
            ## multiple nodes
            g = m.groups()
            prefix = g[0]
            for sub in g[1].split(","):
                if "-" in sub:
                    start, end = re.match(r"(\d+)-(\d+)", sub).groups()
                    fmt = "%%0%dd" % (len(start))
                    for i in range(int(start), int(end) + 1):
                        node = prefix + fmt % i
                        nlist.append(node)
                else:
                    node = prefix + sub
                    nlist.append(node)
    return nlist

def setup_ddp():
    """ "Initialize DDP"""
    if os.getenv("DDSTORE_BACKEND") is not None:
        backend = os.environ["DDSTORE_BACKEND"]
    elif dist.is_nccl_available() and torch.cuda.is_available():
        backend = "nccl"
    elif torch.distributed.is_gloo_available():
        backend = "gloo"
    else:
        raise RuntimeError("No parallel backends available")
    world_size, world_rank = init_comm_size_and_rank()
    ## Default setting
    master_addr = os.getenv("MASTER_ADDR", "127.0.0.1")
    master_port = os.getenv("MASTER_PORT", "8889")
    if os.getenv("LSB_HOSTS") is not None:
        ## source: https://www.olcf.ornl.gov/wp-content/uploads/2019/12/Scaling-DL-on-Summit.pdf
        ## The following is Summit specific
        master_addr = os.environ["LSB_HOSTS"].split()[1]
    elif os.getenv("LSB_MCPU_HOSTS") is not None:
        master_addr = os.environ["LSB_MCPU_HOSTS"].split()[2]
    elif os.getenv("SLURM_NODELIST") is not None:
        ## The following is CADES/Frontier/Perlmutter specific
        master_addr = parse_slurm_nodelist(os.environ["SLURM_NODELIST"])[0]
    try:
        if backend in ["nccl", "gloo"]:
            os.environ["MASTER_ADDR"] = master_addr
            os.environ["MASTER_PORT"] = master_port
            os.environ["WORLD_SIZE"] = str(world_size)
            os.environ["RANK"] = str(world_rank)
        if (backend == "gloo") and ("GLOO_SOCKET_IFNAME" not in os.environ):
            ifname = find_ifname(master_addr)
            if ifname is not None:
                os.environ["GLOO_SOCKET_IFNAME"] = ifname
        print("Distributed data parallel: %s master at %s:%s" % (backend, master_addr, master_port))
        if not dist.is_initialized():
            dist.init_process_group(backend=backend, init_method="env://")
    except KeyError:
        print("DDP has to be initialized within a job - Running in sequential mode")
    return world_size, world_rank




world_size = int(os.environ['SLURM_NTASKS'])
world_rank = int(os.environ['SLURM_PROCID'])
local_rank = int(os.environ['SLURM_LOCALID'])

torch.cuda.set_device(local_rank)
device = torch.cuda.current_device()


#dist.init_process_group('nccl', rank=world_rank, world_size=world_size)

setup_ddp()

print("we are here after dist.init_process_group. world_size ",world_size,flush=True)


total_datasets = 1
my_range_start =  np.clip(ceil(total_datasets,world_size) * world_rank,a_min=None, a_max = total_datasets)

if world_rank != (world_size-1):
    my_range_end = np.clip(ceil(total_datasets,world_size) * (world_rank+1),a_min=None, a_max = total_datasets)
else:
    my_range_end = total_datasets


print("world_rank ",world_rank," my_range_start ",my_range_start," my_range_end ",my_range_end)


channels = 900
views = 5000
rows = 64
for data_idx in range(my_range_start,my_range_end):
    DIR="/lustre/orion/nro108/scratch/xf9/result%03d_Chest_duke1_test_TCM_Chest_23BMI.xcat" % (data_idx)

    print("DIR ",DIR)


    full_sino= torch.zeros(channels,rows,views)

    with open(DIR,'rb') as f:
        full_sino = torch.from_numpy(np.reshape(np.fromfile(f, dtype=np.float32),(channels, rows, views)))

        print("full_sino.shape",full_sino.shape,flush=True)     


    print("full_sino max ",torch.max(full_sino)," min ",torch.min(full_sino), " dtype ",full_sino.dtype, " shape ",full_sino.shape)
    
    #views, channels, rows
    full_sino=torch.permute(full_sino, (2, 0, 1))

    print("full_sino after permute ",full_sino.shape)

    
    weight = torch.exp(-0.2*full_sino)



    outputname=StringIO("/lustre/orion/nro108/scratch/xf9/sino/chest_%03d_proj.sino" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('%d ' % NumRows)
        f.write('%d ' % NumChannels)
        f.write('%d\n' % NumViews)


    with open(outputname.getvalue(),'ab') as f:
        np.array(focal1_Proj.cpu().numpy(),dtype=np.float32).tofile(f)


    outputname=StringIO("/lustre/orion/nro108/scratch/xf9/wght/chest_%03d_weight.wght" % data_idx)

    with open(outputname.getvalue(),'w') as f:
        f.write('%d ' % NumRows)
        f.write('%d ' % NumChannels)
        f.write('%d\n' % NumViews)


    with open(outputname.getvalue(),'ab') as f:
        np.array(weight.cpu().numpy(),dtype=np.float32).tofile(f)
