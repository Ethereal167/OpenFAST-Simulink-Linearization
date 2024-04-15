#  此脚本用于测试所构造的风力机线性模型
import os
import numpy as np
import glob
from linear.linear_models import LinearTurbineModel
from lin_util import interp_plant

'''
Load linear turbine models

Parameters:
-----------
linfile_path: string
    Path to folder containing .lin files
load_parllel: bool
    Load parallel True/False
'''
# Parse openfast linearization filenames
load_parallel = False
linfile_path = r"D:\MatlabProject\OpenFast_v2.4.0-v3.0.0\Linearization\NREL5MW\ws14"


filenames = glob.glob(os.path.join(linfile_path, '*.lin'))
linfiles = [os.path.split(file)[1] for file in filenames]
linroots = np.sort(np.unique([file.split('.')[0] for file in linfiles])).tolist()
linroots[0] = ''.join([linroots[0], '.0'])
linfile_numbers = set([int(file.split('.')[2]) for file in linfiles])

# Load linturb
linturb = LinearTurbineModel(linfile_path, linroots,
                             nlin=max(linfile_numbers), rm_hydro=True, load_parallel=load_parallel)

A_ops = linturb.A_ops[:, :, 0]
B_ops = linturb.B_ops[:, :, 0]
C_ops = linturb.C_ops[:, :, 0]
D_ops = linturb.D_ops[:, :, 0]

u = 14
P = interp_plant(linturb, u, return_scipy=False)
A = P.A
B = P.B
C = P.C
D = P.D

diff = D_ops - D