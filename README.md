# OpenFAST-Simulink-Linearization
OpenFAST与Simulink联合仿真，NREL5MW风电机组在Simulink中的线性化模型

仿真基于OpenFAST v3.4.0
1. 打开ROSCO-2.7.0\Matlab_Toolbox文件夹下的runFAST.m，根据自己电脑的路径重新设置文件路径，点击运行即可
2. 线性化模型\5MW_Land_BD_Linear文件夹下的RunLinear.m用于生成工作点线性化的状态空间矩阵
3. GenerateStateSpace.m函数用于将不同的工作点对应的状态空间矩阵和xop，yop，uop整合在一起，生成StateSpaceList.mat结构体，提供给simulink仿真
线性化相关操作（线性化过程可以看我b站上的视频 https://www.bilibili.com/video/BV1Zh4y1J7zD/?spm_id_from=333.999.0.0&vd_source=a515334638b633a92b8ef6eeb51abee8 ）
