# DMRG
Matlab code of DMRG on 1D Heisenberg Model

1维最近邻海森堡哈密顿量

H=J\sum_{ij}\bm{S}_i\bm{S}_j
  
script_DMRG* 是主程序，v2是改进过的，其余是script* 调用的函数

代码很糙，初学的时候写的。直接上了泡利矩阵，没用升降算符。求基态也是用matlab内置的eigen函数，没用Lanzos。

ps 代码中的中文备注不知道为啥网站上显示乱码，下载再看可能没问题。
