# large-scale-phase-retrieval
This repository contains the code for a large-scale phase retrieval method via Plug-and-play optimization. For more information, please contact Liheng Bian (bian at bit dot edu dot cn).

![image](https://github.com/bianlab/large-scale-phase-retrieval/blob/main/results/Quantitative_results.JPG)

*Table1. Quantitative comparison under the CDP modality (5 modulations).*

## About large-scale phase retrieval and reported method
High-throughput computational imaging requires efficient processing algorithms to retrieve multi-dimensional and multi-scale information. In computational phase imaging, phase retrieval (PR) is required to reconstruct both amplitude and phase in complex space from intensity-only measurements. The existing PR algorithms suffer from the tradeoff among
low computational complexity, robustness to measurement noise and strong generalization on different modalities. In this work, we report an efficient large-scale phase retrieval technique termed as LPR. It extends the plug-and-play generalized-alternating-projection framework
from real space to nonlinear complex space. The alternating projection solver and enhancing
neural network are respectively derived to tackle the measurement formation and statistical
prior regularization. This framework compensates the shortcomings of each operator, so as
to realize high-fidelity phase retrieval with low computational complexity and strong generalization.


## Usage

Please clone this repository by Git or download the zip file firstly. 

### test_CDP

Run `test_CDP.m` file to achieve the phase retrieval results of conventional algorithms (GS, Fienup, WF, TWF, RWF, AF, TAF, RAF, PMAX, CD, KAC)、prDeep algorithm、and ours algorithm (LPR).
 
### Note
This demo code takes nearly half an hour (2K image).

For different noise level and reslution, please adjustment parameters to acquire better results.


## Platform

All the experiments are implemented using MATLAB 2019b with an Intel i7-9700 processor at 3.0GHz and 16GB RAM. 


## Requirements and Dependencies

Notice that prDeep and LPR algorithms require Matconvnet (1.0-beta25), CUDA (10.0) and cuDNN.
