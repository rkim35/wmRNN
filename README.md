# Working Memory Spiking RNN Model

## Overview

This repository provides the code for the model and analyses presented in [this paper](https://www.nature.com/articles/s41593-020-00753-w):

Kim R. & Sejnowski TJ. Strong inhibitory signaling underlies stable temporal dynamics and working memory in spiking neural networks. _Nature Neuroscience_ (2020).

Preprint available [here](https://www.biorxiv.org/content/10.1101/2020.02.11.944751v1).

## Recurrent Neural Network (RNN) model
For this work, we used the method that we previously developed (reported [here](https://www.pnas.org/content/116/45/22811)) to construct spiking RNNs to perform WM tasks. For more details and access to the code, please refer to [this GitHub repository](https://github.com/rkim35/spikeRNN). The trained RNNs analyzed in the study are also available [here](https://osf.io/md4wg/).

## Experimental data
We also analyzed a publicly available dataset collected and generously shared by Dr. Christos Constantinidis's Lab at Wake Forest School of Medicine. Please refer to [this website](http://crcns.org/data-sets/pfc/pfc-3/about-pfc-2) to get acecss and learn more about the dataset.

## Analysis
The code for computing neuronal timescales and analyzing both model and experimental data is implemented in MATLAB (tested in 2016a and 2016b).

### Model Analysis
All the scripts related to the model analysis are located in the `model` folder. The folder also contains a README file showing how to use the included scripts.


### Experimental Data Analysis
All the scripts related to the experimental data analysis are located in the `data` folder. The folder also contains a README file showing how to use the included scripts.

## Citation
If you use this repo for your research, please cite our work:

```
@article {Kim_2020,
  author = {Kim, Robert and Sejnowski, Terrence J.},
  title = {Strong inhibitory signaling underlies stable temporal dynamics and working memory in spiking neural networks},
  year = {2020},
  journal = {bioRxiv}
}

```
