## Getting data
First, download the trained RNNs from [here](https://osf.io/md4wg/). For the DMS RNN model, download and unzip "DMS\_OSF.zip". The folder contains `.mat` files where each file corresponds to one trained RNN for the DMS task. Each `.mat` file contains several important parameters including neuronal timescales (extracted using the method similar to the one developed by [Murray et al.](https://www.nature.com/articles/nn.3862)).

## Extracting all the parameters from the trained networks
The main function for extracting parameters from a group of trained RNNs is called `fnc_get_models.m`.
The data folder (DMS\_OSF) contains both poor and good performance DMS RNNs.


