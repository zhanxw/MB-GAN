# MB-GAN

Microbiome Simulation via Generative Adversarial Network

### Introduction

MB-GAN is a deep learning simulation framework for simulating realistic
microbiome data

### Information provided by this repository

##### Data folder

The folder “data” contains:

  - The real microbiome data used to generate MB-GAN samples analyzed in
    the manuscript

  - The csv files of the MB-GAN samples

  - The processed data used to generate figures and tables in the
    manuscript

##### R code

The folder “Rcode” includes the following R scripts:

  - `generate_NorTA.R`: Simulate data by Normal-To-Anything (NorTA)
  - `nMDS.R`: Perform non-metric multidimensional scaling
  - `ANOSIM.R`: Perform analysis of similarities in Table 1
  - `process_MBGAN.R`: Converts the csv files of MB-GAN samples to
    `Rdata`
  - Code for reproducing the figures in the manuscript

##### Python code

The root folder includes the following Python scripts:

  - `mbgan_train_demo.py`: codes to train a MB-GAN network
  - `mbgan_inference_casectrl.py`: codes to simulate new microbiome 
    abundances using trained model parameters
  
## Contact

Xiaowei Zhan <xiaowei.zhan@utsouthwestern.edu>, 
Quantitative Biomedical Research Center, 
Center for the Genetics of Host Defence,
Department of Population and Data Sciences,
UT Southwestern Medical Center, 
Dallas, TX 75390-8821
