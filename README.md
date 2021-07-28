# ngVLA RFI Impact Simulator

A suite of scripts and apps to evaluate the effect of RFI on the ngVLA telescope (and the VLA)

## RFI_Impact_Calculator
Estimate the fraction of data loss due to RFI for the ngVLA using sample RFI characteristics, several RFI mitigation options and the effect of decorrelation. 
The purpose is to do a cost-benefit analysis of implementing end-to-end RFI mitigation solutions. 

Documentation is available as [ngVLA Memo #70](https://library.nrao.edu/public/memos/ngvla/NGVLA_70.pdf). A summary of the associated browser app and instructions on how to run it are located [here](./RFI_Impact_Calculator/README.md).

## RFI_Decorrelation_And_Attenuation
A tool that implements Eqn 19 of Thompson et al, 1982 to predict the degree of attenuation expected for RFI at different locations in the sky, for the VLA and ngVLA array geometries.

Documentation is available [here](./RFI_Decorrelation_And_Attenuation/README.md). 

## Single_Baseline_Correlator
A script to simulate/implement the inner workings of a single-baseline correlator. Show how complex visibilities are formed and relate to the input voltage streams per antenna, and also how delay and channel bandwidth result in the characteristic Sinc attenuation pattern due to decorrelation.  This is a teaching/learning tool (code from Shrishti Yadav). 

Documentation is available [here](./Single_Baseline_Correlator/README.md).





