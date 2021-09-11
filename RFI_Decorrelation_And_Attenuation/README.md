## Version 2 (Sept 2021)
### Measured SPFD from a LEO satellite and the effects of decorrelation and fringe-winding.

What : Calculation of received power from a LEO satellite and how it is modified by attenuation due to bandwidth decorrelation and fringe washing.

To view/run : [Demo Notebook](./Calculate_RFI_Power_And_Attenuation_Sept21.ipynb) and associated [Python script](./Predict_Decorr_Script_Sept21.py). 

Ack : Thanks to R.Selina & C.De-Pree for discussions.

## Version 1 (July 2021) 
### RFI Attenuation due to decorrelation and fringe-winding

What : A tool to predict the amount of expected attenuation of an RFI signal in an interferometer.

How : Implements Eqn 19 from Thompson et al, 1982 for user-definable phasecenter and RFI locations and a variety of interferometer array configurations. 

To view/run : [Demo Notebook](./Predict_Decorr_Demo.ipynb) and associated [Python script](./Predict_Decorr_Script.py).

Ack : Thanks to R.Selina & and S.Yadav for discussions.

![Example Screenshot for the ngVLA](./demo_ngvla_example.png)
