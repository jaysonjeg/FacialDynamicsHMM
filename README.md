# FacialDynamicsHMM
This repository contains the Matlab-based pipeline to analyse facial affect dynamics with time-frequency representations and hidden Markov model. The code will run this pipeline on the openly available DISFA dataset, and reproduce the healthy control results in the paper "Quantifying dynamic facial expressions under naturalistic conditions".

To reproduce the results, run the following in order:
1) DISFA_Step1_vid2csv.m. Calls OpenFace to generate facial AU time series from DISFA dataset.
2) DISFA_Step2_csv2cubes.m. Puts facial AU time series into nice Matlab arrays
3) DISFA_Step3_cubes2cwt2hmm.m. Runs the pipeline, starting from Matlab arrays containing AU time series. Generates Figures 1 and 2.

If you don't want to reproduce the results on DISFA but rather want to apply it on your own data, first run OpenFace on your data to produce action unit time series. Then make sure your data is in a matlab array called 'cube' which is ntimepoints x nAUs x nSubjects. Using our function get_cwt.m, type in Matlab: abs(get_cwt(cube,YourWebcamFrameRate,true)) to get wavelet transform representation of your data. Then proceed from line 265 in DISFA_Step3_cubes2cwt2hmm onwards.

System:
I used a Windows 10 operating system, 64-bit, 6 cores, with 64GB RAM. I used Matlab R2021a.

Software dependencies:
- DISFA dataset (http://mohammadmahoor.com/disfa/)
- OpenFace version 2.2.0, (https://github.com/TadasBaltrusaitis/OpenFace), commit ad1b3cc45ca05c762b87356c18ad030fcf0f746e, downloaded 8/2/2020 
- HMM-MAR v1.0, (https://github.com/OHBA-analysis/HMM-MAR), commit 710dd64108b578e980e5007d908303e76de9ca57, downloaded 30/6/21
- Matlab 'mafdr': (https://au.mathworks.com/help/bioinfo/ref/mafdr.html)
- Matlab 'panel', optionally to generate figures, (https://au.mathworks.com/matlabcentral/fileexchange/20003-panel)
- Some Matlab panels require 'ds2nfu', (https://au.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion)
- FACSHuman, optionally to make avatar faces in Figure 2, (https://github.com/montybot/FACSHuman)

The HMM-MAR code was modified slightly. These are the changes.
- In initGamma_random.m, changed line 17 to rng(1), (for reproducibility)
- In hmmmar_init.m, lines 109-113, deleted last parameter options.priorOFFvsON (to prevent an error)
- In padGamma.m, line 27, deleted ‘offset=sum(d)’ and ‘Tshifted=T-offset’. Add ‘Tshifted=T’. (to prevent error with padGamma).

My fork of HMM-MAR at https://github.com/jaysonjeg/HMM-MAR/tree/branch1, already has these changes

Citation:
If you use the scripts, please cite the following:
https://elifesciences.org/articles/79581
