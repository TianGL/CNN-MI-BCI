# Introduction

CNN-SAE program for MI-BCI classification. (Based on "A novel deep learning approach for classification of EEG motor imagery signals")

- Please fellow the lisense of rasmusbergpal when use this program.
- CNN-SAE(MI-BCI) is a matlab progam for the classification Motor Imagery EEG signals.
- CNN-SAE(MI-BCI) programed based on rasmusbergpal [DeepLearnToolbox](https://github.com/rasmusbergpalm/DeepLearnToolbox.git)

- The theory of this program based on [Tabar et al-2016-J Neural Eng](https://doi.org/10.1088/1741-2560/14/1/016003). Also, we have make some changes in it to improve the result.

- This performance of this program is based on [BCI Competioion IV dataset 2a](http://www.bbci.de/competition/iv/)(click here for more information).



To improve the performance, we investigate CNN with a form of input from short time Fourier transform (STFT) combining time, frequency and location information. Fisher discriminant analysis-type F-score based on band pass (BP) feature and power spectra density (PSD) feature are employed respectively to select the subjectoptimal frequency bands. In the experiments, typical frequency bands related to motor imagery EEG signals, subject-optimal frequency bands and Extension Frequency Bands are employed respectively as the frequency range of the input image of CNN.

The result of CNN can be found in the "python" file dictionary in "excel" files.

The python codes is programed based on tensorflow 1.6 with GPU acceleration.

The matlab code is smiliar with python codes but need a long time to getting the result.

# Notice

Before you run the program, you should first run ***biosig_installer.m*** in biosig4octmat-3.2.0. This is an opern-source toolbox [BioSig](http://biosig.sourceforge.net/) for loading the GDF file of the dataset.

I also think the result can be improved with carfully tunning.



# Log INFO

2018.12.13

> Upload CNN v2 python codes (/python/MI_EGG_CNN v2)

2018.10.23

@Geliang Tian (tglasd@163.com)

> Upload CNN codes (matbal and python)

2017.12.14 

@Geliang Tian (tglasd@163.com))

> Upload BioSig toolbox
