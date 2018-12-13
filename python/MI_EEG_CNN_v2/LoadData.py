# _*_ coding: utf-8 _*_

'''
read_data.py for read data and labels from BCI competition data set IV 2b
the data set contains left-right hands motor imagery EEG signals form electrodes of C3, Cz, and C4
this file read data from *.mat (https://www.dropbox.com/sh/vkviczhx8ovy51s/AACFBeqSQNPBkX__kqvF-lEza?dl=1)

@self.dataset
|- @subject sub01, sub02, sub03...
    |- @session 0,1,2,3,4
        |- @data 'data'
        |- @label 'label'
'''

import os
import re
import numpy as np
from scipy.io import loadmat


class EEGDataSet:
    # para dataSetDir: the direction of EEG dataset
    def __init__(self, dataset_dir):
        self.subjectNum = 9
        self.sampleRate = 250
        self.dataset = dict()
        # read data file
        data_files = os.listdir(dataset_dir)
        for subject_file in data_files:
            if not re.search(".*\.mat", subject_file):
                continue

            info = re.findall('B(0[0-9])([TE])\.mat', subject_file)
            try:
                filename = dataset_dir + "\\" + subject_file
                print(filename)
                data = loadmat(filename)['data'][0]  # contains train sessions or evaluation sessions

                subject = "sub" + info[0][0]
                if subject not in self.dataset :
                    self.dataset[subject] = list()
                    for i in range(5):
                        self.dataset[subject].append([])
                if info[0][1] == 'T':
                    for i in [0, 1, 2]:
                        trials = data[i]['trial'][0][0]
                        signals = data[i]['X'][0][0]
                        labels = data[i]['y'][0][0]
                        self.dataset[subject][i] = dict()
                        self.dataset[subject][i]['data'] = list()
                        self.dataset[subject][i]['label'] = list()
                        label_idx = 0
                        for trial_start in trials:
                            # read trial's data from C3, Cz and C4 (8 seconds per trial)
                            self.dataset[subject][i]['data'].append(signals[trial_start[0]:trial_start[0] + 250 * 8, 0:3])
                            self.dataset[subject][i]['label'].append(labels[label_idx][0])
                            label_idx += 1
                elif info[0][1] == 'E':
                    for i in [0, 1]:
                        trials = data[i]['trial'][0][0]
                        signals = data[i]['X'][0][0]
                        labels = data[i]['y'][0][0]
                        self.dataset[subject][i+3] = dict()
                        self.dataset[subject][i+3]['data'] = list()
                        self.dataset[subject][i+3]['label'] = list()
                        label_idx = 0
                        for trial_start in trials:
                            # read trial's data from C3, Cz and C4
                            self.dataset[subject][i+3]['data'].append(signals[trial_start[0]:trial_start[0] + 250 * 8, 0:3])
                            self.dataset[subject][i+3]['label'].append(labels[label_idx][0])
                            label_idx += 1
            except Exception as e:
                print(e)
                raise "invalid labels file name"

    # print event type information of EEG data set
    @staticmethod
    def print_type_info():
        print("EEG data set event information and index:")
        print("%12s\t%10s\t%30s" % ("Event Type", "Type Index", "Description"))
        print("%12d\t%10d\t%30s" % (276, 0, "Idling EEG (eyes open)"))
        print("%12d\t%10d\t%30s" % (277, 1, "Idling EEG (eyes closed"))
        print("%12d\t%10d\t%30s" % (768, 2, "Start of a trial"))
        print("%12d\t%10d\t%30s" % (769, 3, "Cue onset left (class 1)"))
        print("%12d\t%10d\t%30s" % (770, 4, "Cue onset right (class 2)"))
        print("%12d\t%10d\t%30s" % (781, 5, "BCI feedback (continuous"))
        print("%12d\t%10d\t%30s" % (783, 6, "Cue unknown"))
        print("%12d\t%10d\t%30s" % (1023, 7, "Rejected trial"))
        print("%12d\t%10d\t%30s" % (1077, 8, "Horizontal eye movement"))
        print("%12d\t%10d\t%30s" % (1078, 9, "Vertical eye movement"))
        print("%12d\t%10d\t%30s" % (1079, 10, "Eye rotation"))
        print("%12d\t%10d\t%30s" % (1081, 11, "Eye blinks"))
        print("%12d\t%10d\t%30s" % (32766, 12, "Start of a new run"))


if __name__=="__main__":
    dataSetDir = r"D:\BCI\EEG dataset\bci competition IV 2b\Graz_mat"
    dataset = EEGDataSet(dataSetDir)
    print("Finish Load Data...")