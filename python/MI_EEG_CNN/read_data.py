'''
read_data.py for read data and labels from BCI competition data set IV 2b
the data set contains left-right hands motor imagery EEG signals form electrodes of C3, Cz, and C4
'''

import numpy as np
import pandas as pd
import mne
import os
import re
from scipy.io import loadmat


# read raw EEG dataset of .gdf file
class RawEEGData:
    # rawData = list()
    # event = pd.DataFrame()
    # channel = list()
    # sample_freq = 0

    def __init__(self, file_name):
        raw = mne.io.read_raw_edf(file_name, preload=True, stim_channel=-1)
        event_type2idx = {276: 0, 277: 1, 768: 2, 769: 3, 770: 4, 781: 5, 783: 6, 1023: 7, 1077: 8, 1078: 9, 1079: 10,
                          1081: 11, 32766: 12}
        self.rawData = raw._data
        self.channel = raw._raw_extras[0]['ch_names']
        self.sample_freq = raw.info['sfreq']
        # self.rawData =
        self.event = pd.DataFrame({
            "length":raw._raw_extras[0]['events'][0],
            "position": raw._raw_extras[0]['events'][1],
            "event type": raw._raw_extras[0]['events'][2],
            "event index": [event_type2idx[event_type] for event_type in raw._raw_extras[0]['events'][2]],
            "duration": raw._raw_extras[0]['events'][4],
            "CHN": raw._raw_extras[0]['events'][3]
        })

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


# arrange data for training and test
def get_data(data_file_dir, labels_file_dir):
    RawEEGData.print_type_info()
    sfreq = 250  # sample frequency of dataset
    # read data file
    data = dict()
    data_files = os.listdir(data_file_dir)
    for data_file in data_files:
        if not re.search(".*\.gdf", data_file):
            continue
        
        info = re.findall('B0([0-9])0([0-9])[TE]\.gdf', data_file)
        try:
            subject = "subject" + info[0][0]
            session = "session" + info[0][1]
            filename = data_file_dir + "\\" + data_file
            print(filename)
            raw_eeg_data = RawEEGData(filename)
            trial_event = raw_eeg_data.event[raw_eeg_data.event['event index'] == 2]
            session_data = dict()
            for event, event_data in trial_event.iterrows():
                trial_data = raw_eeg_data.rawData[:, event_data['position']:event_data['position']+event_data['duration']]
                for idx in range(len(raw_eeg_data.channel)):
                    if raw_eeg_data.channel[idx] not in session_data:
                        session_data[raw_eeg_data.channel[idx]] = list()
                    session_data[raw_eeg_data.channel[idx]].append(trial_data[idx])
            if subject not in data:
                data[subject] = dict()
            data[subject][session] = session_data
        except Exception as e:
            print(e)
            raise ("invalid data file name")

    # read data file
    labels = dict()
    labels_files = os.listdir(labels_file_dir)
    for labels_file in labels_files:
        if not re.search(".*\.mat", labels_file):
            continue

        info = re.findall('B0([0-9])0([0-9])[TE]\.mat', labels_file)
        try:
            subject = "subject" + info[0][0]
            session = "session" + info[0][1]
            filename = labels_file_dir + "\\" + labels_file
            print(filename)
            session_label = loadmat(filename)
            session_label = session_label['classlabel'].astype(np.int8)
            if subject not in labels:
                labels[subject] = dict()
            labels[subject][session] = session_label
        except Exception as e:
            print(e)
            raise ("invalid labels file name")

    return data, labels, sfreq
    # print(data)
    # print(labels)


if __name__ == "__main__":
    # for test 1
    filename = r"D:\BCI\EEG dataset\bci competition IV 2b\dataset\B0101T.gdf"
    d = RawEEGData(filename)
    data = d.event[d.event['event index'] == 2]
    # print(0, "\n", data.shape)

    data_src = r"D:\BCI\EEG dataset\bci competition IV 2b\dataset"
    labels_src = r"D:\BCI\EEG dataset\bci competition IV 2b\dataset\true_labels"
    data, labels, sfreq = get_data(data_src, labels_src)
    print("test: read dataset")