import pandas as pd
import numpy as np
from scipy.signal import butter, lfilter, stft
from scipy import interpolate
import math
from spectrum import pburg
from sklearn import preprocessing
import matplotlib.pyplot as plt

from read_data import get_data


# enrich data for each trail
def preprocess_signal( ori_data, start_time, slide_len, segment_len, num, sfreq ):
    processed_data = list()
    for i in range(num):
        left = int((start_time + i*slide_len)*sfreq)
        right = left + sfreq*segment_len
        # if need to be averaged
        data_foo = ori_data[left:right]
        data_foo = data_foo - np.mean(data_foo)

        processed_data.append(data_foo)
    return processed_data


# butterworth band pass filter design
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


# butterworth band pass filter
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


# combine all zhe data of all trials, session
def combine_processed_data(preprocessed_data, labels):
    combined_data = dict()
    combined_labels = dict()
    for subject in preprocessed_data:
        subject_data = pd.DataFrame()
        subject_labels = list()
        for session in preprocessed_data[subject]:
            subject_data = subject_data.append(preprocessed_data[subject][session])
            subject_labels.extend(labels[subject][session])
        subject_combined_data = pd.DataFrame()
        subject_combined_labels = list()
        labels_flag = True
        labels_idx = 0
        for channel in subject_data:
            channel_data = list()
            for trials_data in subject_data[channel]:
                for segment_data in trials_data:
                    channel_data.append(segment_data)
                    if labels_flag:
                        subject_combined_labels.append(subject_labels[labels_idx])
                labels_idx += 1
            subject_combined_data[channel] = channel_data
            labels_flag = False
        combined_data[subject] = subject_combined_data
        combined_labels[subject] = subject_combined_labels

    return combined_data, combined_labels


# subject optimal frequency bands selection methods based on Band Pass feature
# type = 0, BP features
# type = 1, AR features
def feature_band_selection(data, labels, sfreq, step=1, band_range = (0, 0), band_size=(0, 0, 0),
                              channel=('EEG:C3', 'EEG:Cz', 'EEG:C4'), features_type=0):
    # AR model parameters
    ar_order = 12
    nfft = 1000

    subject_optimal_frequency_bands = dict()

    # mu band selection
    for subject in data:
        if features_type == 1: # compute AR Model PSD based on burg algorithm
            ar_psd = dict()
            freq_flag = True
            for channel_name in channel:
                ar_psd[channel_name] = list()
            for idx in range(len(labels[subject])):
                for channel_name in channel:
                    x = data[subject][channel_name][idx]
                    p = pburg(x, order=ar_order, NFFT=nfft, sampling=sfreq, scale_by_freq=True)
                    if freq_flag:
                        ar_psd['frequency'] = np.array(p.frequencies())
                        freq_flag = False
                    ar_psd[channel_name].append(p.psd)

        f_score = list()
        optimal_band = list()
        for band in band_size:
            for num_windows in range(int((band_range[1]-band_range[0]-band)/step)):
                lowcut = band_range[0] + num_windows * step
                highcut = lowcut + band
                optimal_band.append((lowcut, highcut))
                if features_type == 1:
                    ar_freq = ar_psd['frequency']
                    psd_idx_start = np.where(ar_freq >= lowcut)[0][0]
                    psd_idx_end = np.where(ar_freq >= highcut)[0][0]
                left_features = list()
                right_features = list()
                for idx in range(len(labels[subject])):
                    features = list()
                    for channel_name in channel:
                        # BP features
                        if features_type == 0: # 5th butterworth filter
                            filtered_data = butter_bandpass_filter(data[subject][channel_name][idx], lowcut, highcut,
                                                                   sfreq, order=5)
                        elif features_type == 1: # AR model PSD
                             filtered_data = ar_psd[channel_name][idx][psd_idx_start : psd_idx_end]
                        else:
                            raise Exception("feature type wrong!\n band pass features: features_type=0\n "
                                            "AR PSD features: features_type=1")
                        features.append(math.log10(np.var(filtered_data)))
                    if labels[subject][idx] == 1:
                        left_features.append(features)
                    elif labels[subject][idx] == 2:
                        right_features.append(features)
                left_mean_val = np.mean(left_features, axis=0)
                right_mean_val = np.mean(right_features, axis=0)
                left_var = np.var(left_features, axis=0)
                right_var = np.var(right_features, axis=0)
                f_score.append(sum(np.square(left_mean_val-right_mean_val)) / sum(left_var+right_var))
        # get optimal frequency corresponding to max F-score
        subject_optimal_frequency_bands[subject] = optimal_band[f_score.index(max(f_score))]
        # pause = input("pause")
    return subject_optimal_frequency_bands

# find kth largest number in a 1-d array
def find_kth_largest(arr, k):
    k = k - 1
    lo = 0
    hi = len(arr) - 1
    while lo < hi:
        arr[lo], arr[int((lo+hi)/2)] = arr[int((lo+hi)/2)], arr[lo]
        left = lo; right = hi; pivot=arr[lo];
        while left < right:
            while left < right and arr[right] <= pivot:
                right = right - 1
            arr[left] = arr[right]
            while left < right and arr[left] >= pivot:
                left = left + 1
            arr[right] = arr[left]
        arr[left] = pivot
        if k <= left:
            hi = left - 1
        if k >= left:
            lo = left + 1
    res = arr[k]
    return res


# rescale data
# @percentage: percentage of value considered to be artifact
def recale(data, percentage):
    m, n = data.shape
    arr = data.flatten()
    min_val = np.min(arr)
    max_val = find_kth_largest(arr, int(m*n*percentage))
    for i in range(m):
        for j in range(n):
            if data[i][j] > max_val:
                data[i][j] = 1
            else:
                data[i][j] = (data[i][j] - min_val) / ( max_val - min_val)
    return data


# get the input data of CNN by STFT
def get_input_data(data, mu_band, beta_band, channel=('EEG:C3', 'EEG:Cz', 'EEG:C4')):
    # parameters of stft:
    wlen = 64  # length of the analysis Hamming window
    nfft = 512  # number of FFT points
    fs = 250  # sampling frequency, Hz
    hop = 14  # hop size

    input_data = list()
    num_segments = len(data[channel[0]])
    freq_flag = True
    for idx in range(num_segments):
        input_image = None
        for chn in channel:
            f, t, Fstft = stft(data[chn][idx], fs=fs, window='hamming', nperseg=wlen, noverlap=wlen-hop,
                               nfft=nfft, return_onesided=True, boundary=None, padded=False)
            if freq_flag:  # only need run one time
                mu_left = np.where(f >= mu_band[0])[0][0]
                mu_right = np.where(f >= mu_band[1])[0][0]
                beta_left = np.where(f >= beta_band[0])[0][0]
                beta_right = np.where(f >= beta_band[1])[0][0]
                freq_flag = False
            mu_feature_matrix = np.abs(Fstft[mu_left : mu_right])
            beta_feature_matrix = np.abs(Fstft[beta_left : beta_right])

            # beta band cubic interpolation
            beta_interp = interpolate.interp2d(t, f[beta_left : beta_right], beta_feature_matrix, kind='cubic')
            interNum = len(mu_feature_matrix)
            f_beta = np.arange(beta_band[0], beta_band[1], (beta_band[1]-beta_band[0])/(interNum))
            beta_feature_matrix = beta_interp(t, f_beta)
            # mu_feature_matrix = preprocessing.scale(np.array(mu_feature_matrix), axis=1)
            mu_feature_matrix = recale(mu_feature_matrix, 0.05)
            beta_feature_matrix = recale(beta_feature_matrix, 0.05)
            # mu_feature_matrix = preprocessing.scale(beta_feature_matrix, axis=1)
            # plt.pcolormesh(t, f_beta, beta_feature_matrix, vmin=0)
            # plt.show()
            # pause = input("pause")
            if input_image is None:
                input_image = np.append(mu_feature_matrix, beta_feature_matrix, axis=0)
            else:
                input_image = np.append(input_image, mu_feature_matrix, axis=0)
                input_image = np.append(input_image, beta_feature_matrix, axis=0)
        input_data.append(input_image)
    return input_data



# default run function
# @band_type = 0: band pass optimal frequency bands
# @band_type = 1: AR PSD optimal frequency bands
# @band_type = 2: extend frequency band
def run_sig_processing(data_src, labels_src, band_type):
    # parameters initialization
    start_time = 3
    time_slides = 0.2
    window_length = 2
    segments_num = 11

    data, labels, sfreq = get_data(data_src, labels_src)

    # execute
    preprocessed_data = dict()
    for subject in data:
        if subject not in preprocessed_data:
            preprocessed_data[subject] = dict()
        for session in data[subject]:
            df_trials_data = pd.DataFrame()
            for channel in data[subject][session]:
                session_data = data[subject][session][channel]
                trials_processed_data = list()
                for trial_data in session_data:
                    processed_data = preprocess_signal(trial_data, start_time, time_slides, window_length,
                                                       segments_num, sfreq)
                    trials_processed_data.append(processed_data)
                df_trials_data[channel] = trials_processed_data
            # print(df_trials_data)
            # pause = input("pause: ")
            preprocessed_data[subject][session] = df_trials_data

    if band_type == 0 or band_type == 1:
        combined_data, combined_labels = combine_processed_data(preprocessed_data, labels)
        mu_band = feature_band_selection(combined_data, combined_labels, sfreq, step=1, band_range=(4, 14),
                                          band_size=(4, 5, 6), features_type=band_type)
        beta_band = feature_band_selection(combined_data, combined_labels, sfreq, step=1, band_range=(16, 40),
                                          band_size=(4, 5, 6), features_type=band_type)
    else:
        mu_band = dict()
        beta_band = dict()
        for subject in preprocessed_data:
            mu_band[subject] = (4, 14)
            beta_band[subject] = (16, 32)

    # get input data of CNN, add to column of dataFrame form processed_data[subject][session]
    for subject in preprocessed_data:
        for session in preprocessed_data[subject]:
            preprocessed_data[subject][session]['input data'] \
                = preprocessed_data[subject][session].apply(get_input_data, axis=1,
                                                            mu_band=mu_band[subject], beta_band=beta_band[subject])

    return preprocessed_data, labels



if __name__ == "__main__":
    # for test 2
    data_src = r"D:\BCI\EEG dataset\bci competition IV 2b\temp_dataset"
    labels_src = r"D:\BCI\EEG dataset\bci competition IV 2b\temp_dataset\true_labels"
    data = run_sig_processing(data_src, labels_src, band_type=3)

    print("test:signal processing")