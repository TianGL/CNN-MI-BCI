# _*_ coding: utf-8 _*_
from LoadData import EEGDataSet
from SignalProcess import SigProcess, FeatureExtraction

import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import cohen_kappa_score
import pandas as pd
import tensorflow as tf
from tensorflow import keras


# get train data and labels for network
def arrange_data(data, labels, isChannelImage=False):
    output_data = list()
    output_labels = list()
    for idx in range(len(data)):
        for segment in data[idx]:
            if isChannelImage:
                output_data.append(reshape_channel_image(segment))
            else:
                output_data.append(np.expand_dims(segment, axis=2))
            if labels[idx] == 1:
                output_labels.append(0)
            else:
                output_labels.append(1)
    output_data = np.array(output_data)
    output_labels = np.array(output_labels)
    return output_data, output_labels


# reshape channel_data to be valid input image of cnn
def reshape_channel_image(data):
    data = np.array(data)
    data_shape = data.shape
    reshape_date = np.zeros((data_shape[1], data_shape[2], data_shape[0]))
    for i in range(data_shape[1]):
        for j in range(data_shape[2]):
            for k in range(data_shape[0]):
                reshape_date[i][j][k] = data[k][i][j]
    return reshape_date


# build model
def build_model(size_y, size_x, dim=1):
    # input layer
    img_input = keras.layers.Input(shape=(size_y, size_x, dim))

    # First convolution extracts 30 filters that are (size_y, 3)
    # Convolution is followed by max-pooling layer with a 1x10 window
    x = keras.layers.Conv2D(filters=30, kernel_size=(size_y, 3), activation='relu',
                            kernel_regularizer=keras.regularizers.l2(0))(img_input)
    x = keras.layers.MaxPooling2D(1, 10)(x)

    # Convolution is followed by max-pooling layer with a 1x10 window
    # x = keras.layers.Conv2D(filters=10, kernel_size=(1, 5), activation='relu',
    #                         kernel_regularizer=keras.regularizers.l2(0.05))(x)
    # x = keras.layers.MaxPooling2D(1, 3)(x)

    # Flatten feature map to a 1-dim tensor so we can add fully connected layers
    x = keras.layers.Flatten()(x)

    # Create a fully connected layer with ReLU activation and 512 hidden units
    x = keras.layers.Dense(64, activation='relu', kernel_regularizer=keras.regularizers.l2(0.01))(x)
    # x = keras.layers.MaxPooling1D(2)(x)

    # Create a fully connected layer with ReLU activation and 512 hidden units
    # x = keras.layers.Dense(32, activation='relu', kernel_regularizer=keras.regularizers.l2(0.01))(x)

    # Add a dropout rate of 0.5
    # x = keras.layers.Dropout(0.75)(x)

    # Create a fully connected layer with ReLU activation and 512 hidden units
    # x = keras.layers.Dense(256, activation='sigmoid')(x)

    # Add a dropout rate of 0.5
    # x = keras.layers.Dropout(0.75)(x)

    # Create output layer with a single node and sigmoid activation
    output = keras.layers.Dense(2, activation='sigmoid', kernel_regularizer=keras.regularizers.l2(0))(x)

    # Create model:
    # input = input feature map
    # output = input feature map + stacked convolution/maxpooling layers + fully
    # connected layer + sigmoid output layer
    model = keras.models.Model(img_input, output)
    model.summary()
    # model.compile(loss='categorical_crossentropy',
    #               optimizer=keras.optimizers.RMSprop(lr=0.001),
    #               metrics=['acc'])
    model.compile(loss='categorical_crossentropy',
                  optimizer='adam',
                  metrics=['acc'])
    return model


def trial_evaluate(model, data, labels, isChannelImage=False):
    acc = 0.0
    predict_label = list()
    for idx in range(len(data)):
        test_data, test_label = arrange_data(np.expand_dims(data[idx], axis=0), np.expand_dims(labels[idx], axis=0),
                                             isChannelImage=isChannelImage)
        test_label = keras.utils.to_categorical(test_label, num_classes=2)
        loss, accuracy = model.evaluate(test_data, test_label)
        if accuracy > 0.5:
            predict_label.append(labels[idx])
            acc += 1.0
        else:
            if labels[idx] == 1:
                predict_label.append(2)
            else:
                predict_label.append(1)
    acc = acc/len(data)
    return np.array(predict_label), acc


def run():
    dataset_dir = r"D:\BCI\EEG dataset\bci competition IV 2b\Graz_mat"
    dataset = EEGDataSet(dataset_dir)
    process = SigProcess(dataset, 0, 3)
    process.get_train_data(start=4.0, end=7.0, slide=0.1, window=2)
    process.data_standardization()
    print(process.get_data_info())

    classification_acc = pd.DataFrame()
    for subject in process.data:
        # k-fold cross-validation
        kf = KFold(n_splits=10, shuffle=True)
        subject_acc = list()
        kappa_acc = list()
        # feature extraction
        fe = FeatureExtraction(process.data[subject]['train_data'], process.data[subject]['train_label'])
        fe.stft_process()
        fe.stft_standard()
        print(subject, 'Finish signal processing...')

        input_data = np.array(fe.data['stft_channel'])
        target_labels = np.array(fe.data['train_label'])

        # 10 fold cross-validation
        count = 0
        for train_index, test_index in kf.split(target_labels):
            count += 1
            train_data, train_labels = arrange_data(input_data[train_index], target_labels[train_index], isChannelImage=True)
            test_data, test_labels = arrange_data(input_data[test_index], target_labels[test_index], isChannelImage=True)

            size_y, size_x, dim = train_data[0].shape[0:3]

            print(train_data.shape)
            print(train_labels.shape)
            # train_data_size = train_data.shape[0]
            # test_data_size = test_data.shape[0]

            train_labels = keras.utils.to_categorical(train_labels, num_classes=2)
            test_labels = keras.utils.to_categorical(test_labels, num_classes=2)

            # build model
            model = build_model(size_y, size_x, dim=dim)

            print('Training ------------')
            # train the model
            model.fit(train_data, train_labels, verbose=1, epochs=300, batch_size=40)

            print('\nTesting ------------')
            # Evaluate the model with the metrics we defined earlier
            loss, accuracy = model.evaluate(test_data, test_labels)
            predict_label, trial_acc = trial_evaluate(model, input_data[test_index], target_labels[test_index], isChannelImage=True)
            print(count, subject)
            print('test loss: ', loss)
            print('test accuracy: ', accuracy)
            true_label = target_labels[test_index]
            kappa = cohen_kappa_score(predict_label, true_label)
            subject_acc.append(trial_acc)
            kappa_acc.append(kappa)
            print("trial acc: ", subject_acc)
            print("kappa_acc： ", kappa_acc)
            # input("Finish " + str(count) + " ...")
        print("mean acc: ", np.mean(subject_acc))
        # input("Finish " + str(subject) + " test, input to continue...")
        classification_acc[subject] = subject_acc
        classification_acc[subject + "_kappa_acc"] = kappa_acc
    return classification_acc

if __name__ == "__main__":
    res = run()
    print(res)
    res.to_csv("new_acc_64node_300epos.csv", encoding="utf-8")
    print("cnn classification...")