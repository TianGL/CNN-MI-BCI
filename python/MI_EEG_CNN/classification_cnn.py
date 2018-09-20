import numpy as np
from sklearn.model_selection import KFold
import pickle
import tensorflow as tf
from tensorflow import keras
import pandas as pd


from signalProcessing import run_sig_processing

# get train data and labels for each segment
def arrange_data(data, labels):
    output_data = list()
    output_labels = list()
    for idx in range(len(data)):
        for segment in data[idx]:
            output_data.append(np.expand_dims(segment, axis=2))
            if labels[idx][0] == 1:
                output_labels.append(0)
            else:
                output_labels.append(1)
    output_data = np.array(output_data)
    output_labels = np.array(output_labels)
    return output_data, output_labels


# build model
def build_model(size_y, size_x):
    # input layer
    img_input = keras.layers.Input(shape=(size_y, size_x, 1))

    # First convolution extracts 30 filters that are (size_y, 3)
    # Convolution is followed by max-pooling layer with a 1x10 window
    x = keras.layers.Conv2D(filters=30, kernel_size=(size_y, 3), activation='relu')(img_input)
    x = keras.layers.MaxPooling2D(1, 10)(x)

    # Flatten feature map to a 1-dim tensor so we can add fully connected layers
    x = keras.layers.Flatten()(x)

    # Create a fully connected layer with ReLU activation and 512 hidden units
    # x = keras.layers.Dense(512, activation='relu')(x)

    # Add a dropout rate of 0.5
    # x = keras.layers.Dropout(0.75)(x)

    # Create output layer with a single node and sigmoid activation
    output = keras.layers.Dense(2, activation='sigmoid')(x)

    # Create model:
    # input = input feature map
    # output = input feature map + stacked convolution/maxpooling layers + fully
    # connected layer + sigmoid output layer
    model = keras.models.Model(img_input, output)
    model.summary()
    model.compile(loss='categorical_crossentropy',
                  optimizer=keras.optimizers.RMSprop(lr=0.001),
                  metrics=['acc'])
    return model


# evaluated trial to trial performance
def trial_evaluate(model, data, labels):
    acc = 0.0
    for idx in range(len(data)):
        test_data, test_label = arrange_data(np.expand_dims(data[idx], axis=0), np.expand_dims(labels[idx], axis=0))
        test_label = keras.utils.to_categorical(test_label, num_classes=2)
        loss, accuracy = model.evaluate(test_data, test_label)
        if accuracy > 0.5:
            acc += 1.0
    acc = acc/len(data)
    return acc


# run classification
def run_classification(data, labels, session=(1, 2, 3, 4, 5)):
    kf = KFold(n_splits=10, shuffle=True)
    classification_acc = pd.DataFrame()
    for subject in data:
        # if subject == 'subject1': continue
        subject_acc = list()
        input_data = list()
        target_labels = list()
        # combine trials data of target session
        [input_data.extend(data[subject]["session" + str(idx)]['input data']) for idx in session]
        [target_labels.extend(labels[subject]["session" + str(idx)]) for idx in session]
        input_data = np.array(input_data)
        target_labels = np.array(target_labels)

        # 10 fold cross-validation
        count = 0
        for train_index, test_index in kf.split(input_data):
            count += 1
            train_data, train_labels = arrange_data(input_data[train_index], target_labels[train_index])
            test_data, test_labels = arrange_data(input_data[test_index], target_labels[test_index])

            size_y, size_x = train_data[0].shape[0:2]

            print(train_data.shape)
            # train_data_size = train_data.shape[0]
            # test_data_size = test_data.shape[0]

            train_labels = keras.utils.to_categorical(train_labels, num_classes=2)
            test_labels = keras.utils.to_categorical(test_labels, num_classes=2)


            # build model
            model = build_model(size_y, size_x)

            print('Training ------------')
            # train the model
            model.fit(train_data, train_labels, epochs=300, batch_size=40)

            print('\nTesting ------------')
            # Evaluate the model with the metrics we defined earlier
            loss, accuracy = model.evaluate(test_data, test_labels)

            trial_acc = trial_evaluate(model, input_data[test_index], target_labels[test_index])
            print(count, subject)
            print('test loss: ', loss)
            print('test accuracy: ', accuracy)
            print('trial to trial accuracy: ', trial_acc)
            subject_acc.append(trial_acc)
        classification_acc[subject] = subject_acc
    return classification_acc




if __name__ == '__main__':
    # '''
    data_src = r"D:\BCI\EEG dataset\bci competition IV 2b\dataset"
    labels_src = r"D:\BCI\EEG dataset\bci competition IV 2b\dataset\true_labels"

    #
    # band_type: 0: band pass feature, 1: AR PSD feature, 2: extend band
    #
    data, labels = run_sig_processing(data_src, labels_src, band_type=2)

    # Saving the data and labels:
    # with open('temp_data.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    #     pickle.dump([data, labels], f)
    # '''
    # Getting back the data and labels:
    # with open('temp_data.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    #     data, labels = pickle.load(f)

    res = run_classification(data, labels)
    print(res)
    res.to_csv("BP_acc.csv", encoding="utf-8")
    print("cnn classification")