from __future__ import print_function
import sys, os, math, time
import numpy as np
import keras
from keras.datasets import mnist
from keras.layers import Dense, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.models import Sequential
from ont_fast5_api.fast5_interface import get_fast5_file
from sklearn.model_selection import train_test_split
import matplotlib.pylab as plt

batch_size = 128
num_classes = 10
epochs = 10

# things to change paramters for (alpha, lambda)
# build validation set by splitting train files using train_test_split

class logisticRegression():

    def __init__(self):
        self.count = 1; self.numOfFiles = 0
        self.x = []; self.y = []
        self.x_train = []; self.x_test = []
        self.y_train = []; self.y_test = []
        self.input_shape = ()

    #import data from Fast and multiFast5
    def importMultiFast5(self):

        # read number gives total number of reads
        # sub read number gives total number of matricies created from reads
        # file number counts the number of multiFast5 files read
        subReadNumber = 0; readNumber = 0; y_rng = 0; fileNumber = 1; n = 1024; tempFolder = []; y_arrays2 = []; numPerFile = []

        seconds = time.time()
        local_time = time.ctime(seconds)
        print("Time Before Getting Reads =", local_time)

        for file in os.listdir(sys.argv[1]):
            if file.endswith(".fast5"): # only read files that are in fast5 format
                path = (sys.argv[1] + '/' + file) # format file to be read
                with get_fast5_file(path, mode="r") as fileH: # open multifast5 file
                    for read in fileH.get_reads():
                        raw_data = read.get_raw_data() # get raw data
                        filter_raw = raw_data.tolist() # turn array into list
                        y = min(filter_raw) # get min value of list for curve
                        filter_raw = [i - y for i in filter_raw] # set curve of list to base 0
                        y_rng = n # set Y value for matrix
                        tempFolder = [filter_raw[i * n:(i + 1) * n] for i in range((len(filter_raw) + n - 1) // n )] # devide reads into x number of lists = len 1024
                        if len(tempFolder[:-1]) != n: # get rid of last read if it is not exactly len 1024
                            tempFolder = tempFolder[:-1]
                        for list in tempFolder: # add reads to list we are considering
                            self.x.append(list)
                            subReadNumber += 1
                            y_arrays2.append(fileNumber) # tells Keras Model identity of matrix
                        readNumber += 1
                    fileNumber += 1
                    numPerFile.append(subReadNumber)
                continue
            else:
                continue

        seconds = time.time()
        local_time = time.ctime(seconds)
        print("Time After Getting/Formatting Reads =", local_time)

        # create 32x32 arrays
        dn = 32; storageSpace = []
        for j in self.x:
            l4n = [j[i * dn:(i + 1) * dn] for i in range((len(j) + dn - 1) // dn )] # turn lists into lists of lists 32x32
            storageSpace.append(l4n)
        x_arrays = np.array(storageSpace) # turn into matrix

        # set up image shape
        img_x = dn; img_y = dn
        self.input_shape = (img_x, img_y, 1)

        # get time array (Matrix pic)
        time_array = np.array(range(len(x_arrays)))

        #get Y data (X = 1 array)
        y_arrays = np.ones(len(x_arrays), dtype=int)
        y_arraysIter = np.array(y_arrays2)

        seconds = time.time()
        local_time = time.ctime(seconds)
        print("Time After Formating Data =", local_time)

        # devide data into train and test files
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(x_arrays, y_arraysIter, test_size=.2, random_state=42)

        seconds = time.time()
        local_time = time.ctime(seconds)
        print("Time After All Data Prep/Splits =", local_time)

        self.x_train = np.array(self.x_train); self.y_train = np.array(self.y_train)
        self.x_test = np.array(self.x_test); self.y_test = np.array(self.y_test)

        # reshape image to work with Keras
        self.x_train = self.x_train.reshape(self.x_train.shape[0], img_x, img_y, 1)
        self.x_test = self.x_test.reshape(self.x_test.shape[0], img_x, img_y, 1)


        seconds = time.time()
        local_time = time.ctime(seconds)
        print("Time After Converting Data to Arrays =", local_time)

        return logisticRegression()


lr = logisticRegression()
get = lr.importMultiFast5()

x_train = lr.x_train; x_test = lr.x_test
y_train = lr.y_train; y_test = lr.y_test
input_shape = lr.input_shape

# format X data
x_train = x_train.astype('float32')
x_test = x_test.astype('float32')

# format Y data
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time Before Modeling =", local_time)

# start Keras Modeling
model = Sequential()

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After Sequential Model =", local_time)

model.add(Conv2D(32, kernel_size=(5, 5), strides=(1, 1),
                 activation='relu',
                 input_shape=input_shape)) # use kernal size 5x5 maybe?

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After 2D Conversion =", local_time)

model.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2)))

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After 2D Pooling =", local_time)

model.add(Conv2D(64, (5, 5), activation='relu'))

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After 2nd 2D Conversion =", local_time)

model.add(MaxPooling2D(pool_size=(2, 2)))

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After 2nd 2D Pooling =", local_time)

model.add(Flatten())

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After Flatten =", local_time)

model.add(Dense(1000, activation='relu'))

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After Relu =", local_time)

model.add(Dense(num_classes, activation='softmax'))

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After Dense  =", local_time)

model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])

#time stamp
seconds = time.time()
local_time = time.ctime(seconds)
print("Time After Compile =", local_time)

class AccuracyHistory(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.acc = []

    def on_epoch_end(self, batch, logs={}):
        self.acc.append(logs.get('acc'))

history = AccuracyHistory()


model.fit(x_train, y_train,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test),
          callbacks=[history])
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
"""
plt.plot(range(1, 11), history.acc)
plt.xlabel('Epochs'); plt.xlim([0, epochs])
plt.ylabel('Accuracy'); plt.ylim([0, 1])
plt.show()
"""
