#!/usr/bin/env python3
import matplotlib.pyplot as mplplot, matplotlib.image as mpimg
import numpy as np, statistics as stat
import sys, os, math, time
from decimal import Decimal
from fast5_research import Fast5
from fast5_research import BulkFast5
import argparse as ap
from collections import Counter

'''
This program takes a folder fast5 and multifast5 files from command line and creates a possion
distribution graph that models the read data inside the files. Defualt is fast5 but user can
specify multifast5 by entering the command -m5. The program will also print the mean # of
occurences, the time window, and number of time windows.
'''


###############################################################################

class modelNanopore():

    def __init__(self):
        self.eventChannel = []; self.eventTime = []; self.totalTimes = []
        self.multiTimeStart = []; self.multiTimesEnd = []
        self.totalTime = 1; self.numOfFiles = 0
        self.x_points = []; self.y_points = []


###############################################################################

###############################################################################

    #import data (fast5)
    def importFast5(self):
        startTime = 0
        readlen = []

        for file in os.listdir(sys.argv[1]):
            if file.endswith(".fast5"):
                self.numOfFiles += 1
                path = (sys.argv[1] + '/' + file)
                with Fast5(path) as fileH:
                    tempDict = fileH.summary(scale=False)
                    self.eventChannel.append(tempDict.get("channel", ""))
                    self.eventTime.append(tempDict.get("start_time", ""))
                    self.totalTimes.append(tempDict.get("strand_duration", ""))
                    reads = fileH.get_reads(group=False, raw=True, read_numbers=None) #gives me raw voltage levels in Pico Amps
                    for read in reads:
                        readlen.append(len(read))
                    # print(tempDict.get("strand_duration", ""))
                continue
            else:
                continue

        self.totalTime = max(self.eventTime) + max(self.totalTimes)


###############################################################################

###############################################################################

    #import data (multiFast5) (in progress / not functional)
    def importMultiFast5(self):
        startTime = 0
        readlen = []

        for file in os.listdir(sys.argv[1]):
            if file.endswith(".fast5"):
                path = (sys.argv[1] + '/' + file)
                with Fast5(path) as fileH:
                    reads = fileH.get_reads(self, channel, transitions=False, penultimate_class=True)
                    for read in reads:
                        readlen.append(len(read))
                        self.eventChannel.append(read.get("channel", ""))
                        self.eventTime.append(read.get("event_index_start", ""))
                        self.totalTimes.append(read.get("event_index_end", ""))
                        self.multiTimeStart.append(read.get("read_start", ""))
                        self.multiTimesEnd.append(read.get("read_length", ""))
                        self.numOfFiles += 1
                    waveformTimings = fileH.get_waveform_timings(self)
                continue
            else:
                continue

        self.totalTime = max(self.eventTime) + max(self.totalTimes)


###############################################################################

###############################################################################

    # model fast5 data to possion distribution
    def modelPossion(self):

        x_list = []; y_list=[]

        # get (time windows)
        timeIntervalLen = round(self.numOfFiles * 3.9) # calc # of time windows
        temp_x_list = [self.totalTime for i in range(0,timeIntervalLen)]
        count = 1 / timeIntervalLen
        # subdivie my totalTime into x number of equal sized 'windows' (make this user param later)
        for i in temp_x_list:
            x_value = round(i * count)
            x_list.append(x_value)
            count += 1 / timeIntervalLen
        time_windows = x_list
        time_window_size = time_windows[0] # get size of time window

        # get mean
        mean = sum(self.eventTime) / len(self.eventTime)

        # get variance data (# of instances in the time window)
        temp_eventTime = self.eventTime; temp_y_list = []

        # get number of instances in each time window of event times
        for i in time_windows:
            count = 0
            for j in temp_eventTime:
                if j <= i:
                    temp_eventTime.remove(j)
                    count += 1
                else:
                    continue
            temp_y_list.append(count)

        # calculate # of instances in each window
        instances = temp_y_list # occurences
        meanInstances = sum(instances) / len(instances) # get mean of instances per time window
        instancesNoRepeats = list(dict.fromkeys(instances))
        instancesNoRepeats = sorted(instancesNoRepeats)

        # ugh y
        graphFM = [[x,temp_y_list.count(x)] for x in set(temp_y_list)]
        print(graphFM)
        graphFMextendo = []; graphFMratio = []
        for i in graphFM:
            graphFMextendo.append(i[-1])
        for i in graphFMextendo:
            graphFMratio.append(i / sum(graphFMextendo))
        PMFatZero = graphFMratio[0]
        print(graphFMratio)
        print ("lambda given k is Zero", PMFatZero)
        newMean = np.log(PMFatZero) * -1
        print ("Poisson Average Number of Occurences", newMean)
        print ("Actual Average?", meanInstances)
        print ("Number of Partitions", timeIntervalLen)

        # get Poisson distribution
        for j in instancesNoRepeats:
            Upper =  ( meanInstances ** j ) * ( math.exp(-1 * meanInstances) )
            Lower = math.factorial(j)
            PMF =  Upper / Lower  # PMF = ((lambda**k)*(e**-lambda))/(math.factorial(k))
            y_list.append(PMF)

        # self.x_points = instancesNoRepeats
        # self.y_points = y_list
        self.x_points = instancesNoRepeats
        self.y_points = graphFMratio

###############################################################################

###############################################################################

    # Plot Possion Data
    def graphData(self):

        # set style and create panel
        mplplot.style.use('DBlab.mplstyle')
        fig_width=5; fig_height=3
        mplplot.figure(figsize=(fig_width,fig_height))

        #set panel
        panel1 = mplplot.axes([.1,.15,.85,.8],frameon=True)
        # panel1.set_xlim(0, totalTime); # panel1.set_ylim(0,1)
        panel1.set_xlabel('Number of Occurences'); panel1.set_ylabel('P(K|\u03BB)')

        # plot data
        panel1.plot(self.x_points,self.y_points, color='black', linewidth=.5, zorder=0)
        panel1.plot(self.x_points,self.y_points, color='red', marker='o',markersize=1.5, linewidth=0, zorder=0)

        # getting rid of the ticks and labels on graph
        panel1.tick_params(axis='both',which='both',\
					bottom='on', labelbottom='on',\
					left='on', labelleft='on',\
					right='off', labelright='off',\
					top='off', labeltop='off')

        # save image to file
        mplplot.savefig('modelNanoporeReadsPossionGraph2080.png', dpi=600)


###############################################################################

###############################################################################

    # take user commands (in progress / use sys.argv for now)
    def userParams(self, inpt=None):
        getParams = ap.ArgumentParser(prefix_chars='-')
        #    getParams.add_argument('-5', '--Fast5', default=False, help='file type (multifast5 used post January 2018), set to true to read multiFast5, otherwise program will only take first read in the file')
        #    getParams.add_argument('-m5', '--multiFast5', default=False, help='file type (multifast5 used post January 2018), set to true to read multiFast5, otherwise program will only take first read in the file')
        getParams.add_argument('-u', '--mean', default=False, help='If set to True, print variance')
        getParams.add_argument('-tw', '--time window', default=False, help='If set to True, print time window size and # of time windows')
        #    getParams.add_argument('-iF', '--inputFile', default='', help='Enter folder of files to run fast5 files')
        #    getParams.add_argument('-oF', '--outFile', default='PossionModel.png', help='Enter in format, filename.png to save model to a specific file name')

        if inpt is None:
            args = getParams.parse_args()
        else:
            args = getParams.parse_args(inpt)


###############################################################################

###############################################################################

def main():
    mn = modelNanopore()
#    params = mn.userParams()
    fast5Folder = sys.argv[1]; # outputFile = sys.argv[2]
    if fast5Folder == '':
        print('User must specify a file path to read')
        sys.exit()

#    if params.args.multifast5 == True:
#        multifast5 = importMultiFast5()
#    else:

#    multifast5 = mn.importMultiFast5()
    fast5 = mn.importFast5()
    modelData = mn.modelPossion()
    graphData = mn.graphData()


main()
