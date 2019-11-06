import dessinercestgagner.src.stableRANK as sr
import numpy as np
from dessinercestgagner.src.Shape import Shape

inf = float("inf")
import matplotlib.pyplot as plt

plt.style.use('ggplot')
import _pickle as pickle
from ripser import ripser

list_of_objects = ['apple', 'bat', 'bell']  # , 'bird'
size_of_object = 21
sample_shapes = {(k, y): sr.euc_object(Shape(k, y).points) for k in list_of_objects for y in range(1, size_of_object)}
maxdim = 1
thresh = inf
coeff = 3
distance_matrix = False
do_cocycles = False
metric = "euclidean"
barcodes = {}

# sr.euc_object(Shape('beetle').points).barcode(maxdim, coeff, metric)

for k in sample_shapes.keys():
    barcodes[k] = sample_shapes[k].barcode(maxdim, thresh, coeff, distance_matrix, do_cocycles, metric)

temp = None
contours = {"H0": [[[0], [1]], ("area", inf, inf)]}
border = 1
linewidth = 1.0
i = 0
colours = ['r', 'g', 'b', 'm']
signature_array = []
normelized_sig_array = {}
signature_array_H0 = {}
signature_array_H1 = {}


def plot_signatures(signature_array, homology, title, color_index):
    plt.figure(figsize=(20, 20))
    for k in signature_array:
        k[homology].plot(border, colours[color_index], linewidth)
    plt.title("Signature of " + title + ", " + homology + " homology")
    plt.show()
    normelized_sig_array[title] = normalized_signaturs(signature_array)
    signature_array_H0[title], signature_array_H1[title] = get_H0_H1_signatures(signature_array)
    return


def get_H0_H1_signatures(signature_array):
    temp_H0 = []
    temp_H1 = []
    for k in signature_array:
        temp_H0.append(k['H0'])
        temp_H1.append(k['H1'])
    return temp_H0, temp_H1


def normalized_signaturs(signature_array):
    normelized_sig = []
    for k in signature_array:
        normelized_sig.append(k["H1"] * (k["H0"].invert()))
    return normelized_sig


for k in barcodes.keys():
    if temp != k[0]:
        temp_signature = barcodes[k].stable_rank_max(contours)
        signature_array.append(temp_signature)
    else:
        temp_signature = barcodes[k].stable_rank_max(contours)
        signature_array.append(temp_signature)
        if k[1] == 20:
            plot_signatures(signature_array, 'H0', k[0], i)
            signature_array = []
            i = i + 1
    temp = k[0]

def Average_pcf(pcf_array, title, color_index):
    a = 0
    for pcf in pcf_array:
        a += pcf
    a = a * (1 / len(pcf_array))
    plt.figure(figsize=(10, 10))
    border = 1
    linewidth = 1.0
    a.plot(border, colours[color_index], linewidth)
    plt.title("Average PCF-" + title)
    #plt.show()
    print("Did average PCF")
    return a


def correct_guesses(nom, data_Type, questionType, E, F, H):
    r = 0
    for pcf in nom[data_Type]:
        d_E = abs(pcf - E).integrate(0, 100)
        d_F = abs(pcf - F).integrate(0, 100)
        d_H = abs(pcf - H).integrate(0, 100)
        if questionType == 'apple':
            temp = d_E < d_F and d_E < d_H
            print(temp)
            if (temp):
                r += 1
        elif questionType == 'bat':
            temp = d_F < d_E and d_F < d_H
            print(temp)
            if (temp):
                r += 1
        elif questionType == 'bell':
            temp = d_H < d_F and d_H < d_E
            print(temp)
            if (temp):
                r += 1
        else:
            return 'Error'
    return r / len(nom[data_Type])


def check_H0_H1():
    apple = Average_pcf(normelized_sig_array['apple'], 'apple', 0)
    bat = Average_pcf(normelized_sig_array['bat'], 'bat', 1)
    bell = Average_pcf(normelized_sig_array['bell'], 'bell', 2)
    print("Normalized H1/H0")
    print(correct_guesses(normelized_sig_array, 'apple', 'apple', apple, bat, bell))


def check_H0():
    apple = Average_pcf(signature_array_H0['apple'], 'apple', 0)
    bat = Average_pcf(signature_array_H0['bat'], 'bat', 1)
    bell = Average_pcf(signature_array_H0['bell'], 'bell', 2)
    print("H0 Signatures")
    print(correct_guesses(signature_array_H0, 'apple', 'apple', apple, bat, bell))


def check_H1():
    apple = Average_pcf(signature_array_H1['apple'], 'apple', 0)
    bat = Average_pcf(signature_array_H1['bat'], 'bat', 1)
    bell = Average_pcf(signature_array_H1['bell'], 'bell', 2)
    print("H1 Signatures")
    print(correct_guesses(signature_array_H1, 'apple', 'apple', apple, bat, bell))


check_H0_H1()
check_H0()
check_H1()
