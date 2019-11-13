import dessinercestgagner.src.stableRANK as sr
import numpy as np
from dessinercestgagner.src.Shape import Shape
import operator
import os
import matplotlib.pyplot as plt
plt.style.use('ggplot')
inf = float("inf")
graph_dir = os.path.join(os.path.dirname(__file__), '../../output/')
plt.style.use('ggplot')

output = True
method = ''
list_of_objects = ['apple', 'bat', 'bell']  # , 'bird'
colours_shape = {'apple': 'g', 'bat': 'k', 'bell': 'y'}
size_of_object = 21
perturb_magn = 0.3

maxdim = 1
thresh = inf
coeff = 3
distance_matrix = False
do_cocycles = False
metric = "euclidean"
barcodes = {}


sample_shapes = {(k, y): sr.euc_object(Shape(k, y).perturb(method, perturb_magn).points) for k in list_of_objects
                                                                                         for y in range(1, size_of_object)}

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
signature_array_homology = {}

def normalized_signaturs(signature_array):
    normelized_sig = []
    for k in signature_array:
        normelized_sig.append(k["H1"] * (k["H0"].invert()))
    return normelized_sig

def get_H0_H1_signatures(signature_array):
    return [k['H0'] for k in signature_array], [k['H1'] for k in signature_array]

def plot_signatures(signature_array, homology, title, colour):
    plt.figure(figsize=(20, 20))
    for k in signature_array:
        k[homology].plot(border, colour, linewidth)
    plt.title("Signature of " + title + ", " + homology + " homology")
    plt.show()
    normelized_sig_array[title] = normalized_signaturs(signature_array)
    signature_array_H0[title], signature_array_H1[title] = get_H0_H1_signatures(signature_array)
    return

for k in barcodes.keys():
    if temp != k[0]:
        temp_signature = barcodes[k].stable_rank_max(contours)
        signature_array.append(temp_signature)
    else:
        temp_signature = barcodes[k].stable_rank_max(contours)
        signature_array.append(temp_signature)
        if k[1] == 20:
            plot_signatures(signature_array, 'H0', k[0], colours[i])
            signature_array = []
            i = i + 1
    temp = k[0]


def Average_pcf(pcf_array, title, color, plot = False):
    a = 0
    for pcf in pcf_array:
        a += pcf
    a = a * (1 / len(pcf_array))
    if plot == True:
        plt.figure(figsize=(10, 10))
        border = 1
        linewidth = 1.0
        a.plot(border, color, linewidth)
        plt.title("Average PCF-" + title)
        plt.show()
    print("Did average PCF")
    return a


def correct_guesses(nom, data_Type, types_mean, output):
    d = {}
    res = {fig_type: 0 for fig_type in types_mean.keys()}
    r = 0
    i = 0
    for pcf in nom[data_Type]:
        for type_name in types_mean.keys():
            d[type_name] = abs(pcf - types_mean[type_name]).integrate(0, 100)
        closest = min(d.items(), key=operator.itemgetter(1))[0]
        if output == True:
            plt.figure(figsize=(20,20))
            pcf.plot(1, color = 'b', linewidth = 1)
            for type_name in types_mean.keys():
                types_mean[type_name].plot(1, color = colours_shape[type_name], linewidth = 3)
            plt.title('real type:' + data_Type + ' guessed type:' + closest)
            plt.savefig(graph_dir + 'guess_figure/' + data_Type + '_' + str(i) + '_' + closest + '.pdf')
            plt.show()
            plt.close()
        if data_Type == closest:
            r = r+1
        else:
            print('error:' + closest + ' instead of ' + data_Type)
        res[closest] = res[closest] + 1
        i+=1
    return res


def check_H0_H1(output):
    means = {}
    for object in list_of_objects:
        means[object] = Average_pcf(normelized_sig_array[object], object, colours_shape[object])
    print("Normalized H1/H0")
    for object in list_of_objects:
        print(correct_guesses(normelized_sig_array, object, means, output))

def print_mean(shape_name, colour):
    shape_mean = Average_pcf(normelized_sig_array[shape_name], shape_name, colour)
    plt.figure(figsize=(20, 20))
    for k in normelized_sig_array[shape_name]:
        k.plot(border, colour, linewidth)
    shape_mean.plot(border, 'y', linewidth + 10)
    plt.title("Signature of " + shape_name + ", normalized homology")
    plt.savefig(graph_dir + 'signature_mean/' +shape_name + '_normalized_homology.pdf')
    plt.show()
    plt.close()
    return




def check_H0():
    apple = Average_pcf(signature_array_H0['apple'], 'apple', 0)
    bat = Average_pcf(signature_array_H0['bat'], 'bat', 1)
    bell = Average_pcf(signature_array_H0['bell'], 'bell', 2)
    print("H0 Signatures")
    print(correct_guesses(signature_array_H0, 'apple', 'apple', {'apple':apple, 'bat':bat, 'bell':bell}))


def check_H1():
    apple = Average_pcf(signature_array_H1['apple'], 'apple', 0)
    bat = Average_pcf(signature_array_H1['bat'], 'bat', 1)
    bell = Average_pcf(signature_array_H1['bell'], 'bell', 2)
    print("H1 Signatures")
    for object in list_of_objects:
        print(correct_guesses(signature_array_H1, object, object, {'apple':apple, 'bat':bat, 'bell':bell}))


check_H0_H1(output)

# Print mean and profile H1/H0
for object in list_of_objects:
    print_mean(object, 'b')


#check_H0()
#check_H1()
