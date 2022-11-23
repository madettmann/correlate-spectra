import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

def get_xy(filename, type='sim'):
    x = []
    y = []
    largest = 0
    if type == 'sim':
        dftFile = open(filename)
        i = 0
        for line in dftFile:
            if line[1] != '#':
                newData = line.split(',')
                newData = [float(val) for val in newData]
                x.append(newData[0])
                y.append(newData[1]+newData[2])
                if x[-1] > 20 and y[-1] > largest:
                    largest = y[-1]
            i = i + 1
        dftFile.close()
        return x, y, max(y)
    elif type == 'exp':
        expFile = open(filename)
        i = 0
        flag = 0
        isMev=True
        for line in expFile:
            if i > 0 and line != '\n':
                data = line.split(',')
                if len(data) > 2 and flag == 2:
                    if float(data[0]) < -16:
                        isMev=False
                    if float(data[0]) > 8:
                        x.append(float(data[0]))
                        y.append(float(data[1]))	
                        if x[-1] > 100 and y[-1] > largest:
                            largest = y[-1]
                elif len(data) == 1:
                    flag = flag + 1
            i = i + 1
        expFile.close()
        if isMev:
            x = [val * 8.06554 for val in x]
        return x, y, max(y)
    elif type == 'cleaned':
        file = open(filename)
        x = []
        y = []
        for line in file:
            x.append(float(line.split(' ')[0]))
            y.append(float(line.split(' ')[1]))
        return x, y, sum(y)/len(y)
def getMapping(ref_x, ref_y, sample_x, sample_y):
    # Finds mapping from sample to reference and returns
    if len(ref_x) ==  len(sample_x):
        return sample_x, sample_y
    else:
        new_x = []
        new_y = []
        for rx, ry in zip(ref_x, ref_y):
            nearest_x = 0
            nearest_y = 0
            dist = 999999
            for sx, sy in zip(sample_x, sample_y):
                if abs(sx-rx) < dist:
                    nearest_x = sx
                    nearest_y = sy
                    dist = abs(sx-rx)
            new_x.append(nearest_x)
            new_y.append(nearest_y)
        return new_x, new_y

def normalize(ref_max, sample_y, sample_max):
    return [y * ref_max/sample_max for y in sample_y]

def correlate(file_1, file_type_1, file_2, file_type_2):
    # Computes the correlation between two spectra stored in files file_1 and file_2.
    # file_type_ 1 and 2 are the types of file. Use sim for simulated files and exp for experimental files
    x_1, y_1, max_1 = get_xy(file_1, file_type_1)
    x_2, y_2, max_2 = get_xy(file_2, file_type_2)

    # Determine which set has more points -> longer must be mapped to smaller
    largest_is_1 = len(x_1) > len(x_2)
    if largest_is_1:
        big_x = x_1
        big_y = y_1
        big_max = max_1
        small_x = x_2
        small_y = y_2
        small_max = max_2
    else:
        big_x = x_2
        big_y = y_2
        big_max = max_2
        small_x = x_1
        small_y = y_1
        small_max = max_1
    
    # Normalize Data
    big_y = normalize(small_max, big_y, small_max)

    # Map big x and y values onto small 
    new_big_x, new_big_y = getMapping(small_x, small_y, big_x, big_y)

    df = pd.DataFrame()
    df['small'] = small_y
    df['big'] = new_big_y
    return df.corr().iloc[0][1]

def get_df(file_list, file_type_list, normalize_data=True, return_x=False):
    x_vals = []
    y_vals = []
    max_vals = []
    for file_name, file_type in zip(file_list, file_type_list):
        get_xy(file_name, file_type)
        x, y, maximum = get_xy(file_name, file_type)
        x_vals.append(x)
        y_vals.append(y)
        max_vals.append(maximum)

    min_index = 0
    min_len = np.inf
    for i, y in enumerate(y_vals):
        if len(y) < min_len:
            min_len = len(y)
            min_index = i
    
    df = pd.DataFrame()
    for i in range(len(x_vals)):
        # Normalize all larger data to smaller
        if normalize_data:
            if i != min_index:
                y_vals[i] = normalize(max_vals[min_index], y_vals[i], max_vals[i])
                # Map big x and y values onto small 
                x_vals[i], y_vals[i] = getMapping(x_vals[min_index], y_vals[min_index], x_vals[i], y_vals[i])
        new_name= file_list[i].split('.')[0].replace('/','_').lower()
        df[new_name] = y_vals[i]
    if return_x:
        df['x'] = x_vals[min_index]
    return df

def correlate(file_list, file_type_list, normalize_data=True, return_df=False, return_x=False):
    # Computes the correlation between multiple spectra stored in files.
    # file_type_lists is a list of the types of files. Use sim for simulated files and exp for experimental files
    x_vals = []
    y_vals = []
    max_vals = []
    for file_name, file_type in zip(file_list, file_type_list):
        get_xy(file_name, file_type)
        x, y, maximum = get_xy(file_name, file_type)
        x_vals.append(x)
        y_vals.append(y)
        max_vals.append(maximum)
        

    # Determine which set has more points -> longer must be mapped to smaller
    min_index = 0
    min_len = np.inf
    for i, y in enumerate(y_vals):
        if len(y) < min_len:
            min_len = len(y)
            min_index = i
    
    df = pd.DataFrame()
    for i in range(len(x_vals)):
        # Normalize all larger data to smaller
        if normalize_data:
            if i != min_index:
                y_vals[i] = normalize(max_vals[min_index], y_vals[i], max_vals[i])
                # Map big x and y values onto small 
                x_vals[i], y_vals[i] = getMapping(x_vals[min_index], y_vals[min_index], x_vals[i], y_vals[i])
        new_name= file_list[i].split('.')[0].replace('/','_').lower()
        df[new_name] = y_vals[i]
    if return_x:
        df['x'] = x_vals[min_index]
    if return_df:
        return df, df.corr()
    else:
        return df.corr()

def plot(file_list, file_type_list):
    # Computes the correlation between multiple spectra stored in files.
    # file_type_lists is a list of the types of files. Use sim for simulated files and exp for experimental files
    x_vals = []
    y_vals = []
    max_vals = []
    for file_name, file_type in zip(file_list, file_type_list):
        x, y, maximum = get_xy(file_name, file_type)
        x_vals.append(x)
        y_vals.append(y)
        max_vals.append(maximum)
        

    # Determine which set has more points -> longer must be mapped to smaller
    min_index = 0
    min_len = np.inf
    for i, y in enumerate(y_vals):
        if len(y) < min_len:
            min_len = len(y)
            min_index = i
    
    col_names = []
    df = pd.DataFrame()
    # Normalize all larger data to smaller
    for i in range(len(x_vals)):
        if i != min_index:
            y_vals[i] = normalize(max_vals[min_index], y_vals[i], max_vals[i])
            # Map big x and y values onto small 
            x_vals[i], y_vals[i] = getMapping(x_vals[min_index], y_vals[min_index], x_vals[i], y_vals[i])
        df[file_list[i].split('.')[0]] = y_vals[i]
        col_names.append(file_list[i].split('.')[0])
    df['x'] = x_vals[i]

    max_val = df[col_names].max().max()

    df.plot(x='x', y=col_names)
    plt.xlim([0, 3500])
    plt.ylim([0, max_val*1.1])
    plt.show()
    