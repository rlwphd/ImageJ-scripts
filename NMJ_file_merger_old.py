# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import os

def group_func(s):
    """Add Group Info based on Mouse #"""
    if int(s) in [7488, 8834, 8835, 7482, 7485, 8831]:
        group = "C"
    elif int(s) in [8593, 8836, 8837, 8591, 8832, 8833]:
        group = "S"
    elif int(s) in [8594, 8595, 7591, 7592, 8592]:
        group = "T"
    elif int(s) in [123, 144, 151, 125, 135, 152]:
        group = "TC"
    else:
        group = "Other"
    
    return group

path = r'G:\Tetanus\NMJ Data\raw images\Analysis_Results'
dir_list = os.listdir(path)

filelist = [file for file in dir_list if file.endswith(".csv")]

df = pd.DataFrame()
for file in filelist:
    titlesplit = file.replace('-','_').split('_')
    
    df1 = pd.read_csv(os.path.join(path,file),index_col=0)
    df1.dropna(axis=0, how='all', inplace=True)
    
    titles = {'Mouse':[titlesplit[1] for i in range(len(df1))],
              'Muscle':[titlesplit[2][:-1].upper() for i in range(len(df1))],
              'Slide':[int(titlesplit[2][-1]) for i in range(len(df1))],
              'Slice':[int(titlesplit[3][-1]) for i in range(len(df1))],
              '405':[titlesplit[5].upper() if titlesplit[5][-3:] != '405' else titlesplit[5][:-3].upper() for i in range(len(df1))],
              '488':[titlesplit[6].upper() if titlesplit[6][-3:] != '488' else titlesplit[6][:-3].upper() for i in range(len(df1))],
              '555':[titlesplit[7].upper() if titlesplit[7][-3:] != '555' else titlesplit[7][:-3].upper() for i in range(len(df1))],
              '657':[titlesplit[8] if titlesplit[8][-3:] != '657' else titlesplit[8][:-3] for i in range(len(df1))],
              'Zoom':[titlesplit[9] for i in range(len(df1))]}
    
    df1=df1.assign(**titles)
    
    df = pd.concat([df,df1], ignore_index=True)

df["Innervated"] = df["results_array"]>10

totals = df.groupby(["Mouse","Muscle"])["Innervated"].count().rename("Counts")
sums = df.groupby(["Mouse","Muscle"])["Innervated"].sum()
result = pd.concat([totals,sums], axis=1)
result["Percentage"] = result["Innervated"]/result["Counts"]*100
result["Needs Data"] = result["Counts"] < 75
result.reset_index(inplace=True)
result["Group"] = result["Mouse"].apply(group_func)
result.to_csv(os.path.join(path, "NMJ_Result.csv"), index=False)

df.to_csv(os.path.join(path,"NMJ_Data.csv"), index=False)