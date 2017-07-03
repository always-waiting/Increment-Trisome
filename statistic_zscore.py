#!/usr/bin/python
# -*- coding: UTF-8 -*-
import glob
import pandas as pd
"""
收集输入目录中所有samplename.zscore文件，获取PZ值(即最后一列)，计算均值，方差
作为判断其他样本特性的参考集，
"""

inputdir = "/home/bixichao/Projects/tmp/Increment/mapping/"
ext = "zscore"
retfile = "./zscore"
chromlist = range(1,23)

zscore_list = []

for filename in glob.glob("%s/*.%s"%(inputdir, ext)):
    print filename
    samplename = filename.split("/")[-1].replace(".%s"%ext,"")
    data = pd.read_csv(filename, sep="\t",index_col=0, usecols=(0,3))
    data.rename(columns = {"PZ" : samplename}, inplace=True)
    zscore_list.append(data.squeeze())
zscore = pd.DataFrame(zscore_list)
mean = zscore.mean(axis=0)
std = zscore.std()
zscore.loc['mean'] = mean
zscore.loc['std'] = std
zscore.to_csv(retfile,sep="\t")
