#!/usr/bin/python
# -*- coding: UTF-8 -*-
import glob
import pandas as pd
import numpy as np
import argparse

"""
收集输入目录中所有samplename.zscore文件，获取PZ值(即最后一列)，计算均值，方差
作为判断其他样本特性的参考集，
"""

parse = argparse.ArgumentParser(description="收集输入目录中所有samplename.zscore文件，获取PZ值(即最后1列)，计算均值，方差作为判断其他样本特性的参考集")
parse.add_argument("--indir",'-idir',help="输入目录，存有.zscore文件",required=True)
parse.add_argument("--bin",'-b',help="窗口中bin数", default=1,type=int)
parse.add_argument("--male-list",'-ml',help="男性样本列表，为了处理性染色体",default="01-sample.male.list.txt", dest="malelist")
parse.add_argument("--female-list",'-fl', help="女性样本列表",default="01-sample.female.list.txt",dest="femalelist")
parse.add_argument("--ext","-e",help='查找文件的后缀,有一个%%d参数',default="b%d.zscore.txt")
args = parse.parse_args()

#inputdir = "/home/bixichao/Projects/tmp/Increment/mapping/"
#ext = "zscore"
#retfile = "./zscore"
inputdir = args.indir
win_bin = args.bin
ext = args.ext%args.bin
retfile = "zscore.b%d.txt"%args.bin
malelist = args.malelist
femalelist = args.femalelist
chromlist = range(1,25)

zscore_list = []

for filename in glob.glob("%s/*.%s"%(inputdir, ext)):
    print filename
    samplename = filename.split("/")[-1].replace(".%s"%ext,"")
    data = pd.read_csv(filename, sep="\t",index_col=0, usecols=(0,3))
    data.rename(columns = {"PZ" : samplename}, inplace=True)
    zscore_list.append(data.squeeze())
zscore = pd.DataFrame(zscore_list)

mean = zscore.mean(axis=0).astype(np.str)
std = zscore.std().astype(np.str)

# 处理X染色体
#malelist = "01-c-sample.male.list"
#femalelist = "01-c-sample.female.list"
maleindex = pd.read_csv(malelist,index_col=0, header=None).index
femaleindex = pd.read_csv(femalelist,index_col=0, header=None).index

mean_female_23 = zscore.ix[femaleindex].mean().astype(np.str)[23]
std_female_23 = zscore.ix[femaleindex].std().astype(np.str)[23]
mean_male_23 = zscore.ix[maleindex].mean().astype(np.str)[23]
std_male_23 = zscore.ix[maleindex].std().astype(np.str)[23]

mean[23] = "/".join([mean_female_23,mean_male_23])
std[23] = "/".join([std_female_23,std_male_23])
zscore.loc['mean'] = mean
zscore.loc['std'] = std
zscore.to_csv(retfile,sep="\t")
