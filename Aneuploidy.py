#!/usr/bin/python
# -*- coding:UTF-8 -*-
"""
通过我们的比对文件，以及自己生成的参考集，计算Aneuploidy(非整倍体)情况
输入文件有:
    *.gz.20K.GC.txt
    *.gz.20K.txt
    自制参考集:
        等价与Increment中winRatio_50k参考集
        正常样本计算的染色体PZ值均值与方差
"""

import pandas as pd
import glob
import math
import numpy as np
import argparse

def read_file(
        sample,\
        dirname='./',\
        rd_file_ext='gz.20K.txt',\
        gc_file_ext='gz.20K.GC.txt',\
        win_bin = 1
        ):
    """
    读取比对原始文件，返回需要的变量
    sample:         选取的样本
    dirname:        文件目录
    rd_file_ext:    rd文件后缀
    gc_file_ext:    gc文件后缀
    """
    print "Handling %s....." % sample
    total_loop = 12500/win_bin
    rd_file_regx = "%s/%s*.%s" % (dirname, sample, rd_file_ext)
    try:
        rd_file = glob.glob(rd_file_regx)[0]
    except IndexError,e:
        print rd_file_regx,"没有找到文件", e
        return
    gc_file_regx = "%s/%s*.%s" % (dirname, sample, gc_file_ext)
    try:
        gc_file = glob.glob(gc_file_regx)[0]
    except IndexError,e:
        print gc_file_regx, "没有找到文件", e
        return
    rd = pd.read_csv(rd_file, comment="#", index_col=0, sep="\t", header=None)
    rd_total = sum(rd.sum())
    rd_chr24 = sum(rd.ix['chr24'])
    y_per = (float(rd_chr24)/25652954)/(float(rd_total)/2859017332)
    gender = "M" if y_per > 0.16 else "F" # 0.2过于大了，许多都是0.19xxxx,也有0.188xxx,0.17xxx的
    gc = pd.read_csv(gc_file, comment="#", index_col=0, sep="\t", header=None)
    rd_merge = pd.DataFrame([rd[list(range(1+i*win_bin,1+(i+1)*win_bin))].sum(axis=1) for i in range(total_loop)],index = [range(1,total_loop + 1)]).T
    gc_merge = pd.DataFrame([gc[list(range(1+i*win_bin,1+(i+1)*win_bin))].sum(axis=1) for i in range(total_loop)],index = [range(1,total_loop + 1)]).T
    gc_per = gc_merge/(rd_merge*36)
    gc_per[gc_per.isnull()] = -1
    return rd_merge, gc_per, gender

def read_win_ref(
        ref_file,\
        dirname="./",\
        chromlist=range(1,23)
    ):
    """
    读取自建的increment参考集文件，并生成一个字典，键结构如下:
        chrom:
            start:
                end
                w_gc
                w_coe
                w_sd
    ref_file:   参考集文件名
    dirname:    参考集存放路径
    chromlist:  感兴趣染色体列表
    """
    try:
        refdata = pd.read_csv("%s/%s"%(dirname, ref_file), sep="\t", dtype={"ratio": np.str})
    except:
        print "解析%s/%s时出错"%(dirname,ref_file)
    refdict = {}
    for i,data in refdata.iterrows(): # 按行读取矩阵
        chrom = int(float(data['chr']))
        if (not chrom in chromlist): continue
        start = int(float(data['start']))
        innerdict = {'end': int(float(data['end'])), 'w_gc': data['gc'], 'w_coe': data['ratio'], 'w_sd': data['sd']}
        if refdict.has_key(chrom):
            refdict[chrom][start] = innerdict
        else:
            refdict[chrom] = {start: innerdict}
    return refdict
def read_zscore_ref(ref_file, dirname="./", chromlist=range(1,23)):
    """
    读取zscore参考集文件，获得每个染色体pz值的均值方差
    """
    filename = "%s/%s"%(dirname,ref_file)
    try:
        refdata = pd.read_csv(filename, sep="\t",index_col=0, skiprows=lambda x: (x!=201 and x!=0 and x!=202), dtype={"mean":np.str})
    except:
        print "解析%s/%s时出错"%(dirname,ref_file)
    refdict = {}
    for i,data in refdata.iterrows():
        for key, value  in data.to_dict().iteritems():
            key = int(key)
            if key == 23:
                (fz,mz) = value.split("/")
                if refdict.has_key(key):
                    refdict[key][i] = {"M":float(mz),"F":float(fz)}
                else:
                    refdict[key] = {i:{"M":float(mz),"F":float(fz)}}
            else:
                if refdict.has_key(key):
                    refdict[key][i] = float(value)
                else:
                    refdict[key] = {i:float(value)}
    return refdict

class ParseChromlist(argparse.Action):
    """
    def __init__(self,strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ParseChromlist, self).__init__(strings, dest, **kwargs)
    """
    def __call__(self, parser, namespace, values, option_string=None):
        #print '%r %r %r' % (namespace, values, option_string)
        valueslist = values.split(",")
        chromlist = []
        for item in valueslist:
            if item.find("..") != -1:
                (start, end) = item.split("..");
                chromlist.extend(range(int(start),int(end)+1))
            else:
                chromlist.append(eval(item))
        setattr(namespace, self.dest, set(chromlist))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="用于染色体非整倍体分析，或者生成zscore参考集\n建议大于5的认为是一个非整倍体")
    parser.add_argument("--list", "-l", help="样本列表", required=True)
    parser.add_argument("--inputdir", "-idir", help="比对序列存放路径",required=True)
    parser.add_argument("--refdir","-rdir",help="参考文件存放路径", default="./")
    parser.add_argument("--winref","-wref", help="窗口参考文件", default="ref-all")
    parser.add_argument("--zscoreref","-zref", help="zscore参考文件", default="zscore-all")
    parser.add_argument("--chromlist","-chr",help="感兴趣的染色体",default=set(range(1,25)), action=ParseChromlist)
    parser.add_argument("--bin", '-b', help="窗口bin数",default=1,type=int)
    parser.add_argument("--no-zscore-again",'-nza', help="不应用zscore参考集二次矫正", action="store_false", dest='zscoreflag')
    args = parser.parse_args()

    samplefile = args.list
    samplelist = pd.read_csv(samplefile, header=None)
    input_dir = args.inputdir
    ref_file = args.winref
    ref_dir = args.refdir
    zscore_ref_file = args.zscoreref
    chromlist = range(1,25)
    report_chromlist = args.chromlist
    win_bin = args.bin
    refdict = read_win_ref(ref_file, chromlist = chromlist, dirname=ref_dir)
    if args.zscoreflag:
        zscoreref = read_zscore_ref(zscore_ref_file, ref_dir, chromlist = chromlist)
    #print zscoreref
    for sample in samplelist.values:
        (rd, gc_per, gender) = read_file(sample[0], input_dir, win_bin = win_bin)
        get = (gc_per > 0) & (rd > 0)
        rd_get = rd[get]
        gc_per_get = gc_per[get]
        gc2rd = {}
        for index in rd_get.columns:
            pos = 20000*win_bin*(index-1) + 1
            for chrom in chromlist:
                if refdict[chrom].has_key(pos):
                    if math.isnan(rd_get[index]["chr%d"%chrom]):
                        refdict[chrom][pos]['rc'] = 1
                        refdict[chrom][pos]['gc'] = float("%.3f"%refdict[chrom][pos]['w_gc'])
                    else:
                        key = "%.3f"%gc_per_get[index]["chr%d"%chrom]
                        refdict[chrom][pos]['rc'] = int(rd_get[index]["chr%d"%chrom])
                        refdict[chrom][pos]['gc'] = float(key)
                        if chrom > 22: continue
                        if gc2rd.has_key(key):
                            gc2rd[key].append(int(rd_get[index]["chr%d"%chrom]))
                        else:
                            gc2rd[key] = [int(rd_get[index]["chr%d"%chrom])]
        # remove window when rc less than 11
        for key in gc2rd.keys():
            if len(gc2rd[key]) < 10: gc2rd.pop(key)
        # calculate correction coefficent
        all_median = np.median(sum(gc2rd.values(),[]))
        gc_rectify = {}
        for key,value in gc2rd.iteritems():
            gc_rectify[key] = all_median/np.median(value)
            if gc_rectify[key] > 2 or gc_rectify[key] < 0.4:
                gc_rectify[key] = 1
        # rectify the windows by gc
        total_read = {}
        for chrom in refdict.keys():
            total_read[chrom] = 0
            for pos in refdict[chrom].keys():
                gckey = "%.3f"%refdict[chrom][pos]['gc']
                if gc_rectify.has_key(gckey):
                    refdict[chrom][pos]['rc'] *= gc_rectify[gckey]
                total_read[chrom] += refdict[chrom][pos]['rc']
        # get all rd
        rd_stats = {}
        for chrom in refdict.keys():
            rd_stats[chrom] = {'ratio':[],'sd':[]}
            for pos in refdict[chrom].keys():
                rc = refdict[chrom][pos]['rc']
                if chrom == 23:
                    coe = float(refdict[chrom][pos]['w_coe'].split("/")[0] if gender == "F" else refdict[chrom][pos]['w_coe'].split("/")[1])
                else:
                    coe = float(refdict[chrom][pos]['w_coe'])
                sd = refdict[chrom][pos]['w_sd']
                if (coe == 0): continue
                percent = float("%.3f"%(rc/(sum(total_read.values()[0:22])*coe)))
                rd_stats[chrom]['ratio'].append(percent)
                rd_stats[chrom]['sd'].append(sd)

        all_ratio = sum(map(lambda x: rd_stats[x]['ratio'],rd_stats.keys()),[])
        t_mean = np.mean(all_ratio)
        t_sd = np.std(all_ratio)
        # zscore
        #msb = 0.0
        #mse = 0.0
        #total_number = 0
        with open("%s/%s.b%d.zscore.txt"%(input_dir,sample[0], win_bin),'w') as f:
            f.write("chr\tRatioMean\tSZ\tPZ\n")
            for chrom in report_chromlist:
                ratios = rd_stats[chrom]['ratio']
                sds = rd_stats[chrom]['sd']
                mean = np.median(ratios)
                #std_new = np.std(ratios)
                sd = np.std(sds)
                sd_sds = (sum(np.multiply(sds,sds))/len(sds))**0.5
                #print sd_sds
                sz = (mean - t_mean)/t_sd
                #sz1 = (mean - t_mean)/std_new
                pz = (mean - 1)/sd
                pz_new = (mean -1)/sd_sds

                #msb = msb + (mean - t_mean)**2*len(ratios)
                #mse = mse + np.var(ratios)*len(ratios)
                #total_number = len(ratios) + total_number

                f.write("%d\t%.3f\t%.3f\t%.3f\n"%(chrom,mean,sz,pz))
                if args.zscoreflag:
                    if chrom == 23:
                        print "%d\t%.3f\t%.3f"%(chrom, mean, (pz-zscoreref[chrom]['mean'][gender])/zscoreref[chrom]['std'][gender])
                    else:
                        print "%d\t%.3f\t%.3f"%(chrom, mean, (pz-zscoreref[chrom]['mean'])/zscoreref[chrom]['std'])
                else:
                    print "%d\t%.3f\t%.3f"%(chrom, mean, pz_new)
                    #print "| %d | %.3f | %.3f | %.3f | %.3f | %.3f |"%(chrom,mean, sd, sd_sds, pz, pz_new)
                    #print "| %d | %.3f | %.3f | %.3f | %.3f | %.3f"%(chrom,mean, t_mean, t_sd, sz, std_new)
                    #print "| %d | %.3f | %.3f | %.3f | %.3f | %.3f"%(chrom,mean, t_mean, t_sd, sz1, std_new)


        #print "MSB: %.3f" % (msb/23)
        #print "MSE: %.3f" % (mse/total_number)
        #print "Total: %d" % total_number




