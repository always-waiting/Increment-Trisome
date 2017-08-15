# -*- coding:UTF-8 -*-
"""
Berry-Increment-Trisome:
    用于判断三整倍体的自动分析程序，分析过程如下:
    1.解析CNV分析结果目录，每个样本会生成一个'样本名.zscore'文件
    2.判断分析完成后，会在对应批次目录中生成标志文件，确认分析过这个批次
    3.把'样本名.zscore'文件导入xromate系统
"""
import os
import glob
import datetime
import pandas as pd
import numpy as np
import math
import pymongo
import logging

logger = logging.getLogger("default")
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
formatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
console.setFormatter(formatter)
logger.addHandler(console)

class BerryIncrementTrisomeException(Exception):
    def __init__(self, message):
        super(BerryIncrementTrisomeException, self).__init__(message)

class dirobj(object):
    """
    目录对象
    """
    def __init__(self, path):
        if os.path.isdir(path):
            self.path = path
        else:
            raise BerryIncrementTrisomeException("%s is not a path"%path)

    def filter(self, func):
        """
        按照func过滤目录，返回目录名,以列表形式，没有返回空列表
        """
        fid = []
        dirs = os.listdir(self.path)
        for flowcell in dirs:
            fiddir = os.path.join(self.path,flowcell)
            if func(fiddir):
                fid.append(flowcell)
        return fid

class BerryIncrementTrisomeAuto(object):
    """
    自动分析类,有如下属性:
        indir: 监控的目录，确保目录中的批次都进行了Trisome分析，必要参数
        date: 天数，形式如3,意思是处理距离当前时间3天内数据,默认为3
        endfile: 分析完成标志文件，默认为trisome.end
        samplefile: 每个样本的分析结果文件，默认为'%s.b%d.zscore.txt'
    """
    defaultvalue = {
        "date"              : 3,
        "endfile"           : "trisome.end",
        "samplefile"        : "%s.b%d.zscore.txt",
        "bin"               : 1,
        "rd_ext"            : "gz.20K.txt",
        "gc_ext"            : "gz.20K.GC.txt",
        "refdir"            : os.path.dirname(os.path.abspath(__file__))+"/ref",
        "winref"            : "ref-b%d.txt",
        "zscoreref"         : "zscore.b%d.txt",
        "chromlist"         : set(range(1,25)),
        "report_chromlist"  : set(range(1,25)),
        "za"                : True,
        "mongourl"          : "mongodb://localhost:27017/",
        "db"                : "xromate_mandel",
        "coll"              : "samples",
        "insert_db"         : False,
        "logger"            : logger
    }
    def __init__(self,indir,*args,**kwargs):
        self.indir = dirobj(indir)
        self.flowcell = []
        for key in BerryIncrementTrisomeAuto.defaultvalue.keys():
            if key in kwargs.keys():
                self.__dict__[key] = kwargs[key]
            else:
                self.__dict__[key] = BerryIncrementTrisomeAuto.defaultvalue[key]
        if self.winref.find("%") != -1:
            self.winref = self.winref%self.bin
        if self.zscoreref.find("%") != -1:
            self.zscoreref = self.zscoreref%self.bin
    def detect_flowcell(self):
        """
        扫描indir目录，根据规则，判断要分析的flowcell,返回需要分析的批次
        """
        def __filter_rule(fiddir):
            """
            输入参数为需要判断的目录
            批次选择规则
            1. 时间在self.date内
            2. 有CNV_B1目录 - 保证cnv已经分析完成后再做这个分析
            2. 批次目录中没有self.endfile文件
            3. 批次目录中，不是每个样本都有self.samplefile文件
            """
            today = datetime.date.today()
            ctime = datetime.date.fromtimestamp(os.stat(fiddir).st_mtime)
            if today - ctime > datetime.timedelta(self.date):
                return False
            if not os.path.isdir(os.path.join(fiddir,"CNV_B1")):
                return False
            if os.path.isfile(os.path.join(fiddir, self.endfile)):
                return False
            for filename in glob.glob(os.path.join(fiddir,"*.R1.clean.fastq.gz.txt")):
                samplename = filename.split("/")[-1].split("_")[0]
                samplefile = os.path.join(fiddir, self.samplefile%(samplename, self.bin))
                if(os.path.isfile(samplefile)):
                    pass
                else:
                    return True
            return True
        if len(self.flowcell) != 0:
            self.flowcell = []
        self.flowcell.extend(self.indir.filter(__filter_rule))
        return self.flowcell
    def analyse_flowcell(self):
        for flowcell in self.flowcell:
            self.logger.info("Analyzing %s ...." % flowcell)
            kwargs = self.__dict__.copy()
            kwargs.pop('flowcell')
            obj = BerryIncrementTrisome(flowcell, **kwargs)
            obj.pipeline()

class BerryIncrementTrisome(object):
    """
    单独flowcell的分析程序
    """
    defaultvalue = {
        "indir"             : "/seq_dir/staff/deploy/CNV",
        "rd_ext"            : "gz.20K.txt",
        "gc_ext"            : "gz.20K.GC.txt",
        "samplefile"        : "%s.b%d.zscore.txt",
        "refdir"            : os.path.dirname(os.path.abspath(__file__))+"/ref",
        "winref"            : "ref-finally.txt",
        "zscoreref"         : "zscore.b1.txt",
        "chromlist"         : set(range(1,25)),
        "report_chromlist"  : set(range(1,25)),
        "endfile"           : "trisome.end",
        "bin"               : 1,
        "za"                : True,
        "mongourl"          : "mongodb://localhost:27017/",
        "db"                : "xromate_mandel",
        "coll"              : "samples",
        "insert_db"         : True,
        "logger"            : logger
    }
    def __init__(self,flowcell, *args, **kwargs):
        self.flowcell = flowcell
        for key in BerryIncrementTrisome.defaultvalue.keys():
            if key in kwargs.keys():
                if key == "indir":
                    if isinstance(kwargs[key],dirobj):
                        self.__dict__[key] = kwargs[key]
                    elif isinstance(kwargs[key],str):
                        self.__dict__[key] = dirobj(kwargs[key])
                    else:
                        raise BerryIncrementTrisomeException("indir should be str or dirobj")
                else:
                    self.__dict__[key] = kwargs[key]
            else:
                if key == "indir":
                    self.__dict__[key] = dirobj(BerryIncrementTrisome.defaultvalue[key])
                else:
                    self.__dict__[key] = BerryIncrementTrisome.defaultvalue[key]
        self.__dict__["indir"] = dirobj(self.indir.path + "/" + self.flowcell)
        self.samplelist = []
        self.sampleresult = {}
        self.refdict = {}
        self.zref = {}
    def detect_sample(self):
        self.__detect_sample()

    def __pre_analyze(self):
        self.logger.info("pre analyzing...")
        self.__detect_sample()
        self.__read_win_ref()
        if self.za:
            self.__read_zscore_ref()

    def __read_zscore_ref(self):
        """
        读取zscore参考集文件，获得每个染色体pz值的均值方差
        """
        filename = "%s/%s"%(self.refdir, self.zscoreref)
        try:
            refdata = pd.read_csv(filename, sep="\t",index_col=0, dtype={"mean":np.str})
        except:
            raise BerryIncrementTrisomeException("解析%s时出错"%(filename))
        for i,data in refdata.iterrows():
            for key, value  in data.to_dict().iteritems():
                key = int(key)
                if key == 23:
                    (fz,mz) = value.split("/")
                    if self.zref.has_key(key):
                        self.zref[key][i] = {"M":float(mz),"F":float(fz)}
                    else:
                        self.zref[key] = {i:{"M":float(mz),"F":float(fz)}}
                else:
                    if self.zref.has_key(key):
                        self.zref[key][i] = float(value)
                    else:
                        self.zref[key] = {i:float(value)}
        return self.zref

    def __read_win_ref(self):
        """
        读取自建的increment参考集文件，并生成一个字典，键结构如下:
            chrom:
                start:
                    end,w_gc,w_coe,w_sd
        """
        try:
            refdata = pd.read_csv("%s/%s"%(self.refdir, self.winref), sep="\t", dtype={"ratio": np.str})
        except:
            raise  BerryIncrementTrisomeException("解析%s/%s时出错"%(self.refdir,self.winref))
        for i,data in refdata.iterrows(): # 按行读取矩阵
            chrom = int(float(data['chr']))
            if (not chrom in self.chromlist): continue
            start = int(float(data['start']))
            innerdict = {'end': int(float(data['end'])), 'w_gc': data['gc'], 'w_coe': data['ratio'], 'w_sd': data['sd']}
            if self.refdict.has_key(chrom):
                self.refdict[chrom][start] = innerdict
            else:
                self.refdict[chrom] = {start: innerdict}
        return self.refdict

    def __detect_sample(self):
        if len(self.samplelist) != 0:
            self.samplelist = []
        if len(self.sampleresult) != 0:
            self.sampleresult = {}
        for samplefile in glob.glob(os.path.join(self.indir.path,"*.R1.clean.fastq.gz.txt")):
            samplename = samplefile.split("/")[-1].split("_")[0]
            self.samplelist.append(samplename)
    def __read_file(self, sample):
        """
        读取比对原始文件，返回需要的变量
        """
        self.logger.info("Reading %s"%sample)
        total_loop = 12500/self.bin
        rd_file_regx = "%s/%s*.%s"%(self.indir.path,sample,self.rd_ext)
        try:
            rd_file = glob.glob(rd_file_regx)[0]
        except IndexError,e:
            raise BerryIncrementTrisomeException("%s没有找到"%rd_file_regx)
        gc_file_regx = "%s/%s*.%s"%(self.indir.path,sample,self.gc_ext)
        try:
            gc_file = glob.glob(gc_file_regx)[0]
        except IndexError,e:
            raise BerryIncrementTrisomeException("%s没有找到"%gc_file_regx)
        rd = pd.read_csv(rd_file, comment="#", index_col=0, sep="\t", header=None)
        rd_total = sum(rd.sum())
        rd_chr24 = sum(rd.ix['chr24'])
        y_per = (float(rd_chr24)/25652954)/(float(rd_total)/2859017332)
        gender = "M" if y_per > 0.16 else "F" # 0.2过于大了，许多都是0.19xxxx,也有0.188xxx,0.17xxx的
        gc = pd.read_csv(gc_file, comment="#", index_col=0, sep="\t", header=None)
        rd_merge = pd.DataFrame([rd[list(range(1+i*self.bin,1+(i+1)*self.bin))].sum(axis=1) for i in range(total_loop)],index = [range(1,total_loop + 1)]).T
        gc_merge = pd.DataFrame([gc[list(range(1+i*self.bin,1+(i+1)*self.bin))].sum(axis=1) for i in range(total_loop)],index = [range(1,total_loop + 1)]).T
        gc_per = gc_merge/(rd_merge*36)
        gc_per[gc_per.isnull()] = -1
        return rd_merge, gc_per, gender
    def pipeline(self):
        self.__pre_analyze()
        self.__analyze()
        self.__post_analyze()

    def __analyze(self):
        for sample in self.samplelist:
            self.logger.info("Analysing %s"%sample)
            (rd, gc_per, gender) = self.__read_file(sample)
            #if not rd: continue # 某个样本读取文件出错!
            get = (gc_per > 0) & (rd > 0)
            rd_get = rd[get]
            gc_per_get = gc_per[get]
            gc2rd = {}
            for index in rd_get.columns:
                pos = 20000*self.bin*(index-1) + 1
                for chrom in self.chromlist:
                    if self.refdict[chrom].has_key(pos):
                        if math.isnan(rd_get[index]["chr%d"%chrom]):
                            self.refdict[chrom][pos]['rc'] = 1
                            self.refdict[chrom][pos]['gc'] = float("%.3f"%self.refdict[chrom][pos]['w_gc'])
                        else:
                            key = "%.3f"%gc_per_get[index]["chr%d"%chrom]
                            self.refdict[chrom][pos]['rc'] = int(rd_get[index]["chr%d"%chrom])
                            self.refdict[chrom][pos]['gc'] = float(key)
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
            for chrom in self.refdict.keys():
                total_read[chrom] = 0
                for pos in self.refdict[chrom].keys():
                    gckey = "%.3f"%self.refdict[chrom][pos]['gc']
                    if gc_rectify.has_key(gckey):
                        self.refdict[chrom][pos]['rc'] *= gc_rectify[gckey]
                    total_read[chrom] += self.refdict[chrom][pos]['rc']
            # get all rd
            rd_stats = {}
            for chrom in self.refdict.keys():
                rd_stats[chrom] = {'ratio':[],'sd':[]}
                for pos in self.refdict[chrom].keys():
                    rc = self.refdict[chrom][pos]['rc']
                    if chrom == 23:
                        coe = float(self.refdict[chrom][pos]['w_coe'].split("/")[0] if gender == "F" else self.refdict[chrom][pos]['w_coe'].split("/")[1])
                    else:
                        coe = float(self.refdict[chrom][pos]['w_coe'])
                    sd = self.refdict[chrom][pos]['w_sd']
                    if (coe == 0): continue
                    percent = float("%.3f"%(rc/(sum(total_read.values()[0:22])*coe)))
                    rd_stats[chrom]['ratio'].append(percent)
                    rd_stats[chrom]['sd'].append(sd)

            all_ratio = sum(map(lambda x: rd_stats[x]['ratio'],rd_stats.keys()),[])
            t_mean = np.mean(all_ratio)
            t_sd = np.std(all_ratio)
            self.logger.info("Writing %s result..."%sample)
            with open(os.path.join(self.indir.path,self.samplefile%(sample,self.bin)),'w') as f:
                f.write("chr\tRatioMean\tSZ\tPZ\n")
                for chrom in self.report_chromlist:
                    ratios = rd_stats[chrom]['ratio']
                    sds = rd_stats[chrom]['sd']
                    mean = np.median(ratios)
                    sd = np.std(sds)
                    sz = (mean - t_mean)/t_sd
                    pz = (mean - 1)/sd
                    #f.write("%d\t%.3f\t%.3f\t%.3f\n"%(chrom,mean,sz,pz))
                    if self.za:
                        if chrom == 23:
                            pz_write = (pz-self.zref[chrom]['mean'][gender])/self.zref[chrom]['std'][gender]
                        else:
                            pz_write = (pz-self.zref[chrom]['mean'])/self.zref[chrom]['std']
                    else:
                        pz_write = pz
                    f.write("%d\t%.3f\t%.3f\t%.3f\n"%(chrom,mean,sz,pz_write))
                    if self.sampleresult.has_key(sample):
                        self.sampleresult[sample][chrom] = pz_write
                    else:
                        self.sampleresult[sample] = {chrom : pz_write}

    def __post_analyze(self):
        self.__write_endfile()
        if self.insert_db:
            self.__import_mongo()
    def __write_endfile(self):
        with open(os.path.join(self.indir.path,self.endfile),"w") as f:
            map(lambda x: f.write("%s\n"%x), self.samplelist)
    def import_mongo(self):
        self.__detect_sample()
        for sample in self.samplelist:
            self.sampleresult[sample] = {}
            samplefile = os.path.join(self.indir.path,self.samplefile%(sample,self.bin))
            try:
                with open(samplefile) as f:
                    f.readline()
                    for line in f.readlines():
                        lines = line.strip().split("\t")
                        self.sampleresult[sample][int(lines[0])] = float(lines[3])
            except:
                raise BerryIncrementTrisomeException("Failed when read %s"%samplefile)
        self.__import_mongo()
    def __import_mongo(self):
        client = pymongo.MongoClient(self.mongourl)
        db = client[self.db]
        coll = db[self.coll]
        for sample in self.samplelist:
            self.logger.info("Importing %s"%sample)
            viewdict = {}
            try:
                for chrom in self.sampleresult[sample]:
                    if chrom == 23:
                        viewdict["trisomescore.X"] = "%.3f"%self.sampleresult[sample][chrom]
                    elif chrom == 24:
                        viewdict["trisomescore.Y"] = "%.3f"%self.sampleresult[sample][chrom]
                    else:
                        viewdict["trisomescore.%d"%chrom] = "%.3f"%self.sampleresult[sample][chrom]
                coll.update_one({"name":sample},{"$set": viewdict}, upsert=False)
            except:
                self.logger.info("%s import %s failed"%(sample, self.mongourl))
                self.logger.info(viewdict)
