import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
import seaborn as sns
from  multiprocessing import Pool

def sort_jxfu(input_li:list=[]):
    '''
    带数字的排前面，字符串拍后面
    分别进行排序
    '''
    output_li = sorted([i for i in input_li if type(i)==int])+sorted([i for i in input_li if type(i)==str])
    return output_li

def num_int_list(input_li:list=[]):
    '''
    将列表中的字符串形式数字转化为整数形式
    '''
    output_li = []
    for i in input_li:
        if i in [str(i) for i in range(30)]:
            output_li.append(int(i))
        else:
            output_li.append(i)
    return output_li

def Manhattan_plot(df:pd.DataFrame, chr:str, pos:str, pvalue:str, ChrTag_len:int=5, threshold:float=None, annote:str=None):
    '''
    必须输入的项目df: 矩阵, chr: 染色体, pos: SNP位点, pvalue: 显著性
    '''
    matplotlib.use('agg') # 绘图加速
    plt.rc('font', family='Times New Roman')
    fig = plt.figure(figsize=[12,6], dpi=600)
    ticks_loc = []
    df[chr] = num_int_list(df[chr])
    Chr = sort_jxfu([i for i in set(df[chr]) if len(str(i))<=ChrTag_len])
    colors = sns.color_palette("hls", len(Chr))
    for i in tqdm(Chr):  # 对每条染色体分别进行散点图绘制
        df_tmp = df[df[chr]==i].sort_values(by=[pos])
        site = df_tmp[pos]
        if Chr.index(i)>0:
            ticks_loc.append(len(site)/2+max(x)+0.01*df.shape[0])
            x = [i+max(x)+0.01*df.shape[0] for i in range(len(site))] # 每条染色体之间间隔SNP位点总数的1%
        else:
            x = range(len(site))
            ticks_loc.append(len(site)/2)
        p = df_tmp[pvalue]
        y = -1*np.log10(p)
        # print(f'正在绘制 CHR-{i}')
        # print(f'共有{len(y)}个点')
        color=colors[Chr.index(i)]
        plt.scatter(x, y, alpha=0.2, s=6, color=color)
        
        if threshold != None and max(y)>threshold:
            df_tmp['x'] = x
            df_tmp['y'] = y
            df_annote = df_tmp[df_tmp['y']>threshold]
            plt.scatter(df_annote['x'], df_annote['y'], alpha=1, s=12, color=color)
            if annote == None:
                annote == 'Marker'
            for i in df_annote.index:
                plt.annotate(df_annote[annote].loc[i], (df_tmp['x'].loc[i], df_tmp['y'].loc[i]), xytext=(5, 5), textcoords='offset points')
            
    if threshold != None:
        plt.hlines(y=threshold, xmin=0, xmax=max(x),color='grey', linewidth=1, alpha=0.8)
    plt.xticks(ticks_loc, Chr)
    plt.xlim([0-0.01*max(x),max(x)+0.01*max(x)])
    plt.ylabel('-log$_\mathdefault{10}$(p-value)', fontsize=16, fontweight='bold')
    plt.xlabel('Chromosome', fontsize=16, fontweight='bold')
    return

if __name__ == '__main__':
    path = '/share/org/YZWL/yzbsl_jiaaq/jxfu/Pyscript/github/eQTL/BE_TotalPositionSummary.add_beta.txt'
    data = pd.read_csv(path, delimiter='\t')
    Manhattan_plot(data, 'Chr', 'Site', 'p', 2, threshold=42)
    plt.savefig('曼哈顿图.png', transparent=True)