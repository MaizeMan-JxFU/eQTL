import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
import seaborn as sns
from  multiprocessing import Pool
from matplotlib.gridspec import GridSpec

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

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

def Manhattan_plot(df:pd.DataFrame, chr:str, pos:str, pvalue:str, ChrTag_len:int=5, threshold:float=None, annote:str=None, color_set:list=[]):
    '''
    必须输入的项目df: 矩阵, chr: 染色体, pos: SNP位点, pvalue: 显著性
    '''
    matplotlib.use('agg') # 绘图加速
    plt.rc('font', family='Times New Roman')
    fig = plt.figure(figsize=[12,6], dpi=600)
    gs = GridSpec(12, 1, figure=fig)
    
    ax = fig.add_subplot(gs[0:12,0])
    ticks_loc = []
    df[chr] = num_int_list(df[chr])
    Chr = sort_jxfu([i for i in set(df[chr]) if len(str(i))<=ChrTag_len])
    if color_set == []:
        colors = sns.color_palette("hls", len(Chr))
    else:
        colors = color_set
    for i in tqdm(Chr):  # 对每条染色体分别进行散点图绘制
        df_tmp = df[df[chr]==i].sort_values(by=[pos])
        site = df_tmp[pos]
        if Chr.index(i)>0:
            ticks_loc.append(len(site)/2+max(x)+0.008*df.shape[0])
            x = [i+max(x)+0.008*df.shape[0] for i in range(len(site))] # 每条染色体之间间隔SNP位点总数的1%
        else:
            x = range(len(site))
            ticks_loc.append(len(site)/2)
        p = df_tmp[pvalue]
        y = -1*np.log10(p)
        # print(f'正在绘制 CHR-{i}')
        # print(f'共有{len(y)}个点')
        color=colors[Chr.index(i)]
        ax.scatter(x, y, alpha=1, s=6, color=lighten_color(color, 0.8))
        
        if threshold != None and max(y)>threshold:
            df_tmp['x'] = x
            df_tmp['y'] = y
            df_annote = df_tmp[df_tmp['y']>threshold]
            ax.scatter(df_annote['x'], df_annote['y'], alpha=1, s=12, color=color)
            if annote == None:
                annote == 'Marker'
            for i in df_annote.index:
                ax.annotate(df_annote[annote].loc[i], (df_tmp['x'].loc[i], df_tmp['y'].loc[i]), xytext=(5, 5), textcoords='offset points', rotation=45, fontsize=8)
            
    if threshold != None:
        ax.hlines(y=threshold, xmin=0, xmax=max(x),color='grey', linewidth=1, alpha=0.8, linestyles='--')
    ax.set_xticks(ticks_loc, [f'Chr{i}' for i in Chr])
    ax.set_xlim([0-0.01*max(x),max(x)+0.01*max(x)])
    # ax.set_ylim([min(y)-0.01*max(y), max(y)+0.01*max(y)])
    ax.set_ylabel('-log$_\mathdefault{10}$(p-value)', fontsize=16, fontweight='bold')
    # ax.set_xlabel('Chromosome', fontsize=16, fontweight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    
    # ax2 = fig.add_subplot(gs[14,0], sharex=ax)
    # ax2.spines['right'].set_visible(False)
    # ax2.spines['left'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    return ax

if __name__ == '__main__':
    path = '/share/org/YZWL/yzbsl_jiaaq/jxfu/Pyscript/github/eQTL/BE_TotalPositionSummary.add_beta.txt'
    data = pd.read_csv(path, delimiter='\t')
    Manhattan_plot(data, 'Chr', 'Site', 'p', 2, threshold=42)
    plt.savefig('曼哈顿图.png', transparent=True)