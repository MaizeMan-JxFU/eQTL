import matplotlib.pyplot as plt
import pandas as pd
from tqdm import trange, tqdm
from scipy import stats
from multiprocessing import Pool
import numpy as np
from time import time

def delete_int_inStr(string:str=''):
    for i in range(10):
        string = string.replace(str(i), '')
    return string

class multi_treat0:
    def __init__(self, df:pd.DataFrame, tag:str):
        self.df, self.tag = df, tag
        pass
    def run(self, i):
        df, tag = self.df, self.tag
        df_subset = df[df.gene_tags==i]
        data0 = df_subset[df_subset['group_tags']==tag[0]]['values']
        data1 = df_subset[df_subset['group_tags']==tag[1]]['values']
        t,p_l = stats.levene(data0,data1)
        if p_l > 0.05:
            t,p = stats.ttest_ind(data0, data1)
        else:
            t,p = stats.ttest_ind(data0, data1, equal_var=False)
        p = -np.log10(p)
        fc = np.log2(np.mean(data1)/(np.mean(data0)))
        return [p, fc, i]

class Volcano:
    def __init__(self, df1:pd.DataFrame, group_li:list=[]):
        '''
        df1, 表达矩阵, index是基因编号
        
        '''
        if group_li ==[]:
            values = df1.values.flatten()
            tag = [delete_int_inStr(i)[0:2] for i in df1.columns]
            if len(set(tag)) == 2:
                print('自动检测成功！检测到两组标签:\n',set(tag))
            else:
                print('自动检测失败！请手动输入两组标签(group_li parameter)！')
            tag_cluster = []
            geneTag_cluster = []
            for i in trange(df1.shape[0]):
                tag_cluster += tag
                geneTag_cluster += [df1.index[i] for iteration in range(df1.shape[1])]
        df1_new = pd.DataFrame(np.array([values, geneTag_cluster, tag_cluster]).T, columns=['values', 'gene_tags', 'group_tags'])
        df1_new['values'] = df1_new['values'].astype(float)
        self.tag = list(set(tag))
        self.df = df1_new
        pass
    def Pcalculate(self, multi:int=None):
        '''
        multi, 多核计算运行核心数 (Optional)
        '''
        df = self.df
        tag = self.tag
        m_func = multi_treat0(df, tag)
        '''tag[1]比tag[2]'''
        if type(multi) == int:
            print(f'正在多核运算...\n线程数: {multi}')
            time_start = time()
            p = Pool(processes=multi)
            p_fc_tag_li = p.map(m_func.run, tqdm(set(df.gene_tags)))
            p.close()
            p.join()
            print(f'运算结束，耗时 {round(time()-time_start, 2)} sec')
        elif multi != int:
            print(f'正在单核运算...\n(如需多核运算，修改参数"multi=int")')
            time_start = time()
            p_fc_tag_li = []
            value_add = 0.000001
            for i in tqdm(set(df.gene_tags)):
                df_subset = df[df.gene_tags==i]
                data0 = df_subset[df_subset['group_tags']==tag[0]]['values']
                data1 = df_subset[df_subset['group_tags']==tag[1]]['values']
                t,p = stats.levene(data0,data1) # 方差齐检验
                if p > 0.05:
                    t,p = stats.ttest_ind(data0+value_add, data1+value_add) # 双尾t-test
                else:
                    t,p = stats.ttest_ind(data0+value_add, data1+value_add, equal_var=False) # 方差不齐，双尾t-test
                value_add = 0.000001
                fc = np.log2(round((value_add+np.mean(data1)/(value_add+np.mean(data0))), 4))
                p_fc_tag_li.append([p, fc, i])
            print(f'运算结束，耗时 {round(time()-time_start, 2)} sec')
        p_fc_tag_df = pd.DataFrame(np.array(p_fc_tag_li), columns=['-log10p', 'log2FoldChange', 'Tag']).dropna()
        p_fc_tag_df['-log10p'] = p_fc_tag_df['-log10p'].astype(float)
        p_fc_tag_df['log2FoldChange'] = p_fc_tag_df['log2FoldChange'].astype(float)
        p_fc_tag_df.to_csv('result.csv', index=None, sep='\t')
        return p_fc_tag_df
    def plot(self, p_fc_tag_df:pd.DataFrame, FoldChange:str='log2FoldChange', p:str='-log10p', tag:str='Tag'):
        print(p_fc_tag_df)
        plt.rc('font', family='Times New Roman')
        fig = plt.figure(figsize=[6,8], dpi=600)
        ax = fig.add_subplot()
        p_threashold = -np.log10(0.05)
        fc_threashold = 1
        # 上调基因
        dif_gene_up = p_fc_tag_df[(np.abs(p_fc_tag_df[p])>=p_threashold) & (p_fc_tag_df[FoldChange]>=fc_threashold)]
        ax.scatter(x=dif_gene_up[FoldChange], y=dif_gene_up[p], color='red', alpha=0.8, s=1)
        # 下调基因
        dif_gene_down = p_fc_tag_df[(np.abs(p_fc_tag_df[p])>=p_threashold) & (p_fc_tag_df[FoldChange]<=-fc_threashold)]
        ax.scatter(x=dif_gene_down[FoldChange], y=dif_gene_down[p], color='blue', alpha=0.8, s=1)
        # 不显著基因
        other_gene = p_fc_tag_df.loc[[i for i in p_fc_tag_df.index if i not in dif_gene_up.index and i not in dif_gene_down.index]]
        ax.scatter(x=other_gene[FoldChange], y=other_gene[p], color='grey', alpha=0.8, s=0.5)
        ax.hlines(xmin=min(p_fc_tag_df[FoldChange]), xmax=max(p_fc_tag_df[FoldChange]), linestyles='--', color='grey', y=-np.log10(0.05), linewidth=0.5)
        ax.vlines(ymin=min((p_fc_tag_df[p])), ymax=max(p_fc_tag_df[p]), linestyles='--', color='grey', x=1, linewidth=0.5)
        ax.vlines(ymin=min((p_fc_tag_df[p])), ymax=max(p_fc_tag_df[p]), linestyles='--', color='grey', x=-1, linewidth=0.5)
        ax.set_xlabel('logFC', fontsize=16, fontweight='bold')
        ax.set_ylabel('-log$_\mathdefault{10}$(p-value)', fontsize=16, fontweight='bold')
        for i in dif_gene_up.index:
            ax.annotate(dif_gene_up[tag].loc[i], (dif_gene_up[FoldChange].loc[i], dif_gene_up[p].loc[i]), xytext=(5, 5), textcoords='offset points', rotation=45, fontsize=2)
        for i in dif_gene_down.index:
            ax.annotate(dif_gene_down[tag].loc[i], (dif_gene_down[FoldChange].loc[i], dif_gene_down[p].loc[i]), xytext=(-5, -5), textcoords='offset points', rotation=-45, fontsize=2)

if __name__ == '__main__':
    pass
    # Char1 = np.random.random([10000,10])
    # Char1_df = pd.DataFrame(Char1, columns=['A'+str(i) for i in range(5)]+['B'+str(i) for i in range(5)])
    # volcano = Volcano(Char1_df, )
    # df = volcano.Pcalculate(4)
    # volcano.plot(df)