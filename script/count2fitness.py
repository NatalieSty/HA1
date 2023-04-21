import os
import glob
import pandas as pd
import numpy as np
import argparse

def norm_amp(count_DF,col_name):
    count_DF[col_name]=count_DF[col_name]+1
    count_DF[col_name+'_norm'] = ((count_DF[col_name]) / count_DF.groupby('Amplicon')[col_name].transform('sum')
)
    return count_DF
def cal_fit(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = np.log10((count_DF[column_name])/(count_DF[input]))
    silent_DF=count_DF[count_DF.Mutation.str.contains('silent')]
    silent_mean=silent_DF[column_name+'_enri'].mean()
    silent_median = silent_DF[column_name + '_enri'].median()
    count_DF[column_name + '_enri']=(count_DF[column_name + '_enri']-silent_mean)
    count_DF=count_DF.rename(columns={column_name + '_enri':column_name + '_fitness'})
    return count_DF
def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF

def grab_files_as_Dict(path,suffix):
    #path = 'result/*/*'
    #pattern = '_count.tsv'
    files = glob.glob(path+suffix)
    #filenames = [os.path.basename(i).split(suffix, 1)[0] for i in files]
    file_dict={}
    for file in files:
        filename=os.path.basename(file).split(suffix, 1)[0]
        file_dict[filename]=file
    return file_dict

def main():
    'usage: python count2fitness.py -f 1F04_escape_rep1_GTCGTGAT-AAGGTTCA 1F04_escape_rep2_ACCACTGT-ACTTAGCA -o ../result/Mos99_fit.csv'
 
    parser = argparse.ArgumentParser()
    parser.add_argument("-path", "--path_to_count", default='../result/*/*', type=str,help='path that contatins the count files')
    parser.add_argument("-suf", "--suffix", default='_count.tsv', type=str,help='suffix for the count files')
    parser.add_argument("-o", "--output", default='../result/Mos99_fit.csv', type=str,help='output file name')
    parser.add_argument("-f", "--files", nargs='+', default=["3A10_escape_rep1_AGGTTGGG-TCGATTAG", "3A10_escape_rep2_GTGTGGTG-TGTTCTAG","1F04_escape_rep1_GTCGTGAT-AAGGTTCA", "1F04_escape_rep2_ACCACTGT-ACTTAGCA","Mos99_HK14HA_rep1_TGGGTTTC-GGAACTTA", "Mos99_HK14HA_rep2_TGGTCACA-TAGGTCTA"],help='list the files that are needed to calculate fitness')
    parser.add_argument("-n", "--names", nargs='+', default=["WT_plasmid_GTGTGGTG-TGTTCTAG","Input_AGGTTGGG-TCGATTAG"],help='file name for WT and Input')
    args = parser.parse_args()


    file_dict=grab_files_as_Dict(args.path_to_count,args.suffix)
    DF=pd.DataFrame(columns=['Mutation'])

    for filename,file in file_dict.items():
        df= pd.read_csv(file, header=0, sep= '\t')
        df.columns=['Mutation',filename]
        DF=pd.merge(DF,df,how='outer',on='Mutation')
    #DF.to_csv('result/Mos99_all_count.csv')
    # split amplicon and mutation
    
    DF[['Amplicon','Mutation']] = DF.Mutation.str.split("|", n=1, expand=True)
    DF=DF.fillna(0)
    #normalization
    norm_ls = args.names
    for col in norm_ls:DF = norm_amp(DF, col)
    # only one mutation & filter low input
    onemut_df = DF[(~DF["Mutation"].str.contains('-')) & (DF[args.names[1]] >= 10)]
    # filter out WT high error mutants
    onemut_df = onemut_df[onemut_df[args.names[1]+'_norm'] > onemut_df[args.names[0]+'_norm'] * 6]
    #calculate fitness
    fit_ls=args.files

    for c in fit_ls:
        onemut_df=norm_amp(onemut_df,c)
        onemut_df=cal_fit(onemut_df,c+'_norm',args.names[1]+'_norm')
    #cal_mean(onemut_df,'Fitness','Rep1_norm_fitness','Rep2_norm_fitness')
    #DF.to_csv('result/Mos99_fit_raw.csv')
    # classify mutation type(silent;missense;nonsense)
    onemut_df['mutation_type'] = 'missense'
    onemut_df['mutation_type'][onemut_df.Mutation.str.contains('silent')] = 'silent'
    onemut_df['mutation_type'][(~onemut_df.Mutation.str.contains('silent')) & onemut_df.Mutation.str.contains('_')] = 'nonsense'

    pos_df = onemut_df.Mutation.str.extract('(\d+)')
    one_mut_sortdf = onemut_df.join(pos_df, lsuffix='_caller', rsuffix='_other') \
        .set_index(0) \
        .fillna(0) \
        .sort_index()
    one_mut_sortdf.to_csv(args.output)


if __name__ == "__main__":
    main()

