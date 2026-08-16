import os,sys,re
import pandas as pd
import argparse

__author__ = 'zhangyu'
__mail__ = 'zhangyu@genomeprecision.com'
__date__ = '2021-08-05'

def get_group(compare_name):
    group_tre, group_ck = re.split('vs',compare_name)
    return group_tre, group_ck

def set_list(list_a, list_b):
    overlap_ele = list(set(list_a).intersection(set(list_b)))
    return overlap_ele

def judge_compare(pd_condition, group_a, group_b):
    location_a = (pd_condition == group_a).any()
    location_a = location_a.index[location_a][0]
    group_a_sample_pd = pd_condition[location_a].where(pd_condition[location_a] == group_a).dropna()
    group_a_sample = group_a_sample_pd.index._data

    location_b = (pd_condition == group_b).any()
    location_b = location_b.index[location_b][0]
    group_b_sample_pd = pd_condition[location_b].where(pd_condition[location_b] == group_b).dropna()
    group_b_sample = group_b_sample_pd.index._data
    #
    # print(group_a_sample, location_a)
    # print(group_b_sample, location_b)

    duplicate_elements = set_list(group_a_sample, group_b_sample)

    if duplicate_elements != []:
        print('ERROR:group_a:{0} and group_b:{2} have duplicate elements:{1},pleace check it!!!'.format(group_a,duplicate_elements,group_b))
        sys.exit(1)

    condition_pd = pd.DataFrame(pd.concat([group_a_sample_pd,group_b_sample_pd]))
    condition_pd.rename(columns={0:'group'},inplace = True)
    condition_pd.columns = ['group',]
    
    return condition_pd

def write_to_file(condition_pd, condition_file):
    outdir = os.path.dirname(os.path.realpath(condition_file))
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except FileExistsError as e:
            print('dir {0} is exists,not duplicated creating'.format(outdir))
        except:
            os.system('mkdir -p {0}'.format(outdir))
    condition_pd.to_csv(condition_file, sep='\t', mode='w')

def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__, __mail__,__date__))

    parser.add_argument('-c', '--compare', help='compare name',dest='compare', required=True)
    parser.add_argument('-o', '--outfile', help='output file', metavar='',dest='outfile',required=True)
    parser.add_argument('-i', '--infile', help='condition file', metavar='',dest='infile',required=True)

    args = parser.parse_args()

    complex_condition = pd.read_csv(args.infile, header=0, sep='\t', index_col=0)
    group_tre, group_ck = get_group(args.compare)
    condition_pd = judge_compare(complex_condition, group_tre, group_ck)
    write_to_file(condition_pd, args.outfile)

if __name__ == '__main__':
    if len(sys.argv) < 1:
        print(__doc__)
        sys.exit(1)

    main()
