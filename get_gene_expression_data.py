#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import pickle
import optparse
import glob
import os



class DataQC:
 
  def __init__(self,exp_type,ehk,strict):
    self.exp_type = exp_type
    self.ehk = ehk
    self.strict =strict
    if strict=='tcga':
      self.QC= '/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/RNAseq_tools/ALLPlates_MAY21_1_23_rna_seq_stats_extra.tsv'
    elif strict =='new':
      self.QC = '/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/RNAseq_tools/ALLPlates_Sept4_1_25_rna_seq_stats_new.tsv'
    else:  
      self.QC = '/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/RNAseq_tools/ALLPlates_Sept4_1_25_rna_seq_stats.tsv'
    self.counts = '/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/QC/ALLPlates_Sept4_1_25_counts.tsv'
    self.tpm = '/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/QC/ALLPlates_Sept4_1_25_tpm.tsv'
    self.logCPM ='/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/QC/ALLPlates_Sept4_1_25_logCPM.tsv'#Mar13_logCPM_plates13_22.tsv
    self.read_count = 10
    
  def get_data(self):
    '''
    Returns a df of expression type requested that pass QC at the
    given threshhold and filter out low count genes
    For counts data has the additional filtering of low count genes (gene
    with a mean count < read_count)
    '''
    print "in bin"
    if self.exp_type == 'TPM':
      df = pd.DataFrame.from_csv(self.tpm, sep=('\t'),header=0)
      df = self.get_adjusted_sample_names(df)
    elif self.exp_type =='counts':  
      df = pd.DataFrame.from_csv(self.counts, sep=('\t'),header=0)
      df = self.filter_low_count_genes(df)
      df = self.get_adjusted_sample_names(df)
    elif self.exp_type == 'logCPM':
      df = pd.DataFrame.from_csv(self.logCPM, sep=('\t'),header=0)
      df = self.get_adjusted_sample_names(df)
    passed = self.get_passed_QC()
    df = df.filter(passed, axis=1)
    return df
    
  
  def get_adjusted_sample_names(self,df):
    names = list(map(lambda x: x[:-3], df.columns.values))
    df.columns = names
    return df
  
  def get_passed_QC(self):
    df = pd.DataFrame.from_csv(self.QC, sep='\t', header=0)
    if self.strict != "N":
      df = df[df["Eisenberg Strict Score"] >= self.ehk]
    else:  
      df = df[df.eisenberg_score >= self.ehk]
    return list(df.index.values)
  
  def filter_low_count_genes(self,df):
    '''
     filter rows with mean counts less than read_count
    '''  
    df_counts = pd.DataFrame.from_csv(self.counts, sep=('\t'),header=0)
    count = df_counts.mean(axis=1,skipna=True)
    df_counts['count'] = count
    df_counts = df_counts[df_counts['count'] > self.read_count]
    new_df = df.filter(list(df_counts.index.values),axis=0)
    return new_df

  def filter_zero_rows_cols(self,df):
    '''
    remove rows with all zeroes and columns with all zeroes   
    '''
    df = df.dropna(how='all',axis=0)
    df = df.dropna(how='all',axis=1)
    #df = df.replace(np.nan,100000)
    df = df.replace(0.0,np.nan)
    df = df.dropna(how='all',axis=0)
    df = df.dropna(how='all',axis=1)
    df = df.replace(np.nan,0.0)
    return df
  
def option_parser():
  '''
  Required flags for getting gene expression data
  '''
  description = 'Get gene expression data'
  usage = "usage: %prog -o <output file name> [options] ALL OPTIONS ARE REQUIRED"
  parser = optparse.OptionParser(usage=usage, description=description)
  parser.add_option('-o', help='Path and name of output file', dest='o',
                    action='store', metavar="<output file>")
  parser.add_option('-g', help='path to a file with a list of genes (HUGO Names) one per line DEFAULT ALL', dest='g',
                    default="ALL",action='store', metavar="<gene list>")
  parser.add_option('-t', help='Type of data, counts or TPM or logCPM', dest='t',
                    action='store', metavar="<Data_type")
  parser.add_option('-e', help='Eisenburg threshold default = 7', dest='e',
                    default= 7,action='store', metavar="<EHK threshold>")
  parser.add_option('-E', help='Eisenburg strict threshold tcga, new or N  Default = N', dest='E',
                    default="N",action='store', metavar="<EHK Strict threshold>")
  parser.add_option('-p', help='Protein coding genes only? Y or N  Default = Y', dest='p',
                    default="Y",action='store', metavar="<protein coding>")
  parser.add_option('-i', help='Filter on indication e.g CRC choose [SRC KID BLAD HNSC GYN HEP MEL CRC PNET PDAC LUNG SI GALL]', dest='i',
                    default="",action='store', metavar="<indication>")
  parser.add_option('-c', help='Filter on compartment  e.g. live choose [ live  myeloid  treg  tcell  stroma  tumor ]', dest='c',
                    default="",action='store', metavar="<compartment>")
  parser.add_option('-m', help='Do you want to mean-centerand scale the data? Y or N Default N' , dest='m',
                    default="N",action='store', metavar="<mean_center>")
  parser.add_option('-s', help='Filter on list of samples', dest='s',
                     default="",action='store', metavar="<samples>")
  parser.add_option('-T', help='Transpose output to samples as rows: Y or N Default N', dest='T',
                     default="N",action='store', metavar="<transpose>")
  
  return parser


def initialize():
  # Get input arguments
  parser = option_parser()
  (opt, args) = parser.parse_args()
  outfile = opt.o
  gene_file = opt.g
  exp_type = opt.t
  ehk = int(opt.e)
  strict = opt.E
  pc = opt.p
  ind = opt.i
  comp = opt.c
  mcs = opt.m
  samples = opt.s
  transpose =opt.T
  if not (outfile and gene_file and exp_type):
    print '\n'
    parser.print_help()
    print '\n'
    sys.exit('Error: Not all required inputs have been provided')
    # print '\n'        
    # return
 
  return(outfile,gene_file,exp_type,ehk,strict,pc,ind,comp,mcs,samples,transpose)


def filter_zero_rows_cols(df):
    '''
    remove rows with all zeroes and columns with all zeroes   
    '''
    df = df.dropna(how='all',axis=0)
    df = df.dropna(how='all',axis=1)
    #df = df.replace(np.nan,100000)
    df = df.replace(0.0,np.nan)
    df = df.dropna(how='all',axis=0)
    df = df.dropna(how='all',axis=1)
    df = df.replace(np.nan,0.0)
    return df

def get_gene_list(f_name):
  '''
  Read the gene file and return
  gene_list = Hugo name sof genes
  ens_list = ensembl name sof genes
  '''
  h_e_dict = pickle.load(open("/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/RNAseq_tools/h_e_dict.pkl"))
  if f_name != 'ALL':
    gene_list = [line.rstrip() for line in open(f_name)]
  else:
    gene_list = h_e_dict.keys()
  ens_list = list(map(lambda x: h_e_dict[x],gene_list))
  return (ens_list,gene_list)

def get_sample_list(f_name):
  '''
  Read the sample file and return
  sample_list 
  '''
  return [line.rstrip() for line in open(f_name)]

def filter_pc_no_mt_genes(df):
  df_pc_genes = pd.read_csv('/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/RNAseq_tools/pc_genes_wo_mt.txt', sep='\t', header=0)
  df_pc_genes = df_pc_genes.squeeze()#make it a series
  return df.filter(df_pc_genes, axis='index')

def main():
  
  outfile,gene_file,exp_type,ehk,strict,pc,ind,comp,mcs,samples,transpose = initialize()
  ens_list, gene_list = get_gene_list(gene_file)
  qc_data = DataQC(exp_type,ehk,strict)
  df = qc_data.get_data()
  print df.shape
  if comp:
    df = df.filter(regex=comp, axis=1)
    print 'There are %d samples in the %s compartment that have an EHK > %d'%(df.shape[1],comp,ehk)
    print df.shape
  
  if ind:
    df = df.filter(regex=ind, axis=1)
    print 'There are %d samples in the %s indication'%(df.shape[1],ind)
    
  if pc == "Y":
    df = filter_pc_no_mt_genes(df)
    print df.shape
    print 'There are %d protein coding genes'%(df.shape[0])
  
  if samples:
    sample_list = get_sample_list(samples)
    
    df = df.filter(sample_list, axis=1)
  # print df.shape
  # df = df.sample(n=100,axis=0)
  # ens_list = list(df.index.values)
  # df = df.from_csv('/Users/bushrasamad/Documents/Immunoprofiler/RNAseq/Mar13_tcell_logCPM_ehk_9_10.tsv',sep='\t',header=0)
  # print df.shape
  # print df
  df = df.filter(ens_list, axis=0)
  # df = qc_data.filter_low_count_genes(df)
  print 'There are %d genes that you requested'% (df.shape[0])
  

  
  if mcs == 'Y':
    df= df.subtract(df.mean(axis=1,skipna=True), axis=0)
    df = df.divide(df.std(axis=1),axis=0)
    print 'Data is mean-centered and scaled'
  
  df = df.T.rename(index=str, columns= dict(zip(ens_list, gene_list)))
  #df = df.T
  print 'Printing..... %s'% outfile
  if transpose == "Y":
    df.to_csv(outfile,sep='\t',header=True)
  else:  
    df.T.to_csv(outfile,sep='\t',header=True)
  print'Done'

if __name__ == "__main__":
    main()  
