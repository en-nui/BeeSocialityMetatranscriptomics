#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 14:22:41 2022

@author: whirling-in-rags
"""

import os
import itertools
import pandas as pd
import csv
#import vcf
os.chdir('/home/whirling-in-rags/bioinfo/projects/dolezal_virus_2021/variant_calls/narna/dataframes/fst/')

#dataframes to work with
all_snp_df = pd.read_csv("all.csv")
all_reads_df = pd.read_csv("read_depths.csv")
unique_sites_df = pd.read_csv("unique_site_counts.csv")
variant_sites_by_pop = pd.read_csv("variant_sites.csv")
unique_counts_df= pd.read_csv("unique_counts_df.csv")
#end dataframes to work with

#unique_variant_list
unique_variant_list = []
for i in variant_sites_by_pop["0"]:
    unique_variant_list.append(i)
newline_strip = [i.replace("\n",'') for i in unique_variant_list]
newline_split = []
for i in newline_strip:
    i = i.split()
    newline_split.append(i)
    del i[0]
iterator_a = 0
while iterator_a < len(newline_split):
    del newline_split[iterator_a][-1]
    iterator_a += 1
final_variant_list = [[int(float(j)) for j in i] for i in newline_split]
#end unique_variant_list
pair_order_list = [tuple(i) for i in final_variant_list]
pair_order_list = itertools.permutations(final_variant_list,2)
pair_order_list = list(pair_order_list)
unique_pairwise_combinations = []
for i in pair_order_list:
    a = set(i[0])
    b = set(i[1])
    c = list((set(a) | set(b) - set(a) & set(b)))
    unique_pairwise_combinations.append(c)
    
#pandas stuff
pair_order_list2 = []
pair_order_list3 = []
test_list = []
for i in unique_counts_df["pos"]:
    test_list.append(i)
    
despair = all_reads_df[all_reads_df["pos"].isin(test_list)]
holder_list2 = []
iterator_base = 0
iterator_adj = 1
while iterator_base < (len(final_variant_list)):
         iterator_adj = 1 + iterator_base
         while iterator_adj < len(final_variant_list):
             x = [final_variant_list[(iterator_base)],final_variant_list[(iterator_adj)]]
             pair_order_list2.append(x)
             
             iterator_adj += 1
         iterator_base += 1
    
for i in pair_order_list2:
    a = set(i[0])
    b = set(i[1])
    c = list((set(a) | set(b) - set(a) & set(b)))
    pair_order_list3.append(c)
    
pop_BER1G4 = unique_counts_df[(unique_counts_df['population'] == 'BER1-G4')]
pop_BER2A1 = unique_counts_df[(unique_counts_df['population'] == 'BER2-A1')]
pop_BER2A5 = unique_counts_df[(unique_counts_df['population'] == 'BER2-A5')]
pop_BER2C12 = unique_counts_df[(unique_counts_df['population'] == 'BER2-C12')]
pop_BER2F1 = unique_counts_df[(unique_counts_df['population'] == 'BER2-F1')]
pop_BER2F2 = unique_counts_df[(unique_counts_df['population'] == 'BER2-F2')]
pop_BER2G10 = unique_counts_df[(unique_counts_df['population'] == 'BER2-G10')]
pop_BER3A11 = unique_counts_df[(unique_counts_df['population'] == 'BER3-A11')]
pop_BER3A5 = unique_counts_df[(unique_counts_df['population'] == 'BER3-A5')]
pop_BER3B5 = unique_counts_df[(unique_counts_df['population'] == 'BER3-B5')]
pop_BER4A9 = unique_counts_df[(unique_counts_df['population'] == 'BER4-A9')]
pop_BER4B4 = unique_counts_df[(unique_counts_df['population'] == 'BER4-B4')]
pop_BER4B6 = unique_counts_df[(unique_counts_df['population'] == 'BER4-B6')]
pop_BER4C2 = unique_counts_df[(unique_counts_df['population'] == 'BER4-C2')]
pop_BER4D4 = unique_counts_df[(unique_counts_df['population'] == 'BER4-D4')]
pop_BER4D5 = unique_counts_df[(unique_counts_df['population'] == 'BER4-D5')]
pop_BER4H6 = unique_counts_df[(unique_counts_df['population'] == 'BER4-H6')]
pop_BER5E9 = unique_counts_df[(unique_counts_df['population'] == 'BER5-E9')]

dataframe_list = [pop_BER1G4,pop_BER2A1,pop_BER2A5,pop_BER2C12,pop_BER2F1,pop_BER2F2,pop_BER2G10,\
                  pop_BER3A11,pop_BER3A5,pop_BER3B5,pop_BER4A9,pop_BER4B4,pop_BER4B6,\
                      pop_BER4C2,pop_BER4D4,pop_BER4D5,pop_BER4H6,pop_BER5E9]

# =============================================================================
# dataframe_test = [pop_BER1G4,pop_BER2A1]    
# tester_1 = unique_pairwise_combinations[0]
pop_1_pos = []
pop_1_count = []
pop_1_call = []
pop_1_allele_freq = []
pop_1_a_freq = []
pop_1_c_freq = []
pop_1_t_freq = []
pop_1_g_freq = []
# 
pop_2_pos = []
pop_2_count = []
pop_2_call = []
pop_2_allele_freq = []
pop_2_a_freq = []
pop_2_c_freq = []
pop_2_t_freq = []
pop_2_g_freq = []
depth_list = []
# 
# fst_list = []
# iterator_test = 0
# 
# 
# 
# while iterator_test < len(dataframe_test):
#     for i in dataframe_test[iterator_test]["pos"]:
#         if i in (tester_1):
#             pop_1_pos.append(i)
#     for j in pop_1_pos:
#         pop_1_count.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'count'].iloc[0])
#         pop_1_call.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_count'].iloc[0])
#         pop_1_a_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_a'].iloc[0])
#         pop_1_c_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_c'].iloc[0])
#         pop_1_t_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_t'].iloc[0])
#         pop_1_g_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_g'].iloc[0])
#         pop_1_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
#         
#     iterator_test += 1
#     for i in dataframe_test[iterator_test]["pos"]:
#         if i in (tester_1):
#             pop_2_pos.append(i)
#     for j in pop_2_pos:
#         pop_2_count.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'count'].iloc[0])
#         pop_2_call.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_count'].iloc[0])
#         pop_2_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
#         pop_2_a_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_a'].iloc[0])
#         pop_2_c_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_c'].iloc[0])
#         pop_2_t_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_t'].iloc[0])
#         pop_2_g_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_g'].iloc[0])
#     iterator_test += 1
# 
# 
pop_1_numerator_list = []
pop_1_call_squared_list = []
pop_1_ref_squared_list = []
pop_1_freq_sum_list = []
pop_1_freq_sum_dff_list = []
pop_1_numerator_over_freq_sum_dff_list = []
pop_1_hs_list = []    
pop_2_numerator_list = []
pop_2_call_squared_list = []
pop_2_ref_squared_list = []
pop_2_freq_sum_list = []
pop_2_freq_sum_dff_list = []
pop_2_numerator_over_freq_sum_dff_list = []
pop_2_hs_list = []    
pop_total_numerator_list = []
pop_total_sum_list = []
pop_total_allele_sum_list = []
pop_total_allele_freq_list = []
pop_total_call_squared_list = []
pop_total_ref_squared_list = []
pop_total_freq_sum_list = []
pop_total_freq_sum_dff_list = []
pop_total_numerator_over_freq_sum_dff_list = []
pop_total_ht_list = []  
hs_list = []
iterator_calc = 0
# 
# 
# while iterator_calc < len(pop_1_count):
#     #pop1_h
#     pop_1_numerator_list.append(pop_1_count[iterator_calc]/(pop_1_count[iterator_calc]-1))
#     pop_1_call_squared_list.append(pop_1_allele_freq[iterator_calc]**2)
#     pop_1_ref_squared_list.append((1-pop_1_allele_freq[iterator_calc])**2)
#     pop_1_freq_sum_list.append(pop_1_call_squared_list[iterator_calc]+pop_1_ref_squared_list[iterator_calc])
#     pop_1_freq_sum_dff_list.append(1-pop_1_freq_sum_list[iterator_calc])
#     pop_1_numerator_over_freq_sum_dff_list.append(pop_1_numerator_list[iterator_calc]*pop_1_freq_sum_dff_list[iterator_calc])
#     pop_1_hs_list.append(pop_1_numerator_over_freq_sum_dff_list[iterator_calc])
#     #pop2_h
#     
#     pop_2_numerator_list.append(pop_2_count[iterator_calc]/(pop_2_count[iterator_calc]-1))
#     pop_2_call_squared_list.append(pop_2_allele_freq[iterator_calc]**2)
#     pop_2_ref_squared_list.append((1-pop_2_allele_freq[iterator_calc])**2)
#     pop_2_freq_sum_list.append(pop_2_call_squared_list[iterator_calc]+pop_2_ref_squared_list[iterator_calc])
#     pop_2_freq_sum_dff_list.append(1-pop_2_freq_sum_list[iterator_calc])
#     pop_2_numerator_over_freq_sum_dff_list.append(pop_2_numerator_list[iterator_calc]*pop_2_freq_sum_dff_list[iterator_calc])
#     pop_2_hs_list.append(pop_2_numerator_over_freq_sum_dff_list[iterator_calc])
#     
#     #pop_total
#     pop_total_numerator_list.append((pop_1_count[iterator_calc]+pop_2_count[iterator_calc])/((pop_1_count[iterator_calc]-1)+(pop_2_count[iterator_calc]-1)))
#     pop_total_sum_list.append(pop_1_count[iterator_calc]+pop_2_count[iterator_calc])
#     pop_total_allele_sum_list.append(pop_1_call[iterator_calc]+pop_2_call[iterator_calc])
#     pop_total_allele_freq_list.append(pop_total_allele_sum_list[iterator_calc]/pop_total_sum_list[iterator_calc])
#     pop_total_call_squared_list.append(pop_total_allele_freq_list[iterator_calc]**2)
#     pop_total_ref_squared_list.append((1-pop_total_allele_freq_list[iterator_calc])**2)
#     pop_total_freq_sum_list.append(pop_total_call_squared_list[iterator_calc]+pop_total_ref_squared_list[iterator_calc])
#     pop_total_freq_sum_dff_list.append(1-pop_total_freq_sum_list[iterator_calc])
#     pop_total_numerator_over_freq_sum_dff_list.append(pop_total_numerator_list[iterator_calc]*pop_total_freq_sum_dff_list[iterator_calc])
#     pop_total_ht_list.append(pop_total_numerator_over_freq_sum_dff_list[iterator_calc])
#     
#     
#     #hs
#     hs_list.append((pop_1_hs_list[iterator_calc]+pop_2_hs_list[iterator_calc])/2)
#     
#     iterator_calc += 1
#     #fst for pairwise comparison
# ht = sum(pop_total_ht_list)
# hs = sum(hs_list)
# fst = (ht-hs)/ht
# =============================================================================
    

### adding snp unique to pop comparison
pop_1_pos = []
pop_2_pos = []
iterator_pairwise = 0
iterator_adj = 1
iterator_base = 0
iterator_test = 0
fst_list = []
while iterator_pairwise < len(pair_order_list3):

    while iterator_base < len(dataframe_list):
        iterator_adj = 1 + iterator_base
        while iterator_adj < len(dataframe_list):
            dataframe_test = [dataframe_list[(iterator_base)],dataframe_list[(iterator_adj)]]
            
            while iterator_test < len(dataframe_test):

                for i in dataframe_test[iterator_test]["pos"]:
                    if i in (pair_order_list3[iterator_pairwise]):
                        pop_1_pos.append(i)
                for j in pop_1_pos:
                    pop_1_count.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'count'].iloc[0])
                    pop_1_call.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_count'].iloc[0])
                    pop_1_a_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_a'].iloc[0])
                    pop_1_c_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_c'].iloc[0])
                    pop_1_t_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_t'].iloc[0])
                    pop_1_g_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_g'].iloc[0])
                    pop_1_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
                    
                iterator_test += 1
                for i in dataframe_test[iterator_test]["pos"]:
                    if i in (pair_order_list3[iterator_pairwise]):
                        pop_2_pos.append(i)
                for j in pop_2_pos:
                    pop_2_count.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'count'].iloc[0])
                    pop_2_call.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_count'].iloc[0])
                    pop_2_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
                    pop_2_a_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_a'].iloc[0])
                    pop_2_c_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_c'].iloc[0])
                    pop_2_t_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_t'].iloc[0])
                    pop_2_g_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'call_freq_g'].iloc[0])
                iterator_test += 1

    
            while iterator_calc < len(pop_1_count):
                #pop1_h
                pop_1_numerator_list.append(pop_1_count[iterator_calc]/(pop_1_count[iterator_calc]-1))
                pop_1_call_squared_list.append(pop_1_allele_freq[iterator_calc]**2)
                pop_1_ref_squared_list.append((1-pop_1_allele_freq[iterator_calc])**2)
                pop_1_freq_sum_list.append(pop_1_call_squared_list[iterator_calc]+pop_1_ref_squared_list[iterator_calc])
                pop_1_freq_sum_dff_list.append(1-pop_1_freq_sum_list[iterator_calc])
                pop_1_numerator_over_freq_sum_dff_list.append(pop_1_numerator_list[iterator_calc]*pop_1_freq_sum_dff_list[iterator_calc])
                pop_1_hs_list.append(pop_1_numerator_over_freq_sum_dff_list[iterator_calc])
                #pop2_h
                
                pop_2_numerator_list.append(pop_2_count[iterator_calc]/(pop_2_count[iterator_calc]-1))
                pop_2_call_squared_list.append(pop_2_allele_freq[iterator_calc]**2)
                pop_2_ref_squared_list.append((1-pop_2_allele_freq[iterator_calc])**2)
                pop_2_freq_sum_list.append(pop_2_call_squared_list[iterator_calc]+pop_2_ref_squared_list[iterator_calc])
                pop_2_freq_sum_dff_list.append(1-pop_2_freq_sum_list[iterator_calc])
                pop_2_numerator_over_freq_sum_dff_list.append(pop_2_numerator_list[iterator_calc]*pop_2_freq_sum_dff_list[iterator_calc])
                pop_2_hs_list.append(pop_2_numerator_over_freq_sum_dff_list[iterator_calc])
                
                #pop_total
                pop_total_numerator_list.append((pop_1_count[iterator_calc]+pop_2_count[iterator_calc])/((pop_1_count[iterator_calc]-1)+(pop_2_count[iterator_calc]-1)))
                pop_total_sum_list.append(pop_1_count[iterator_calc]+pop_2_count[iterator_calc])
                pop_total_allele_sum_list.append(pop_1_call[iterator_calc]+pop_2_call[iterator_calc])
                pop_total_allele_freq_list.append(pop_total_allele_sum_list[iterator_calc]/pop_total_sum_list[iterator_calc])
                pop_total_call_squared_list.append(pop_total_allele_freq_list[iterator_calc]**2)
                pop_total_ref_squared_list.append((1-pop_total_allele_freq_list[iterator_calc])**2)
                pop_total_freq_sum_list.append(pop_total_call_squared_list[iterator_calc]+pop_total_ref_squared_list[iterator_calc])
                pop_total_freq_sum_dff_list.append(1-pop_total_freq_sum_list[iterator_calc])
                pop_total_numerator_over_freq_sum_dff_list.append(pop_total_numerator_list[iterator_calc]*pop_total_freq_sum_dff_list[iterator_calc])
                pop_total_ht_list.append(pop_total_numerator_over_freq_sum_dff_list[iterator_calc])
                
                
                #hs
                hs_list.append((pop_1_hs_list[iterator_calc]+pop_2_hs_list[iterator_calc])/2)
                iterator_calc += 1
            #fst for pairwise comparison
            iterator_calc = 0
            ht = sum(pop_total_ht_list)
            hs = sum(hs_list)
            fst_list.append((ht-hs)/ht)
            pop_1_numerator_list = []
            pop_1_call_squared_list = []
            pop_1_ref_squared_list = []
            pop_1_freq_sum_list = []
            pop_1_freq_sum_dff_list = []
            pop_1_numerator_over_freq_sum_dff_list = []
            pop_1_hs_list = []    
            pop_2_numerator_list = []
            pop_2_call_squared_list = []
            pop_2_ref_squared_list = []
            pop_2_freq_sum_list = []
            pop_2_freq_sum_dff_list = []
            pop_2_numerator_over_freq_sum_dff_list = []
            pop_2_hs_list = []    
            pop_total_numerator_list = []
            pop_total_sum_list = []
            pop_total_allele_sum_list = []
            pop_total_allele_freq_list = []
            pop_total_call_squared_list = []
            pop_total_ref_squared_list = []
            pop_total_freq_sum_list = []
            pop_total_freq_sum_dff_list = []
            pop_total_numerator_over_freq_sum_dff_list = []
            pop_total_ht_list = []  
            hs_list = []
            pop_1_pos = []
            pop_1_count = []
            pop_1_call = []
            pop_1_allele_freq = []
            pop_1_a_freq = []
            pop_1_c_freq = []
            pop_1_t_freq = []
            pop_1_g_freq = []
            pop_2_pos = []
            pop_2_count = []
            pop_2_call = []
            pop_2_allele_freq = []
            pop_2_a_freq = []
            pop_2_c_freq = []
            pop_2_t_freq = []
            pop_2_g_freq = []
            depth_list = []
            iterator_calc = 0
            iterator_test = 0
            iterator_adj += 1
            iterator_pairwise += 1
            
            
            
            
            
        iterator_base += 1
        
print(fst_list)
print((sum(fst_list))/(len(fst_list)))

 