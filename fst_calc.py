#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:57:24 2022

@author: whirling-in-rags
"""

import os
import pandas as pd
from collections import Counter
import statistics
#import vcf
os.chdir('/home/whirling-in-rags/bioinfo/projects/dolezal_virus_2021/variant_calls/narna/dataframes/fst/')
#vcf_reader = vcf.Reader(open('BER5-E9.sai_sorted.vcf','r'))
CHROM = []
POS = []
REF = []
ALT = []
temp_list = []
INFO = []
depth = []
allele_freq = []
temp_list2 = []
temp_list3 = []
vcf_df = pd.read_csv("BER1-G4_fst.csv")
unique = pd.read_csv("unique_site_counts.csv")
unique_counts_df= pd.read_csv("unique_counts_df.csv")


unique_set_list = []
unique_set = unique.set_index('pos').to_dict()['count']
for i in unique_set:
    unique_set_list.append(i)

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
# =============================================================================
# pop_1 = pop_BER1G4.set_index('POS').to_dict()['num_alt_bases']
# pop_2 = pop_BER2A1.set_index('POS').to_dict()['num_alt_bases']
# pop_3 = pop_BER2A5.set_index('POS').to_dict()['num_alt_bases']
# pop_4 = pop_BER2C12.set_index('POS').to_dict()['num_alt_bases']
# pop_5 = pop_BER2F1.set_index('POS').to_dict()['num_alt_bases']
# pop_6 = pop_BER2F2.set_index('POS').to_dict()['num_alt_bases']
# pop_7 = pop_BER2G10.set_index('POS').to_dict()['num_alt_bases']
# pop_8 = pop_BER3A11.set_index('POS').to_dict()['num_alt_bases']
# pop_9 = pop_BER3A5.set_index('POS').to_dict()['num_alt_bases']
# pop_10 = pop_BER3B5.set_index('POS').to_dict()['num_alt_bases'] 
# pop_11 = pop_BER4A9.set_index('POS').to_dict()['num_alt_bases']
# pop_12 = pop_BER4B4.set_index('POS').to_dict()['num_alt_bases']
# pop_13 = pop_BER4B6.set_index('POS').to_dict()['num_alt_bases']
# pop_14 = pop_BER4C2.set_index('POS').to_dict()['num_alt_bases']
# pop_15 = pop_BER4D4.set_index('POS').to_dict()['num_alt_bases']
# pop_16 = pop_BER4D5.set_index('POS').to_dict()['num_alt_bases']
# pop_17 = pop_BER4H6.set_index('POS').to_dict()['num_alt_bases']
# pop_18 = pop_BER5E9.set_index('POS').to_dict()['num_alt_bases']
# =============================================================================




dataframe_list = [pop_BER1G4,pop_BER2A1,pop_BER2A5,pop_BER2C12,pop_BER2F1,pop_BER2F2,pop_BER2G10,\
                  pop_BER3A11,pop_BER3A5,pop_BER3B5,pop_BER4A9,pop_BER4B4,pop_BER4B6,\
                      pop_BER4C2,pop_BER4D4,pop_BER4D5,pop_BER4H6,pop_BER5E9]

iterator_test = 0
hs_list = []
hs_sum_list = []
for i in dataframe_list:
    pop_count = i["count"]
    pop_call = i["call_count"]
    pop_a_freq = i["call_freq_a"]
    pop_c_freq = i["call_freq_c"]
    pop_t_freq = i["call_freq_t"]
    pop_a_freq = i["call_freq_a"]
    pop_allele_freq = i["allele_freq"]
    pop_numerator = (pop_count/(pop_count-1))
    pop_call_squared = pop_allele_freq **2
    pop_ref_squared = (1-pop_allele_freq)**2
    pop_freq_sum = pop_call_squared+pop_ref_squared
    pop_freq_sum_dff = (1-pop_freq_sum)
    pop_numerator_over_freq_sum_dff = (pop_numerator*pop_freq_sum_dff)
    hs_list.append(pop_numerator_over_freq_sum_dff)
for i in hs_list:
    hs_sum_list.append(sum(i))

global_hs = sum(hs_sum_list)/len(hs_sum_list)
global_ht = 8.582551779843925 #from another calculation, thanks poorva
average_pairwise_fst = 0.3228934996979474 #from another calculation, pairwise
pairwise_fst_list = [0.092932113317148, 0.3296526606044559, 0.0987222156685444, 0.009897595806228593, 0.19483658065429818, 0.28008959317718696, 0.3446474789117823, 0.23957222306602102, 0.25880744620920654, 0.3903412808749252, 0.3381415604000391, 0.16145021793960976, 0.2916380111543695, 0.3458832210900801, 0.12819988176587224, 0.11837391125195433, 0.24818942893918147, 0.3210421208636089, 0.13886042010359587, 0.048847691038408275, 0.1960224195706627, 0.23366696068056614, 0.3439644025144304, 0.24644221043372722, 0.2827531983538485, 0.3538503077483568, 0.3296602347105665, 0.14466759336505833, 0.2502364580457452, 0.3243522340204993, 0.13239550962456656, 0.09567463276160414, 0.2774845347833843, 0.38855925182565376, 0.3428624196164213, 0.3697349979255497, 0.4305450727437385, 0.4877122438520943, 0.33176838613188003, 0.4633354526242613, 0.37005830613051816, 0.36843803402469966, 0.32443447383702556, 0.484227781922485, 0.4188603516655711, 0.3875096190416454, 0.363308171080791, 0.4547947221244085, 0.1091577998737632, 0.24380421050407486, 0.32877981827727587, 0.39826457830437884, 0.2986256076634695, 0.3344478762320096, 0.4475944971403791, 0.38935449681693185, 0.22330954624979715, 0.37095600941535345, 0.31572111059813257, 0.21407439780598744, 0.20421095431007663, 0.3544875953843029, 0.20171333835730698, 0.2837358533227154, 0.360990423968229, 0.23480668376847302, 0.29574148238933057, 0.4233503866325697, 0.3483000592677442, 0.16880686496367367, 0.34671574332175137, 0.36596743987437336, 0.1513200126861356, 0.14658725022843724, 0.3070762869385577, 0.2851420629016528, 0.34065268148493066, 0.3057587371593589, 0.2933206129211669, 0.3843474942407187, 0.3462460383827534, 0.17320898078695482, 0.3319689217170865, 0.3501626326589936, 0.23557685338340922, 0.22049164249396624, 0.3134712185746332, 0.34941022581437764, 0.32586294791965953, 0.32534722062986055, 0.42230921669692617, 0.41159940074468093, 0.24986867849121494, 0.31374623494423937, 0.4102447136089093, 0.3025182163317942, 0.26901114538915233, 0.28828060393143, 0.3580019695354513, 0.37056843678841567, 0.48480574861888187, 0.45907711648328614, 0.3146239866469216, 0.4529537562521888, 0.46962451309572967, 0.3756469767236306, 0.35209156202139497, 0.3744640740224604, 0.3948696882128112, 0.36971799713884007, 0.3448109445186599, 0.26826284592450406, 0.40186657477121296, 0.3476522356685863, 0.29549685618306276, 0.2806468809085453, 0.2953291728733297, 0.47867789138389555, 0.44580835464078833, 0.2964127771697796, 0.4366419357666809, 0.4548649310830558, 0.3286813356585746, 0.31751630364512096, 0.4069177093648291, 0.43481575960684704, 0.37257349951529045, 0.5187835412564287, 0.44909440382903865, 0.43959343964122116, 0.4245024168838011, 0.4874899047100945, 0.31240093975721317, 0.4546207537047027, 0.4031048925244095, 0.37546682249534424, 0.35671573204330753, 0.4432462913090591, 0.32515505722784577, 0.3320977168449243, 0.21689920600530402, 0.1918456628474087, 0.2884145547472158, 0.48281000301789406, 0.38477636096557993, 0.3743296943712391, 0.34436641703813137, 0.3946195328023129, 0.3748465888762661, 0.45580709214413295, 0.21642148722707807, 0.3456172481272222, 0.3323654196666123]
hs_distribution_by_pop = []
global_fst = (global_ht-global_hs)/global_ht



# =============================================================================
#     pop_1_numerator_list.append(pop_1_count[iterator_calc]/(pop_1_count[iterator_calc]-1))
#     pop_1_call_squared_list.append(pop_1_allele_freq[iterator_calc]**2)
#     pop_1_ref_squared_list.append((1-pop_1_allele_freq[iterator_calc])**2)
#     pop_1_freq_sum_list.append(pop_1_call_squared_list[iterator_calc]+pop_1_ref_squared_list[iterator_calc])
#     pop_1_freq_sum_dff_list.append(1-pop_1_freq_sum_list[iterator_calc])
#     pop_1_numerator_over_freq_sum_dff_list.append(pop_1_numerator_list[iterator_calc]*pop_1_freq_sum_dff_list[iterator_calc])
#     pop_1_hs_list.append(pop_1_numerator_over_freq_sum_dff_list[iterator_calc])
# =============================================================================
    


# =============================================================================
# while iterator_test < len(dataframe_list):
#     pop_1_count = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'count'].iloc[0])
#     pop_1_call = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'call_count'].iloc[0])
#     pop_1_a_freq = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'call_freq_a'].iloc[0])
#     pop_1_c_freq = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'call_freq_c'].iloc[0])
#     pop_1_t_freq = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'call_freq_t'].iloc[0])
#     pop_1_g_freq = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'call_freq_g'].iloc[0])
#     pop_1_allele_freq = (dataframe_list[iterator_test].loc[dataframe_list[iterator_test]['pos']==j,'allele_freq'].iloc[0])
#     iterator_test = iterator_test + 1
# =============================================================================





pos_iterator = 0

while pos_iterator < len(vcf_df['POS']):
    n_over_n_minus_one = ((vcf_df['DEPTH'][pos_iterator])/(vcf_df['DEPTH'][pos_iterator]-1))
    alt_freq_squared = ((vcf_df['ALLELE_FREQ'][pos_iterator])**2)
    ref_freq_squared = (((1-vcf_df['ALLELE_FREQ'][pos_iterator])**2))
    freq_sum = (alt_freq_squared + ref_freq_squared)
    freq_sum_diff = (1-freq_sum)
    n_over_freq_sum_diff = (n_over_n_minus_one*freq_sum_diff)
    hs = (n_over_freq_sum_diff)
    hs_list.append(hs)
    pos_iterator += 1

    
#test_0 = (vcf_df['DEPTH']/(vcf_df['DEPTH']-1))
#test_1 = (vcf_df['ALLELE_FREQ']**2) 
#test_2 = ((1-vcf_df['ALLELE_FREQ'])**2)
#test_3 = (test_1 + test_2)
#test_4 = (1 - test_3)
#test_5 = (test_0 * test_4)
#test_6 = (sum(test_5)/len(test_5))




os.chdir('/home/whirling-in-rags/bioinfo/projects/dolezal_virus_2021/variant_calls/narna/dataframes/')
big_df = pd.read_csv("all.csv")
pos_list = []
j_list = []
unique_j_list = []
overlap_list = []
overlap_counts = []
populations = big_df.groupby('population')
big_df_pos = populations.apply(lambda x: x['POS'].unique())

for i in big_df_pos:
    pos_list.append(list(i))
    
for i in pos_list:
    for j in i:
        j_list.append(j)
unique_j_list = list(set(j_list))
flat_pos_list = [j for i in pos_list for j in i]
iterator_1 = 0

while iterator_1 < len(unique_j_list):
    if flat_pos_list.count(unique_j_list[iterator_1]) >= 1:
        overlap_counts.append(flat_pos_list.count(unique_j_list[iterator_1]))
        overlap_list.append(unique_j_list[iterator_1])
    iterator_1 = iterator_1 + 1

#os.chdir('/home/whirling-in-rags/bioinfo/projects/dolezal_virus_2021/variant_calls/narna/dataframes/depths/')

#read_depths_df = pd.read_csv("read_depths.csv")

#y = read_depths_df[read_depths_df['pos'].isin(unique_j_list)]
#x = y.groupby("pos")["count"].sum()


#x.to_csv("unique_site_counts.csv")
