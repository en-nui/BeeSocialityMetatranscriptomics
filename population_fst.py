#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:55:47 2022

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
pop_BER4B4 = unique_counts_df[(unique_counts_df['population'] == 'BER4-B4')]#actually ber5b4
pop_BER4B6 = unique_counts_df[(unique_counts_df['population'] == 'BER4-B6')]
pop_BER4C2 = unique_counts_df[(unique_counts_df['population'] == 'BER4-C2')]
pop_BER4D4 = unique_counts_df[(unique_counts_df['population'] == 'BER4-D4')]
pop_BER4D5 = unique_counts_df[(unique_counts_df['population'] == 'BER4-D5')]
pop_BER4H6 = unique_counts_df[(unique_counts_df['population'] == 'BER4-H6')]
pop_BER5E9 = unique_counts_df[(unique_counts_df['population'] == 'BER5-E9')]
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

COB = [pop_BER1G4,pop_BER2F1,pop_BER2F2]
cob_pos_list = []
cob_count_list = []
cob_call_list = []
cob_call_count_list = []
WTH = [pop_BER2A1,pop_BER3A5,pop_BER3B5]
wth_pos_list = []
wth_count_list = []
wth_call_list = []
wth_call_count_list = []
SOS = [pop_BER2A5,pop_BER2C12,pop_BER4C2,pop_BER4D5,pop_BER4H6,pop_BER4B4]
sos_pos_list = []
sos_count_list = []
sos_call_list = []
sos_call_count_list = []
PEE = [pop_BER2G10,pop_BER3A11,pop_BER4D4]
pee_pos_list = []
pee_count_list = []
pee_call_list = []
pee_call_count_list = []
CRN = [pop_BER4A9,pop_BER4B6,pop_BER5E9]
crn_pos_list = []
crn_count_list = []
crn_call_list = []
crn_call_count_list = []

pops = [COB,WTH,SOS,PEE,CRN]
iterator_z = 0
iterator_y = 0

for i in COB:
    while iterator_y < len(i["pos"]):
        cob_pos_list.append(i["pos"].iloc[iterator_y])
        cob_count_list.append(i["count"].iloc[iterator_y])
        cob_call_list.append(i["call"].iloc[iterator_y])
        cob_call_count_list.append(i["call_count"].iloc[iterator_y])
        iterator_y +=1
    iterator_y = 0       

for i in WTH:
    while iterator_y < len(i["pos"]):
        wth_pos_list.append(i["pos"].iloc[iterator_y])
        wth_count_list.append(i["count"].iloc[iterator_y])
        wth_call_list.append(i["call"].iloc[iterator_y])
        wth_call_count_list.append(i["call_count"].iloc[iterator_y])
        iterator_y +=1
    iterator_y = 0

for i in SOS:
    while iterator_y < len(i["pos"]):
        sos_pos_list.append(i["pos"].iloc[iterator_y])
        sos_count_list.append(i["count"].iloc[iterator_y])
        sos_call_list.append(i["call"].iloc[iterator_y])
        sos_call_count_list.append(i["call_count"].iloc[iterator_y])
        iterator_y +=1
    iterator_y = 0       


for i in PEE:
    while iterator_y < len(i["pos"]):
        pee_pos_list.append(i["pos"].iloc[iterator_y])
        pee_count_list.append(i["count"].iloc[iterator_y])
        pee_call_list.append(i["call"].iloc[iterator_y])
        pee_call_count_list.append(i["call_count"].iloc[iterator_y])
        iterator_y +=1
    iterator_y = 0       

iterator_y = 0       
for i in CRN:
    while iterator_y < len(i["pos"]):
        crn_pos_list.append(i["pos"].iloc[iterator_y])
        crn_count_list.append(i["count"].iloc[iterator_y])
        crn_call_list.append(i["call"].iloc[iterator_y])
        crn_call_count_list.append(i["call_count"].iloc[iterator_y])
        iterator_y +=1
    iterator_y = 0       

cob_pos_three = []
cob_count_three = []
cob_call_three = []
cob_call_count_three = []
for i in range(0,len(cob_pos_list),495):
    x = i
    cob_pos_three.append(cob_pos_list[x:x+495])
    cob_count_three.append(cob_count_list[x:x+495])
    cob_call_three.append(cob_call_list[x:x+495])
    cob_call_count_three.append(cob_call_count_list[x:x+495])

wth_pos_three = []
wth_count_three = []
wth_call_three = []
wth_call_count_three = []
for i in range(0,len(cob_pos_list),495):
    x = i
    wth_pos_three.append(wth_pos_list[x:x+495])
    wth_count_three.append(wth_count_list[x:x+495])
    wth_call_three.append(wth_call_list[x:x+495])
    wth_call_count_three.append(wth_call_count_list[x:x+495])
sos_pos_three = []
sos_count_three = []
sos_call_three = []
sos_call_count_three = []   
for i in range(0,len(sos_pos_list),495):
    x = i
    sos_pos_three.append(sos_pos_list[x:x+495])
    sos_count_three.append(sos_count_list[x:x+495])
    sos_call_three.append(sos_call_list[x:x+495])
    sos_call_count_three.append(sos_call_count_list[x:x+495])
pee_pos_three = []
pee_count_three = []
pee_call_three = []
pee_call_count_three = []
for i in range(0,len(cob_pos_list),495):
    x = i
    pee_pos_three.append(pee_pos_list[x:x+495])
    pee_count_three.append(pee_count_list[x:x+495])
    pee_call_three.append(pee_call_list[x:x+495])
    pee_call_count_three.append(pee_call_count_list[x:x+495])
crn_pos_three = []
crn_count_three = []
crn_call_three = []
crn_call_count_three = []
for i in range(0,len(cob_pos_list),495):
    x = i
    crn_pos_three.append(crn_pos_list[x:x+495])
    crn_count_three.append(crn_count_list[x:x+495])
    crn_call_three.append(crn_call_list[x:x+495])
    crn_call_count_three.append(crn_call_count_list[x:x+495])
    
cob_counts = []
wth_counts = []
sos_counts = []
pee_counts = []
crn_counts = []

cob_calls = []
wth_calls = []
sos_calls = []
pee_calls = []
crn_calls = []

cob_call_counts = []
wth_call_counts = []
sos_call_counts = []
pee_call_counts = []
crn_call_counts = []
iterator_x = 0
iterator_y = 0
while iterator_x < len(cob_pos_three[0]):
    while iterator_y < len(cob_pos_three):
        cob_counts.append(cob_count_three[iterator_y][iterator_x])
        wth_counts.append(wth_count_three[iterator_y][iterator_x])
        pee_counts.append(pee_count_three[iterator_y][iterator_x])
        crn_counts.append(crn_count_three[iterator_y][iterator_x])
        cob_calls.append(cob_call_three[iterator_y][iterator_x])
        wth_calls.append(wth_call_three[iterator_y][iterator_x])
        pee_calls.append(pee_call_three[iterator_y][iterator_x])
        crn_calls.append(crn_call_three[iterator_y][iterator_x])
        cob_call_counts.append(cob_call_count_three[iterator_y][iterator_x])
        wth_call_counts.append(wth_call_count_three[iterator_y][iterator_x])
        pee_call_counts.append(pee_call_count_three[iterator_y][iterator_x])
        crn_call_counts.append(crn_call_count_three[iterator_y][iterator_x])
        iterator_y += 1
    iterator_x += 1
    iterator_y = 0
iterator_x = 0
while iterator_x < len(sos_pos_three[0]):
    while iterator_y < len(sos_pos_three):
        sos_counts.append(sos_count_three[iterator_y][iterator_x])
        sos_calls.append(sos_call_three[iterator_y][iterator_x])
        sos_call_counts.append(sos_call_count_three[iterator_y][iterator_x])
        iterator_y +=1
    iterator_x +=1
    iterator_y = 0
cob_count = []
wth_count = []
sos_count = []
pee_count = []
crn_count = []

cob_call = []
wth_call = []
sos_call = []
pee_call = []
crn_call = []

cob_call_count = []
wth_call_count = []
sos_call_count = []
pee_call_count = []
crn_call_count = []



for i in range(0,len(cob_counts),3):
    x = i
    cob_count.append(cob_counts[x:x+3])
    wth_count.append(wth_counts[x:x+3])
    pee_count.append(pee_counts[x:x+3])
    crn_count.append(crn_counts[x:x+3])
    
    cob_call.append(cob_calls[x:x+3])
    wth_call.append(wth_calls[x:x+3])
    pee_call.append(pee_calls[x:x+3])
    crn_call.append(crn_calls[x:x+3])
    
    cob_call_count.append(cob_call_counts[x:x+3])
    wth_call_count.append(wth_call_counts[x:x+3])
    pee_call_count.append(pee_call_counts[x:x+3])
    crn_call_count.append(crn_call_counts[x:x+3])
    
for i in range(0,len(sos_counts),6):
    x = i
    sos_count.append(sos_counts[x:x+6])
    sos_call.append(sos_calls[x:x+6])
    sos_call_count.append(sos_call_counts[x:x+6])
cob_count_sum = []
wth_count_sum = []
sos_count_sum = []
pee_count_sum = []
crn_count_sum = []

for i in cob_count:
    cob_count_sum.append(sum(i))
for i in wth_count:
    wth_count_sum.append(sum(i))
for i in sos_count:
    sos_count_sum.append(sum(i))
for i in pee_count:
    pee_count_sum.append(sum(i))
for i in crn_count:
    crn_count_sum.append(sum(i))

count_sums = [cob_count_sum,wth_count_sum,sos_count_sum,pee_count_sum,crn_count_sum]

cob_call_freqs = []
wth_call_freqs = []
sos_call_freqs = []
pee_call_freqs = []
crn_call_freqs = []
iterator_w = 0
iterator_v = 0

while iterator_w < len(cob_call_count):
    while iterator_v < len(cob_call_count[iterator_w]):
        cob_call_freqs.append(cob_call_count[iterator_w][iterator_v]/cob_count_sum[iterator_w])
        wth_call_freqs.append(wth_call_count[iterator_w][iterator_v]/wth_count_sum[iterator_w])
        pee_call_freqs.append(pee_call_count[iterator_w][iterator_v]/pee_count_sum[iterator_w])
        crn_call_freqs.append(crn_call_count[iterator_w][iterator_v]/crn_count_sum[iterator_w])
        iterator_v += 1
    iterator_w += 1
    iterator_v = 0
iterator_w = 0
while iterator_w < len(sos_call_count):
    while iterator_v < len(sos_call_count[iterator_w]):
        sos_call_freqs.append(sos_call_count[iterator_w][iterator_v]/sos_count_sum[iterator_w])
        iterator_v += 1
    iterator_w += 1
    iterator_v = 0
cob_call_freq = []
wth_call_freq = []
sos_call_freq = []
pee_call_freq = []
crn_call_freq = []   
for i in range(0,len(cob_call_freqs),3):
    x = i
    cob_call_freq.append(cob_call_freqs[x:x+3])
    wth_call_freq.append(wth_call_freqs[x:x+3])
    pee_call_freq.append(pee_call_freqs[x:x+3])
    crn_call_freq.append(crn_call_freqs[x:x+3])
for i in range(0,len(sos_call_freqs),6):
    x = i
    sos_call_freq.append(sos_call_freqs[x:x+6])

cob_call_freq_final = []
wth_call_freq_final = []
sos_call_freq_final = []
pee_call_freq_final = []
crn_call_freq_final = []  

for i in cob_call_freq:
    cob_call_freq_final.append(sum(i))
for i in wth_call_freq:
    wth_call_freq_final.append(sum(i))
for i in sos_call_freq:
    sos_call_freq_final.append(sum(i))
for i in pee_call_freq:
    pee_call_freq_final.append(sum(i))
for i in crn_call_freq:
    crn_call_freq_final.append(sum(i))
  
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


cob_variants = final_variant_list[0] + final_variant_list[4] + final_variant_list[5]
wth_variants = final_variant_list[1]+ final_variant_list[8]+final_variant_list[9]
sos_variants = final_variant_list[2] + final_variant_list[3] + final_variant_list[12] + final_variant_list[14] + final_variant_list[15] + final_variant_list[16]
pee_variants = final_variant_list[6] + final_variant_list[7] + final_variant_list[13]
crn_variants = final_variant_list[10] + final_variant_list[11] + final_variant_list[17]

cob_variant = set(cob_variants)
wth_variant = set(wth_variants)
sos_variant = set(sos_variants)
pee_variant = set(pee_variants)
crn_variant = set(crn_variants)

variant_list = [cob_variant,wth_variant,sos_variant,pee_variant,crn_variant]
iterator_base = 0
iterator_adj = 1
pair_order_list2 = []
pair_order_list3 = []
while iterator_base < (len(variant_list)):
         iterator_adj = 1 + iterator_base
         while iterator_adj < len(variant_list):
             x = [variant_list[(iterator_base)],variant_list[(iterator_adj)]]
             pair_order_list2.append(x)
             
             iterator_adj += 1
         iterator_base += 1
    
for i in pair_order_list2:
    a = set(i[0])
    b = set(i[1])
    c = list((set(a) | set(b) - set(a) & set(b)))
    pair_order_list3.append(c)

cob = pd.read_csv("cob.csv")
wth = pd.read_csv("wth.csv")
sos = pd.read_csv("sos.csv")
pee = pd.read_csv("pee.csv")
crn = pd.read_csv("crn.csv")

dataframe_list = [cob,wth,sos,pee,crn]
# =============================================================================

pop_1_pos = []
pop_1_count = []
pop_1_allele_freq = []
# 
pop_2_pos = []
pop_2_count = []
pop_2_allele_freq = []

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
                    pop_1_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
                    
                iterator_test += 1
                for i in dataframe_test[iterator_test]["pos"]:
                    if i in (pair_order_list3[iterator_pairwise]):
                        pop_2_pos.append(i)
                for j in pop_2_pos:
                    pop_2_count.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'count'].iloc[0])
                    pop_2_allele_freq.append(dataframe_test[iterator_test].loc[dataframe_test[iterator_test]['pos']==j,'allele_freq'].iloc[0])
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
                pop_total_allele_sum_list.append((pop_1_allele_freq[iterator_calc]*pop_1_count[iterator_calc])+(pop_2_allele_freq[iterator_calc]*pop_2_count[iterator_calc]))
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
            print(ht)
            print(hs)
            print(iterator_base)
            print(iterator_adj)
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

print(len(fst_list))
print((sum(fst_list))/(len(fst_list)))
df = pd.DataFrame(fst_list, columns=["colummn"])
df.to_csv('population_fst_list.csv', index=False)