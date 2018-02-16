from __future__ import division
import pandas as pd
from math import sqrt,log
from sklearn.metrics import auc,precision_recall_curve
from sklearn.metrics import roc_curve
import numpy as np
from sklearn.metrics import accuracy_score
import sys
from copy import deepcopy as dp

# Author: Yichao Li
# Problem: Positive Negative Partial Set K-Cover 
# The source code is an implementation of the Greedy_PNPSKC algorithm in "A Novel Discriminative Set Multi-Cover Model for Discovering DNA Motifs and Motif Pairs".

# Algotithm 1
def greedy_WSKC(df,k,W):
	# 1. C = empty
	my_selected_sets = []
	# 2. t_i = k
	cover_index = {}
	uncovered_elements = df.index.tolist()
	for x in uncovered_elements:
		cover_index[x] = k
	# for the for loop below
	unused_sets = df.columns.tolist()
	# init
	cost = 0
	# for calculating P_j
	covered_elements = []		
	while True:	
		min_price = 9999999
		min_price_set = ""
		# 4. for loop to find s_j
		for s in unused_sets:
			# 3. P_j is shown below
			# covered_elements_dict[s] is pre-computed for speed issue
			CE_s = len(set(covered_elements_dict[s]).difference(covered_elements))
			if CE_s == 0:
				continue		
			W_s = W[s]
			P_s = W_s/CE_s
			if P_s < min_price:
				min_price = P_s
				min_price_set = s
		# NOTE: if use t_i for termination conditions, when the elements cannot be covered even if all the sets are used, this can create an infinite loop. So we need to use the change of min_price to terminate the while loop.
		if min_price == 9999999:
			return my_selected_sets,uncovered_elements,cost
		# 5. add the set
		my_selected_sets.append(min_price_set)
		unused_sets.remove(min_price_set)
		covered_elements = check_kCover(df[my_selected_sets],k)
		new_uncovered_elements = list(set(uncovered_elements).difference(covered_elements))
		# 6. update cover index
		for x in covered_elements_dict[min_price_set]:
			if cover_index[x] > 0:
				cover_index[x] -= 1
		uncovered_elements = new_uncovered_elements
		cost+= W[min_price_set]

# Algotithm 2
def greedy_RBSKC(df,k):
	# convert to WSKC instance
	df_for_wskc = df[df[df.columns[-1]] == 1][df.columns[:-1]]
	weight = {}
	# creating weights
	df_for_weights = df[df[df.columns[-1]] != 1][df.columns[:-1]]
	for s in df_for_weights.columns.tolist():
		weight[s] = df_for_weights[s].sum()
	my_selected_sets,uncovered_elements,cost = greedy_WSKC(df_for_wskc,k,weight)
	return my_selected_sets,uncovered_elements,cost
	
# Algotithm 3
def LOW_DEG_RBSKC(df,k,X):
	# 1. Discard from S the sets with more than X red elements
	kept_sets = []
	
	for s in negative_data.columns.tolist():
		if negative_data[s].sum() <= X:
			kept_sets.append(s)
	# 2. if the new df can't cover all elements, return all sets as a solution
	if len(check_kCover(positive_df[kept_sets],k)) < positive_size:
		print X,"didn't go through"
		return df.columns.tolist()[:-1]
	df_kept_sets = df[kept_sets+[df.columns[-1]]]
	# 3 & 4. high_red element is pre-computed for speed issue
	
	# 5. drop high red
	df_kept_sets_drop_high_red = df_kept_sets.drop(high_red)
	# 6. apply greedy_RBSKC
	my_selected_sets,uncovered_elements,cost = greedy_RBSKC(df_kept_sets_drop_high_red,k)
	return my_selected_sets
	
# Algotithm 4
def LOW_DEG_RBSKC2(df,k):
	# initialize the range for X
	# for X from min covered number to max covered number
	temp = negative_data.sum().sort_values().tolist()
	min_covered_num_neg = int(temp[0])
	if min_covered_num_neg == 0:
		min_covered_num_neg = 1
	max_covered_num_neg = int(temp[-1])
	print max_covered_num_neg,min_covered_num_neg
	my_solution = ""
	min_cost = 99999999999999
	X = min_covered_num_neg
	# find the solution with the lowest cost
	# if two costs are the same, use the one with smaller number of sets
	while X <= max_covered_num_neg:
		print "X =",X
		selected_sets = LOW_DEG_RBSKC(df,k,X)
		cost = calculate_RBSKC_cost(df,k,selected_sets)
		if cost < min_cost:
			my_solution = selected_sets
			min_cost = cost
		else:
			if cost == min_cost:
				if len(selected_sets) < len(my_solution):
					my_solution = selected_sets

		print "Number of selected motifs:",len(my_solution),"Cost:",min_cost
		X += 1
	return my_solution

# Algorithm 5
def backward_search(df,k,selected_sets,old_cost):
	# eliminate a set if it can reduce the cost
	# 1. the set to be removed
	removed_set = ""
	for s in selected_sets:
		# 2. Q_j tilte
		temp = dp(selected_sets)
		temp.remove(s)
		if len(temp) == 0:
			continue
		# 3. new cost
		cost = calculate_cost(df,k,temp)
		# print "new cost ",cost
		# 4 & 5 update old cost and set
		if cost <= old_cost or cost == 0:
			old_cost = cost
			removed_set = s
			
	return removed_set,old_cost

# Algorithm 6
def backward_elimination(df,k,selected_sets):
	old_cost = calculate_cost(df,k,selected_sets)
	print "backward_elimination old cost ",old_cost
	new_cost = 0
	while True:
		# 1. run backward search
		removed_set,new_cost = backward_search(df,k,selected_sets,old_cost)
		# 2. if nothing to remove
		if removed_set == "":
			print "nothing to eliminate"
			return selected_sets
		# 3. update
		selected_sets.remove(removed_set)
		# 4. update
		old_cost = new_cost
		
	return selected_sets

def check_kCover(df,k):
	return df.index[df.sum(axis=1) >= k].tolist()	
	
# clean the dataset, make sure if all the sets are used, all the elements can be k-covered
def clean_data(df,k):
	df_for_clean = df[df[df.columns[-1]] == 1][df.columns[:-1]]
	uncovered_k_elements = df_for_clean[df_for_clean.sum(axis=1)<k].index.tolist()
	print "cleaning for: ",len(uncovered_k_elements)," elements by: ",k
	return df.drop(uncovered_k_elements)

# below are some evaluation functions
# even though auROC and auPRC are included, they are not used for the PNPSKC	
def summerize_result(df,k,selected_sets):
	# TP,FP,TN,FN,cost,coverage_for_each_element,ACC,auPRC,auROC
	y_true = df[df.columns[-1]].astype(int).tolist()
	# print "y_true",y_true
	y_pred = (df[selected_sets].sum(axis=1) >=k ).astype(int).tolist()
	# print "y_pred",y_pred
	y_score = df[df.columns[:-1]].sum(axis=1).tolist()
	# print "y_score",y_score
	se,sp,acc,tp,tn,cost = calculate_SE_SP_ACC(y_pred,y_true)
	auROC,auPRC = calculate_auROC_auPRC(y_true,y_score)
	return se,sp,acc,tp,tn,cost,auROC,auPRC
	
def calculate_cost(df,k,selected_sets):
	y_true = df[df.columns[-1]].tolist()
	y_pred = (df[selected_sets].sum(axis=1) >=k ).astype(int).tolist()
	return len(y_true) - accuracy_score(y_true, y_pred, normalize=False)
	
def calculate_RBSKC_cost(df,k,selected_sets):
	y_true = df[df.columns[-1]].tolist()
	y_pred = (df[selected_sets].sum(axis=1) >=k ).astype(int).tolist()
	count = 0
	for i in range(len(y_true)):
		y_i = y_true[i]
		if y_i == 1:
			continue
		y_p = y_pred[i]
		if y_p == 1:
			count+= 1
	return count
		
# functions for calculating machine learning metrics
def calculate_SE_SP_ACC(y_pred,y_true,wrt=1):
	tp = 0.0
	tn = 0.0
	correct = 0.0
	if len(y_pred) != len(y_true):
		print "len(y_pred) != len(y_true)!"
		exit()
	for i in range(len(y_pred)):
		if y_pred[i] == y_true[i]:
			if y_pred[i] == wrt:
				tp += 1
			else:
				tn += 1
			correct += 1
	total = len(y_true)
	T = y_true.count(wrt)
	N = total - T
	if T * N == 0:
		print "T or N is zero!"
		return "error"
	cost = total - correct
	return tp/T,tn/N,correct/total,tp,tn,cost

def calculate_auROC_auPRC(y,scores):
	fpr, tpr, _ = roc_curve(y, scores, pos_label=1)
	auROC = auc(fpr,tpr)
	precision, recall, _ = precision_recall_curve(y,scores)
	auPRC = auc(recall, precision)
	return auROC,auPRC

data = pd.read_csv(sys.argv[1],index_col=0)
k=int(sys.argv[2])
greedy_result = open("Greedy."+str(k)+".result","wb")

# preprocessing
df = clean_data(data,k)
Total_sets = data.columns[:-1]
negative_data = data[data[data.columns[-1]] != 1][Total_sets]
positive_df = df[df[df.columns[-1]] == 1]
positive_size = positive_df.shape[0]
background_size = negative_data.shape[0]
Y = sqrt((df.shape[1]-1)/log(positive_size))
high_red = negative_data[negative_data.sum(axis=1)>Y].index.tolist()
covered_elements_dict = {}
for c in data.columns[:-1]:
	covered_elements_dict[c] = positive_df[positive_df[c] >=1].index.tolist()

# Greedy_PNPSKC
selected_sets = LOW_DEG_RBSKC2(df,k)
selected_sets_back = backward_elimination(data,k,selected_sets)

# output result
outresult = summerize_result(data,k,selected_sets_back)
print "Length of the selected set:",len(selected_sets_back)
print "Fraction of the positives being k-covered:",outresult[0]
print "Fraction of the negatives being k-covered:",1-outresult[1]
print "Number of positives being k-covered:",outresult[3]
print "Number of negatives that are covered less than k times:",outresult[4]
print "The selected motifs are shown below:"
print "\n".join(selected_sets_back)
print >>greedy_result,"Length of the selected set:",len(selected_sets_back)
print >>greedy_result,"Fraction of the positives being k-covered:",outresult[0]
print >>greedy_result,"Fraction of the negatives being k-covered:",1-outresult[1]
print >>greedy_result,"Number of positives being k-covered:",outresult[3]
print >>greedy_result,"Number of negatives that are covered less than k times:",outresult[4]
print >>greedy_result,"The selected motifs are shown below:"
print >>greedy_result,"\t".join(selected_sets_back)


