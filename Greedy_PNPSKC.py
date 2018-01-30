from __future__ import division
import pandas as pd
from math import sqrt,log
from sklearn.metrics import auc,precision_recall_curve
from sklearn.metrics import roc_curve
import time
import numpy as np
from sklearn.metrics import accuracy_score
import sys

# Author: Yichao Li
# The source code is an implementation of the Greedy_PNPSKC algorithm in "A Novel Discriminative Set Multi-Cover Model for Discovering DNA Motifs and Motif Pairs".




def greedy_WSKC(df,k,W):
	my_selected_sets = []
	uncovered_elements = df.index.tolist()
	unused_sets = df.columns.tolist()
	cost = 0
	covered_elements = []
	cover_index = {}
	count = 0
	for x in uncovered_elements:
		cover_index[x] = k
	while True:
		if len(unused_sets) == 0:
			return my_selected_sets,uncovered_elements,cost		
		count += 1
		current_df = df.ix[uncovered_elements]
		min_price = 9999999
		min_price_set = ""
		time1 = time.time()
		for s in unused_sets:			
			CE_s = len(set(covered_elements_dict[s]).difference(covered_elements))
			if CE_s == 0:
				continue		
			W_s = W[s]
			P_s = W_s/CE_s
			if P_s < min_price:
				min_price = P_s
				min_price_set = s
		if min_price == 9999999:
			return my_selected_sets,uncovered_elements,cost
		my_selected_sets.append(min_price_set)
		unused_sets.remove(min_price_set)
		covered_elements = check_kCover(df[my_selected_sets],k)
		new_uncovered_elements = list(set(uncovered_elements).difference(covered_elements))
		for x in covered_elements_dict[min_price_set]:
			if cover_index[x] > 0:
				cover_index[x] -= 1
		uncovered_elements = new_uncovered_elements
		cost+= W[min_price_set]
		
def check_kCover(df,k):
	return df.index[df.sum(axis=1) >= k].tolist()


def greedy_RBSKC(df,k):
	# df has a class columns , df.columns[-1]
	# This algorithm is a modification for Greedy_RB
	# subsetting df
	
	df_for_wskc = df[df[df.columns[-1]] == 1][df.columns[:-1]]
	weight = {}
	# creating weights
	df_for_weights = df[df[df.columns[-1]] != 1][df.columns[:-1]]
	for s in df_for_weights.columns.tolist():
		weight[s] = df_for_weights[s].sum()
	my_selected_sets,uncovered_elements,cost = greedy_WSKC(df_for_wskc,k,weight)
	return my_selected_sets,uncovered_elements,cost
	

def LOW_DEG_RBSKC(df,k,X):
	# discard from S the sets with more than X red elements
	kept_sets = []
	
	for s in negative_data.columns.tolist():
		if negative_data[s].sum() <= X:
			kept_sets.append(s)
	# if the new df can't cover all elements, return all sets as a solution
	if len(check_kCover(positive_df[kept_sets],k)) < positive_size:
		print X,"didn't go through"
		return df.columns.tolist()[:-1],True
	df_kept_sets = df[kept_sets+[df.columns[-1]]]
	# find high degree red elements = row sum
	
	# drop high red
	df_kept_sets_drop_high_red = df_kept_sets.drop(high_red)
	# apply greedy_RBSKC
	my_selected_sets,uncovered_elements,cost = greedy_RBSKC(df_kept_sets_drop_high_red,k)
	return my_selected_sets,False
	


def LOW_DEG_RBSKC2(df,k):
	# for X from min covered number to max covered number
	temp = negative_data.sum().sort_values().tolist()
	min_covered_num_neg = int(temp[0])
	if min_covered_num_neg == 0:
		min_covered_num_neg = 1
	max_covered_num_neg = int(temp[-1])
	print max_covered_num_neg,min_covered_num_neg
	negative_size = negative_data.shape[0]
	my_solution = ""
	min_cost = 99999999999999
	num_iter_not_improved = 10
	counter = 0
	count_didnt_go_through_number = 0
	X = min_covered_num_neg
	while X <= max_covered_num_neg:
		print "negative iteration:",X
		selected_sets,flag = LOW_DEG_RBSKC(df,k,X)
		cost = calculate_RBSKC_cost(df,k,selected_sets)
		if cost < min_cost:
			counter = 0
			my_solution = selected_sets
			min_cost = cost
		else:
			if cost == min_cost:
				if len(selected_sets) < len(my_solution):
					my_solution = selected_sets
			counter += 1

		print len(my_solution),min_cost
		X += 1
	return my_solution


def clean_data(df,k):
	df_for_clean = df[df[df.columns[-1]] == 1][df.columns[:-1]]
	uncovered_k_elements = df_for_clean[df_for_clean.sum(axis=1)<k].index.tolist()
	print "cleaning for: ",len(uncovered_k_elements)," elements by: ",k
	return df.drop(uncovered_k_elements)

def summerize_result(df,k,selected_sets):
	# TP,FP,TN,FN,cost,coverage_for_each_element,ACC,auPRC,auROC
	# df = df[selected_sets+[df.columns[-1]]]
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
	# cost
	# df = df[selected_sets+[df.columns[-1]]]
	y_true = df[df.columns[-1]].tolist()
	# print "y_true",y_true
	y_pred = (df[selected_sets].sum(axis=1) >=k ).astype(int).tolist()
	# print "y_pred",y_pred
	# y_score = df[df.columns[:-1]].sum(axis=1).tolist()
	# print "y_score",y_score
	# se,sp,acc,tp,tn,cost = calculate_SE_SP_ACC(y_pred,y_true)
	# auROC,auPRC = calculate_auROC_auPRC(y_true,y_score)
	return len(y_true) - accuracy_score(y_true, y_pred, normalize=False)
	
def calculate_RBSKC_cost(df,k,selected_sets):
	# cost
	# df = df[selected_sets+[df.columns[-1]]]
	y_true = df[df.columns[-1]].tolist()
	# print "y_true",y_true
	y_pred = (df[selected_sets].sum(axis=1) >=k ).astype(int).tolist()
	count = 0
	for i in range(len(y_true)):
		y_i = y_true[i]
		if y_i == 1:
			continue
		y_p = y_pred[i]
		if y_p == 1:
			count+= 1
		
	# print "y_pred",y_pred
	# y_score = df[df.columns[:-1]].sum(axis=1).tolist()
	# print "y_score",y_score
	# se,sp,acc,tp,tn,cost = calculate_SE_SP_ACC(y_pred,y_true)
	# auROC,auPRC = calculate_auROC_auPRC(y_true,y_score)
	return count
		
	
# functions for calculating machine learning metrics
def calculate_SE_SP_ACC(y_pred,y_true,wrt=1):
	# diabetes are 1
	#print y_pred,y_true
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
	# print "TP:",tp
	# print "TN:",tn
	T = y_true.count(wrt)
	N = total - T
	#print T,N,total
	if T * N == 0:
		print "T or N is zero!"
		return "error"
	cost = total - correct
	#print total
	return tp/T,tn/N,correct/total,tp,tn,cost

def calculate_auROC_auPRC(y,scores):
	fpr, tpr, _ = roc_curve(y, scores, pos_label=1)
	auROC = auc(fpr,tpr)
	precision, recall, _ = precision_recall_curve(y,scores)
	auPRC = auc(recall, precision)
	return auROC,auPRC
from copy import deepcopy as dp
def backward_search(df,k,selected_sets,old_cost):
	# eliminate a set if it can reduce the cost
	removed_set = ""
	for s in selected_sets:
		temp = dp(selected_sets)
		# print temp
		temp.remove(s)
		# print "remvoed:",s," ",temp
		if len(temp) == 0:
			continue
		cost = calculate_cost(df,k,temp)
		print "new cost ",cost
		if cost <= old_cost or cost == 0:
			# print "go in for:", s 
			old_cost = cost
			removed_set = s
			
	return removed_set,old_cost
# forward search is not meaningful	
def forward_search(df,k,selected_sets,old_cost):
	# eliminate a set if it can reduce the cost
	all_sets = df.columns[:-1].tolist()
	added_set = ""
	unused_set = list(set(all_sets).difference(selected_sets))
	for s in unused_set:
		cost = calculate_cost(df,k,selected_sets + [s])
		if cost < old_cost:
			old_cost = cost
			added_set = s
	return added_set,old_cost

def backward_elimination(df,k,selected_sets):
	old_cost = calculate_cost(df,k,selected_sets)
	print "backward_elimination old cost ",old_cost
	new_cost = 0
	# old_set = selected_sets
	while new_cost <= old_cost or new_cost == 0:
		removed_set,new_cost = backward_search(df,k,selected_sets,old_cost)
		if removed_set == "":
			print "nothing to eliminate"
			return selected_sets
		# print "backward_elimination:",removed_set
		old_set = dp(selected_sets)
		selected_sets.remove(removed_set)
		# se,sp,acc,tp,tn,cost,auROC,auPRC = summerize_result(df,k,selected_sets)
		# if se < 0.9:
			# return old_set
		
	return selected_sets

def forward_addition(df,k,selected_sets):
	old_cost = calculate_cost(df,k,selected_sets)
	new_cost = 0
	while new_cost < old_cost:
		added_set,new_cost = forward_search(df,k,selected_sets,old_cost)
		if added_set == "":
			return selected_sets
		selected_sets.append(added_set)
	return selected_sets

	





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
print "selected_sets_back:",len(selected_sets_back)
print outresult
print selected_sets_back
print >>greedy_result,len(selected_sets_back)
print >>greedy_result,outresult[0]
print >>greedy_result,1-outresult[1]
print >>greedy_result,outresult[3]
print >>greedy_result,outresult[4]
print >>greedy_result,"\t".join(selected_sets_back)
print >>greedy_result,"\t".join(map(lambda x:str(x),outresult))


