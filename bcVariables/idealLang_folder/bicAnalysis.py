import time
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import numpy as np
import math
import os
import sys

data_set = 'idealEley_fitting_folder/'
print(data_set)
k = 4

#data_set = 'idealLang_fitting_folder/'
#k = 6

#data_set = 'idealLang_fittingBoth_folder/'
#k = 8

n = 6000

df1CO = pd.read_csv('./flux_data/CO.csv',header=None)
df2CO = pd.read_csv('./../'+data_set+'/flux_data/CO.csv',header=None)
df1O2 = pd.read_csv('./flux_data/O2.csv',header=None)
df2O2 = pd.read_csv('./../'+data_set+'/flux_data/O2.csv',header=None)
df1CO2 = pd.read_csv('./flux_data/CO2.csv',header=None)
df2CO2 = pd.read_csv('./../'+data_set+'/flux_data/CO2.csv',header=None)

df3 =   (df1CO2.iloc[:,1] - df2CO2.iloc[:,1])**2 + (df1CO.iloc[:,1] - df2CO.iloc[:,1])**2 + (df1O2.iloc[:,1] - df2O2.iloc[:,1])**2 

test = df3.sum()/n


print(n*math.log(test/n) + 2*k)
print(k*math.log(n) + n*math.log(test))
print(test)


#data_set = 'idealEley_fitting_folder/'
#k = 4

data_set = 'idealLang_fitting_folder/'
print(data_set)
k = 6

#data_set = 'idealLang_fittingBoth_folder/'
#k = 8

n = 6000

df1CO = pd.read_csv('./flux_data/CO.csv',header=None)
df2CO = pd.read_csv('./../'+data_set+'/flux_data/CO.csv',header=None)
df1O2 = pd.read_csv('./flux_data/O2.csv',header=None)
df2O2 = pd.read_csv('./../'+data_set+'/flux_data/O2.csv',header=None)
df1CO2 = pd.read_csv('./flux_data/CO2.csv',header=None)
df2CO2 = pd.read_csv('./../'+data_set+'/flux_data/CO2.csv',header=None)

df3 =   (df1CO2.iloc[:,1] - df2CO2.iloc[:,1])**2 + (df1CO.iloc[:,1] - df2CO.iloc[:,1])**2 + (df1O2.iloc[:,1] - df2O2.iloc[:,1])**2 

test = df3.sum()/n


print(n*math.log(test/n) + 2*k)
print(k*math.log(n) + n*math.log(test))
print(test)



#data_set = 'idealEley_fitting_folder/'
#k = 4

#data_set = 'idealLang_fitting_folder/'
#k = 6

data_set = 'idealLang_fittingBoth_folder/'
print(data_set)
k = 8

n = 6000

df1CO = pd.read_csv('./flux_data/CO.csv',header=None)
df2CO = pd.read_csv('./../'+data_set+'/flux_data/CO.csv',header=None)
df1O2 = pd.read_csv('./flux_data/O2.csv',header=None)
df2O2 = pd.read_csv('./../'+data_set+'/flux_data/O2.csv',header=None)
df1CO2 = pd.read_csv('./flux_data/CO2.csv',header=None)
df2CO2 = pd.read_csv('./../'+data_set+'/flux_data/CO2.csv',header=None)

df3 =   (df1CO2.iloc[:,1] - df2CO2.iloc[:,1])**2 + (df1CO.iloc[:,1] - df2CO.iloc[:,1])**2 + (df1O2.iloc[:,1] - df2O2.iloc[:,1])**2 

test = df3.sum()/n


print(n*math.log(test/n) + 2*k)
print(k*math.log(n) + n*math.log(test))
print(test)



#total_value = 1
#
#for k in range(0,len(df1CO2.iloc[:,1].tolist())):
#	val = (1/math.sqrt(2*3.14159*std)) * math.exp(-((df1CO2.iloc[k,1] - df2CO2.iloc[k,1])**2)/(2*std))
#	#print(val)
#	total_value *= val
#
#print(mean)
#print(std)
#print(total_value)
#sys.exit()
#
#for k in range(0,len(df1CO2.iloc[:,1].tolist())):
#	val = (1/math.sqrt(2*3.14159*std))*math.exp(-((df1CO.iloc[k,1] - df2CO.iloc[k,1])**2)/(2*std))
#	total_value *= val
#
#for k in range(0,len(df1CO2.iloc[:,1].tolist())):
#	val = (1/math.sqrt(2*3.14159*std))*math.exp(-((df1O2.iloc[k,1] - df2O2.iloc[k,1])**2)/(2*std))
#	total_value *= val
#
#print(total_value)

plt.show()