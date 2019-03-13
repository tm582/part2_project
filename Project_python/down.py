#Import relevant functions that are used in the script
import numpy as np
import pandas as pd
import random
from statistics import mean

#Import csv matrix file as data frame using pandas
x = pd.read_csv('gtex_counts.csv', index_col=0)
orix=pd.read_csv('gtex_counts.csv', index_col=0)


#Calculate column median & total average of the medians
means=x.mean()
print("The median of each sample is\n",means)

avemean=mean(means)


#Check number of columns and rows
ncol = (len(x.columns))


nrow = (len(x.index))

change = 0

while change < 33:
    col = random.randint(0,(ncol-1))
    row = random.randint(0,(nrow-1))

    oricount = x.iloc[row,col]

    if oricount>0:
        newcount=oricount-1
    else:
        newcount=oricount
    
    print(x.iat[row,col])
    x.iat[row,col]=newcount
    print(x.iat[row,col])
    means=x.mean()

    newavemeans=mean(means)
    #print('mean')
    print(newavemeans)
    change=((avemean-newavemeans)/avemean)*100
    #print('change')
    #print(change)
    if change==33:
        break


x.to_csv('gtex_down33.csv', sep='\t')
