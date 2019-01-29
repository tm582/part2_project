#Import relevant functions that are used in the script
import numpy as np
import pandas as pd
import random
from statistics import mean

#Import csv matrix file as data frame using pandas
x = pd.read_csv('test.csv', index_col=0)
orix=pd.read_csv('test.csv', index_col=0)


#Calculate column median & total average of the medians
medians=x.median()
print("The median of each sample is\n",medians)

avemed=mean(medians)
print("The average of medians is\n",avemed)


#Check number of columns and rows
ncol = (len(x.columns))
print("Total #columns in the data frame = ", ncol)

nrow = (len(x.index))
print("Total #rows in the data frame = ", nrow)

change = 0


#Begin while loop where each iteration one count is sleected and reduced by one. Repeat until the mean of the medians of each sample has be reduced by 30%
while change < 30:

    #Generate random column number and random row number
    col = random.randint(0,(ncol-1))
    row = random.randint(0,(nrow-1))

    print("Coordinate of count to be changed is ROW = ", row,"and COLUMN = ",col)


    #Use coordinates to find random count
    oricount = x.iloc[row,col]
    print("The original count is", oricount)

    #Reduce count by 1
    if oricount>0:
        print("Original count is greater than 0. Can proceed to -1")
        newcount=oricount-1
        print(newcount)
    else:
        print("Original count = 0. Therefore do not adjust and repeat on new random count.")
        newcount=oricount

    #Check new count
    print(newcount)

    #Replace old count for new count
    x.iat[row,col]=newcount

    print(x)

    #Calculate NEW column medians and total average of medians
    medians=x.median()
    print("The NEW median of each sample is\n",medians)

    newavemed=mean(medians)
    print("The NEW average of medians is\n",avemed)

    change=((avemed-newavemed)/avemed)*100

    print("Percentage of original =", change,"%")
    
    if change==30:
        break


print("The ORIGINAL matrix of counts\n", orix)
print("The NEW matrix of counts - has been reduced by 30%\n", x)
