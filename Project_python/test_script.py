#cmd+S = save
#ctrl+` = opens terminal
#Shortcut to run script =ctrl+alt+n

#Any functions are called using ()
#All text in '' or ""
#Triple quoates can be used to format long strings of text '''

#Print fucntion is used to print messages
print('hello')

#Simple maths can be done in python
print(2+2)
print(1>3)

#Able to set variables
counts = 100
print(counts)

#Python uses 0 index thereofre the first character/number in an object is labelled 0 EXAMPLE below
matrix_a = [2, 4, 6, 8, 10]
#Length
print(len(matrix_a))
#First and last character respectively
print(matrix_a[0])
print(matrix_a[-1])

#\ in python strings allow python to ignore the following character if it would affect the function
#\n causes new line
print('Python \'Programming')

#joining expressions
first = 'Tara'
last = 'Morrison'
full = f'{first} {last}'

#Integrating user input
x = input('x:')
print(type(x))

y = int(x) + 1 #int(x) converts x into an integer
#float(x) - converts to floating number (real value will include decimal place)
#str(x) - converts to string (text)
#bool(x) - returns true or false

#
import numpy as np
import pandas as pd

with open('test.txt', 'r+') as matrix:
    array = matrix.read().split()
    print(array)

import pandas as pd
with open('test.txt', 'r+') as matrix:
    data_table=pd.read_csv(matrix, index_col=0)
print(data_table.iloc[:5, :5])




####################################################









