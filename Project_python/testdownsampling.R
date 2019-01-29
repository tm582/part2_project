test=read.table('test.txt', header = T, row.names = 1)

randompoint=(test[sample(nrow(test),1),sample(ncol(test),1)])-1
randompoint
newpoint=randompoint-1



randomcol=sample(ncol(test),1)
randomcol
randomrow=sample(nrow(test),1)
randomrow

randompoint=test[randomrow,randomcol]

if (randompoint>0){
  print('Original Count > 0. Therefore minus 1 from count and replace')
  newpoint=test[randomrow, randomcol]-1
} else{
  print('Original Count = 0. Will therefore remain unchanged')
  newpoint=randompoint
}

randompoint
newpoint

test[randomrow,randomcol]=newpoint