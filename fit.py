iimport scipy
import numpy
import numpy as np

l=3
n=l*l
arr=[]

po = open("permutations.txt", "wb")
test = open("test.txt", "wb")


for i in range(2 ** n):
	#Generate permutations in a long string
    a = format(i, '0{}b'.format(n))
    k=0
    #split up string adn write to file
    for j in range(0,l):
       po.write(a[k:k+l]+"\n")
       for m in range(0, l):
       		arr.append(a[k+m])
       k=k+l

    po.write("#\n")
    


arr2=[]
count = 0
#Add rows to 2d array to form a plane type structure
for i in xrange(0, len(arr), l):
	arr2.append([])
	for j in range(0, l):
		arr2[count].append(arr[i+j])
	count+=1

print arr2
po.close()
