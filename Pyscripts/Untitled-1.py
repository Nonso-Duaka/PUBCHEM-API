def get(items, x):
    return items[x]

def counts(r, t):
    count= 0
    for i in r:
        if i>t:
            count+=1
    print(count) 




def can_pair(item_quantities):
	for i in range(len(item_quantities)):
          if item_quantities[i] % 2 != 0:
              print (False)
          else:
              print (True)



item_quantities = [2, 4, 6, 8]
can_pair(item_quantities)
