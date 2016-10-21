import sys

first = True
for line in sys.stdin:
    elems = line.split('\t')
    if (elems[0] == 'U'):
        if first:
            first = False
        else:
            print("")
        print(elems[1], end=" ")
    else:
        print(elems[2], end=" ")
