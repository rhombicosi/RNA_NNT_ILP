from fnmatch import fnmatchcase
import sys
import math
import copy
            
from score_matrix import bp_df, stack_df 

with open('1etg', 'r') as file:
    MOTIF = file.read().rstrip()    
    print (MOTIF)

length = len(MOTIF)
print (length)

listP = ""      # listP will sum the edge-to-edge pairs variables
listQ = ""      # listQ will sum the face-to-face stacking variables
listS = ""

# OUT = open('2oe6-motifs-dl.lp', "w")
OUT = open('1etg-motifs-dl.lp', "w")
# OUT = open('1lnt-motifs-dl.lp', "w")
# OUT = open('2lx1-motifs-dl.lp', "w")
# OUT = open('3dvz-motifs-dl.lp', "w")
# OUT = open('5jf4-motifs-dl.lp', "w")
OUT.write ("Maximize \n")
OUT.write ("E \n")

length_for_range = length + 1
weights = [[0 for x in range(1, length_for_range + 1)] for y in range(1, length_for_range + 1)]

# strands
a_start = 1
a_end = 4

b_start = 5
b_end = 9

A = range(a_start,a_end)
B = range(b_start,b_end)

# stacks floors
STRAND_A = MOTIF[a_start-1:a_end]
STRAND_A = [x for x in STRAND_A]
STRAND_A.reverse()
STRAND_B = MOTIF[b_start-1:b_end]
STRAND_B = [x for x in STRAND_B]
# STRAND_B.reverse()
print(STRAND_A)
print(len(STRAND_A))

print(STRAND_B)
print(len(STRAND_B))

FLOOR_1 = {}
FLOOR_2 = {}
FLOOR_3 = {}
FLOOR_4 = {}

m = len(STRAND_A)-2
n = len(STRAND_B)-2
print(m)
print(n)

# total number of floors is less than or equal to 
# max possible number of floors
fs = min(m,n)

# build 1st floor
f1=1

FLOOR_1[(0,0)] = []
for i in range(f1,m-fs+f1+1):
    for j in range(f1,n-fs+f1+1):
        FLOOR_1[(0,0)].append((i,j))

# all other floors
def next_floors(floor_prev, m, n, fn, fs):
    floor_next={}
    for f in floor_prev:
        floor_next[f] = []
        for i in range(f[0]+1,m-fs+fn+1):
            for j in range(f[1]+1,n-fs+fn+1):
                floor_next[f].append((i,j))
    return floor_next


FLOOR_2 = next_floors(FLOOR_1[(0,0)], m, n, 2, fs)
# FLOOR_3 = next_floors(FLOOR_2[(1,1)], m, n, 3, fs)
# FLOOR_4 = next_floors(FLOOR_3[(2,2)], m, n, 4, fs)

for k,v in FLOOR_1.items():
    print(f'{k}:{v}')

for k,v in FLOOR_2.items():
    print(f'{k}:{v}')

# for k,v in FLOOR_3.items():
#     print(f'{k}:{v}')

# for k,v in FLOOR_4.items():
#     print(f'{k}:{v}')

# build first possible floors
poss_floors_1 = []
for f in FLOOR_1[(0,0)]:
    poss_floors_1.append([f])

for pf in poss_floors_1:
    print(pf)

def all_poss_floors(poss_floors, floor_next):

    all_floors = []
    fff = []
    for ff in poss_floors:
        fff = copy.deepcopy(ff)
        for i in floor_next[ff[-1]]:
            fff.append(i)
            all_floors.append(copy.deepcopy(fff))
            fff.pop()

    return all_floors

poss_floors_2 = all_poss_floors(poss_floors_1, FLOOR_2)
# poss_floors_3 = all_poss_floors(poss_floors_2, FLOOR_3)
# poss_floors_4 = all_poss_floors(poss_floors_3, FLOOR_4)

for pf in poss_floors_2:
    print(pf)