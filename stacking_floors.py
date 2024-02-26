#%%
from fnmatch import fnmatchcase
import copy
            
from score_matrix import bp_df, stack_df 

def read_sqnc(sqncname):
    with open(sqncname, 'r') as file:
        MOTIF = file.read().rstrip()    
        # print (MOTIF)

    # length = len(MOTIF)
    # print (length)

    return MOTIF

def internal_stacks_sizes(MOTIF, a_start, a_end, b_start, b_end):

    # cut sequence into two strands
    STRAND_A = MOTIF[a_start-1:a_end]
    STRAND_A = [x for x in STRAND_A]
    STRAND_A.reverse()
    STRAND_B = MOTIF[b_start-1:b_end]
    STRAND_B = [x for x in STRAND_B]

    # print(STRAND_A)
    # print(len(STRAND_A))

    # print(STRAND_B)
    # print(len(STRAND_B))

    m = len(STRAND_A)
    n = len(STRAND_B)
    # print(m)
    # print(n)

    return m, n

# max stacking size (total number of floors) 
# fs = min(m,n)
# fn - floor number
# fs - floors total size
def build_possible_floors(m, n, floor, block):   
    if block > min(m,n):
        print("Error: total number of floors should be less than number of nts in the smallest strand")
    elif floor > block:
        print("Error: number of floors (stacking size) should be less than possible total number of floors")
    else:
        floor_next={}
        # initialization: build 1st floor
        if floor == 1:

            floor_next[(0,0)] = []
            for i in range(floor,m-block+floor+1):
                for j in range(floor,n-block+floor+1):
                    floor_next[(0,0)].append((i,j))

            return floor_next

        else:

            for f in build_possible_floors(m, n, floor-1, block)[(floor-2,floor-2)]:
                floor_next[f] = []
                for i in range(f[0]+1,m-block+floor+1):
                    for j in range(f[1]+1,n-block+floor+1):
                        floor_next[f].append((i,j))
        
            return floor_next

# print("rec start")
# floor = build_possible_floors(m, n, 3, 3)
# for k,v in floor.items():    
#     print(f'{k}:{v}')
# print("rec end")

def build_possible_stackings(m, n, block, stacking_size):
    if stacking_size > min(m,n):
        print("Error: total number of floors should be less than number of nts in the smallest strand")
    elif block > stacking_size:
        print("Error: number of floors (stacking size) should be less than possible total number of floors")
    else:
    
        all_floors = []
        inter_f = []

        if block == 1:
            all_floors = []
            for f in build_possible_floors(m, n, block, stacking_size)[(0,0)]:
                all_floors.append([f])
            return all_floors
        else:        
            floor_next = build_possible_floors(m, n, block, stacking_size)            
            for f in build_possible_stackings(m, n, block-1, stacking_size):
                inter_f = copy.deepcopy(f)
                for i in floor_next[f[-1]]:
                    inter_f.append(i)
                    all_floors.append(copy.deepcopy(inter_f))
                    inter_f.pop()

            return all_floors
#%%
sqnc_name = '3dvz'
MOTIF = read_sqnc(sqnc_name)
length = len(MOTIF)

# strands nts indices
a_start = 1
a_end = 7

b_start = a_end + 1
b_end = length

m, n = internal_stacks_sizes(MOTIF, a_start, a_end, b_start, b_end)

# set floor size to maximum
fs = min(m, n)
# print(fs)
# fs=5
floors = build_possible_stackings(m, n, fs, fs)
# print("rec start")
for pf in floors:
    print(pf)
# print("rec end")
# print(len(floors))

#%%
# convert floor to stacks
length = len(MOTIF)

def floor2stacks(floor, length):
    stacks_A_idx = [i[0] for i in floor]
    stacks_B_idx = [i[1] for i in floor]

    stacks_A = [(i, j) for i, j in zip(stacks_A_idx, stacks_A_idx[1:])]
    stacks_B = [(abs(i-length-1), abs(j-length-1)) for i, j in zip(stacks_B_idx, stacks_B_idx[1:])]

    return stacks_A, stacks_B

stacks_A, stacks_B = floor2stacks(floors[0], length)

# print(stacks_A)
# print(stacks_B)

# #%%
# # draw floors
# import networkx as nx
# from networkx.algorithms import bipartite
# import matplotlib.pyplot as plt


# for i in range(len(floors)):
    
#     G = nx.Graph()

#     PART_A = [x for x in range(a_start,a_end+1)]
#     PART_B = [x for x in range(b_start,b_end+1)]

#     # print(PART_A)
#     # print(PART_B)

#     G.add_nodes_from(PART_A, bipartite=0)
#     G.add_nodes_from(PART_B, bipartite=1)

#     edges = []

#     # print(floors[i])

#     for t in floors[i]:
#         edges.append((PART_A[t[0]-1], PART_B[t[1]-1]))
    
#     # when stacks without first and last nts 
#     # edges.append((PART_A[0], PART_B[0]))
#     # edges.append((PART_A[-1], PART_B[-1]))

#     G.add_edges_from(edges)

#     bipartite.is_bipartite(G)

#     nodes_A = [x for x in range(a_start,a_end+1)]
#     nodes_B = [x for x in range(b_start,b_end+1)]

#     labels_A = {k: MOTIF[k-1] for k in nodes_A}   
#     labels_B = {k: MOTIF[abs(k-length)+m] for k in nodes_B} # abs(k-1-length)+m+1? when stacks without first and last nts 
    
#     labels_A.update(labels_B)
#     labels=labels_A

#     plt.figure(i)
#     nx.draw_networkx(G, pos = nx.drawing.layout.bipartite_layout(G, PART_A), labels=labels, width = 1)
    
#     fig_name = sqnc_name
#     plt.savefig(f'{fig_name}_floors_{fs}_{i+1}.png')  


# %%
