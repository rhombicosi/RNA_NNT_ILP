import re
import sys

from score_matrix import stack_df

# first argument is a full path to a solution file 
fn = sys.argv[1]
#second argument is a file w/o extension with rna sequence
sqnc = sys.argv[2]

filepath = fn
pattern = "\((.*?)\)"
p_pattern = "^[P]\s"
q_pattern = "^[Q]\s"

p_score = 0

INFILE = open(sqnc, "r")

for RNA in INFILE:
    RNA = RNA.rstrip()
    print ("You input the sequence: ", RNA)
    print (RNA)

length = len(RNA)
print (length)
length_for_range = length + 1

q_weights = [0 for x in range(1, length_for_range + 1)]

# def find_stacks_score(RNA):

with open(filepath) as fp:
   line = fp.readline()
   cnt = 1
   while line:
       
       if "P" in line:
           pp = re.search(p_pattern, line)
           if pp!= None:
               p_line = line.split(' ')
               p_score = int(p_line[1])
               # print(p_line[1])  
       if "Q" in line:
           qq = re.search(q_pattern, line)
           if qq!= None:
               q_line = line.split(' ')
               q_score = int(q_line[1])
       
           
       # print("Line {}: {}".format(cnt, line.strip()))
       if " 1" in line and "Q(" in line:
            substring = re.search(pattern, line).group(1)
            indices = substring.split(",")
           
            i = int(indices[0]) - 1
            j = int(indices[1]) - 1
           
            q_weights[i] = stack_df.loc['>>', RNA[i] + RNA[i+1]]
            q_weights[j] = stack_df.loc['>>', RNA[j-1] + RNA[j]]
           
            print("Line {}: {}".format(cnt, line.strip()))
                     
       line = fp.readline()
       cnt += 1

print(f'Total pairs score: {p_score}')
print(f'Total stacks score: {q_score}')
print(f'Total score: {p_score + q_score}')
    # return sum(q_weights)
