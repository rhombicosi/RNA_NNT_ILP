import re
import sys

# first argument is a full path to a solution file 
fn = sys.argv[1]
#second argument is a rna sequence length
lngth = int(sys.argv[2])

filepath = fn
pattern = "\((.*?)\)"

fold = ["." for x in range(lngth)]

def find1(line, cnt, p_str):
    cnt = cnt
    while line:       
           
       if " 1" in line and p_str in line:
           print("Line {}: {}".format(cnt, line.strip()))
           substring = re.search(pattern, line).group(1)
           
           print(substring)
           indices = substring.split(",")
           
           i = int(indices[0])
           j = int(indices[1])
           fold[i - 1] = "("
           fold[j - 1] = ")"
           
       # if " 1" in line and "Q(" in line:
       #     print("Line {}: {}".format(cnt, line.strip()))
       line = fp.readline()
       # p_str = line.find("P(")
       cnt += 1
       
       
with open(filepath) as fp:
    line = fp.readline()
   
   # w_str = line.find("W(")
   # h_str = line.find("H(")
   # s_str = line.find("S(")
      
    cnt = 1
      
    while line:       
           
        if " 1" in line and ("W(" in line or "H(" in line or "S(" in line or "F(" in line or "K(" in line):
            print("Line {}: {}".format(cnt, line.strip()))
            substring = re.search(pattern, line).group(1)
            indices = substring.split(",")
           
            i = int(indices[0])
            j = int(indices[1])
            # fold[i - 1] = "("
            # fold[j - 1] = ")"
           
        # if " 1" in line and "W(" in line:
        #     print("Line {}: {}".format(cnt, line.strip()))
        #     substring = re.search(pattern, line).group(1)
        #     indices = substring.split(",")
           
        #     i = int(indices[0])
        #     j = int(indices[1])
        #     fold[i - 1] = "("
        #     fold[j - 1] = ")"
            
        line = fp.readline()
        # p_str = line.find("P(")
        cnt += 1

fold = ''.join([str(elem) for elem in fold])
print(fold)
