import re
from stacking_floors import *
from motifs_lvls_444 import SCORE

sqnc_name = "5jf4"

# stacks
MOTIF = read_sqnc(sqnc_name)
length = len(MOTIF)

# strands nts indices
a_start = 1
a_end = 5

b_start = a_end + 1
b_end = length

m, n = internal_stacks_sizes(MOTIF, a_start, a_end, b_start, b_end)
fs = min(m, n)
floors = build_possible_stackings(m, n, fs, fs)
print(f'# of floors: {len(floors)}')

f_num = 2
stacks_A, stacks_B = floor2stacks(floors[f_num], length)
# print(stacks_A)
# print(stacks_B)

# file with solution
# fn = f'{sqnc_name}-motifs-dl.sol'

for i in range(1, 11):

    print(f'solution {i}')
    soln = i
    fn = f'{sqnc_name}-plist_{soln}.sol'

    # input sequence length
    with open(sqnc_name, 'r') as file:
        MOTIF = file.read().rstrip()
    lngth = len(MOTIF)

    filepath = fn
    pattern = "\((.*?)\)"

    pattern2 = r'\D\D\(\d+,\d+\)\s[1]'

    # fold = '.' * lngth
    fold = '(...()......)'

    with open(filepath) as fp:
        line = fp.readline()
        bp_list = []
        while line:

            substring = re.findall(pattern2, line)       

            if substring != []:
                bp_list.append(substring)

                bp = substring[0][0:2]
                indices = substring[0][2:-2]

                print(f'{bp}::{indices}')
                
            line = fp.readline()

    model_name = "motifs-bp"
    file_name = f'{sqnc_name}-{model_name}-{f_num}-{soln}-plist-stack.png'
    color = "#0000FF"
    stericity = "cis"
    color_string = f',color={color}'

    VARNA = open(f'{sqnc_name}-{model_name}-{soln}-plist-stack.bat', "w")

    string1 = f'java -cp ./VARNA.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "{MOTIF}"\
    -structureDBN "{fold}" -o {file_name} -auxBPs "'

    VARNA.write(string1)

    # write pairs
    for l in bp_list:
        edge3 = "wc" if l[0][0] == "W" else l[0][0].lower()
        edge5 = "wc" if l[0][1] == "W" else l[0][1].lower()
        indices = l[0][2:-2]
        
        if((l[0][0] == 'W' and l[0][1] == 'H') or (l[0][0] == 'H' and l[0][1] == 'S')):
            stericity = 'trans'
            print(indices[1:-1].split(','))
            index = indices[1:-1].split(',')
            i1 = int(index[0])
            i2 = int(index[1])
            if(i1<i2):                
                edge = edge3
                edge3 = edge5
                edge5 = edge
               
                # index = i1
                # i1 = i2
                # i2 = i1  
                # indices = ''.join(index)   
        else:
            stericity = 'cis'


        VARNA.write(f'{indices}:edge3={edge3},edge5={edge5}')
        VARNA.write(color_string)

        VARNA.write(f',stericity={stericity};')


    # write stacks
    for s in stacks_A:
        VARNA.write(f'({s[0]},{s[1]}):edge3=>>,edge5=>>,color=#E06666;')

    for s in stacks_B:
        VARNA.write(f'({s[1]},{s[0]}):edge3=>>,edge5=>>,color=#E06666;')

    VARNA.write('"')
    VARNA.close()