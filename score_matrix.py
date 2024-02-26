import pandas as pd 

'''
each nucleotide has botom face, up face
'''


def sum_df_rows_columns(df):
    df = df.reset_index()
    df['Sum'] = df.sum(axis=1)

    column_list = list(df)
    
    column_list.remove('Desc')
    
    column_sum = []
    for c in column_list:
        df_l = df[c].drop_duplicates()
        column_sum.append(df_l.sum())
    
    column_sum_list = ['Sum']
    column_sum_list.extend(column_sum)
    
    df.loc[len(df)] = column_sum_list

# base pairings
bp_data = { 'Desc':['WW_cis','WW_tran','WH_cis','HW_cis','WH_tran','HW_tran','WS_cis','SW_cis','WS_tran',
                    'SW_tran','HH_cis','HH_tran','HS_cis','SH_cis','HS_tran','SH_tran','SS_cis','SS_tran'],
        'AA':[190, 674, 18, 18, 664, 664, 404, 404, 75, 75, 4, 1325, 315, 315, 447, 447, 6, 4],
        'AC':[753, 96, 18, 6, 14, 934, 713, 278, 52, 42, 3, 76, 11, 110, 210, 116, 5, None],
        'CA':[753, 96, 6, 18, 934, 14, 278, 713, 42, 52, 3, 76, 110, 11, 116, 210, 5, None],
        'AG':[1088, 16, 125, 91, 25, 26, 269, 11, 1184, 11, 4, 65, 102, 1, 7449, 1, 989, 3664],
        'GA':[1088, 16, 91, 125, 26, 25, 11, 269, 11, 1184, 4, 65, 1, 102, 1, 7449, 989, 3664],
        'AU':[18988, 527, 4, 849, 12, 3541, 195, 252, 3, 200, 8, 162, 119, 2, 91, 1, 5, None],
        'UA':[18988, 527, 849, 4, 3541, 12, 252, 195, 200, 3, 8, 162, 2, 119, 1, 91, 5, None],
        'CC':[628, 18, 128, 128, 77, 77, 271, 271, 29, 29, None, 2, 50, 50, 158, 158, None, None],
        'CG':[57896, 406, 65, 5, 90, 8, 79, 83, 156, 41, 5, 151, 6, None, 2, 1, 161, 20],
        'GC':[57896, 406, 5, 65, 8, 90, 83, 79, 41, 156, 5, 151, None, 6, 1, 2, 161, 20],
        'CU':[329, 59, 1, 31, 4, 3, 63, 19, 1, 4, None, 19, 182, 1, 237, 1, None, None],
        'UC':[329, 59, 31, 1, 3, 4, 19, 63, 4, 1, None, 19, 1, 182, 1, 237, None, None],
        'GG':[33, 25, 460, 460, 214, 214, 92, 92, 12, 12, None, 2, 95, 95, 213, 213, 72, 402],
        'GU':[8028, 116, 2, 141, 24, 11, 195, 121, 253, 120, None, None, None, 1179, 1, 307, 219, 50],
        'UG':[8028, 116, 141, 2, 11, 24, 121, 195, 120, 253, None, None, 1179, None, 307, 1, 219, 50],
        'UU':[1137, 112, 33, 33, 124, 124, 20, 20, 6, 6, None, None, 12, 12, 1, 1, None, None]}

bp_df = pd.DataFrame(bp_data)
bp_df = bp_df.set_index('Desc')

# sum_df_rows_columns(bp_df);

# print(bp_df)

# stacking
stack_data = {'Desc':['<<','<>','><','>>'],
        'AA':[13927, 6968, 1170, 13927],
        'AC':[7747, 2933, 1132, 13968],
        'CA':[13968, 2933, 1132, 7747],
        'AG':[12380, 10895, 1729, 16609],
        'GA':[16609, 10895, 1729, 12380],
        'AU':[3708, 2579, 985, 7141],
        'UA':[7141, 2579, 985, 3708],
        'CC':[16955, 203, 524, 16955],
        'CG':[19351, 9509, 943, 12464],
        'GC':[12464, 9509, 943, 19351],
        'CU':[8886, 313, 642, 11091],
        'UC':[11091, 313, 642, 8886],
        'GG':[27269, 10491, 2234, 27269],
        'GU':[8268, 3506, 1456, 13742],
        'UG':[13742, 3506, 1456, 8268],
        'UU':[4949, 301, 232, 4949]}

stack_df = pd.DataFrame(stack_data)
stack_df = stack_df.set_index('Desc')

# sum_df_rows_columns(stack_df)

# print(stack_df)

# base-phosphate
bph_data = {'Desc':['H_0BPh','S_1BPh','SW_2BPh','W_345BPh','W_6BPh','H_789BPh'],
        'AA':[256, None, 167, None, 1854, 254],
        'AC':[220, None, 113, None, 503, 298],
        'CA':[40, None, None, None, 300, 2907],
        'AG':[255, None, 842, None, 323, 250],
        'GA':[326, 839, None, 6126, None, None],
        'AU':[189, None, 324, None, 416, 206],
        'UA':[41, None, None, 746, None, 719],
        'CC':[67, None, None, None, 83, 1177],
        'CG':[130, None, None, None, 288, 1014],
        'GC':[369, 602, None, 916, None, None],
        'CU':[26, None, None, None, 75, 1510],
        'UC':[44, None, None, 430, None, 378],
        'GG':[251, 1024, None, 1289, None, None],
        'GU':[67, 500, None, 2251, None, None],
        'UG':[48, None, None, 950, None, 234],
        'UU':[17, None, None, 495, None, 317]}

bph_df = pd.DataFrame(bph_data)

# sum_df_rows_columns(bph_df)

# print(bph_df)

# base-ribose
br_data = {'Desc':['H_0BR','S_1BR','SW_2BR','W_345BR','W_6BR','H_789BR'],
        'AA':[1022, None, 1292, None, 938, 173],
        'AC':[846, None, 1048, None, 578, 335],
        'CA':[373, None, None, None, 85, 686],
        'AG':[2332, None, 1732, None, 5990, 1400],
        'GA':[662, 999, None, 252, None, None],
        'AU':[1058, None, 1337, None, 484, 355],
        'UA':[279, None, None, 75, None, 353],
        'CC':[442, None, None, None, 101, 416],
        'CG':[416, None, None, None, 214, 781],
        'GC':[499, 478, None, 368, None, None],
        'CU':[624, None, None, None, 160, 542],
        'UC':[360, None, None, 212, None, 261],
        'GG':[593, 1101, None, 856, None, None],
        'GU':[821, 583, None, 506, None, None],
        'UG':[731, None, None, 173, None, 536],
        'UU':[342, None, None, 129, None, 355]}

br_df = pd.DataFrame(br_data)

# sum_df_rows_columns(br_df)

# print(br_df)

