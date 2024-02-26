from score_matrix import bp_df, stack_df


score1 = stack_df.loc['>>','CG'] + stack_df.loc['>>','AA'] + stack_df.loc['>>','AA'] \
    + stack_df.loc['<>','GG']
    
# print(stack_df.loc['>>','AA'])
# print(stack_df.loc['<>','GG'])

score2 = stack_df.loc['>>','CU'] + stack_df.loc['>>','CC'] + stack_df.loc['>>','GU'] \
    + stack_df.loc['>>','UG'] + stack_df.loc['<>','CU'] + bp_df.loc['WW_cis','UU']
    

loop6_s = stack_df.loc['>>','UU'] + stack_df.loc['>>','AA'] + stack_df.loc['>>','UC'] \
    + stack_df.loc['>>','GU'] + stack_df.loc['>>','AG'] + stack_df.loc['>>','CG'] \
        + stack_df.loc['>>','CA']
        
loop6_q = stack_df.loc['<>','UA'] + stack_df.loc['>>','AC']

loop6_p = bp_df.loc['WW_cis','UU'] + bp_df.loc['SS_cis','AG']
 
    
print(loop6_s)
print(loop6_q)
print(loop6_p)
print(loop6_q + loop6_p)
print(loop6_q + loop6_p + loop6_s)