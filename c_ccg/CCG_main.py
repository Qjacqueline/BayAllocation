# -*- coding: utf-8 -*-
# @Time    : 2025/4/20 3:10 PM
# @Author  : JacQ
# @File    : CCG_main.py
from itertools import permutations

from a_data_process.read_data import read_data
from c_ccg.CCG_multi import CCG_multi
from c_ccg.CCG_single import CCG_single


def generate_permutations(sequence, swapped=None):
    all_perms = permutations(sequence)
    valid_perms = [perm for perm in all_perms if valid_permutation(sequence, perm)]
    return valid_perms


def valid_permutation(sequence, perm):
    for i, elem in enumerate(perm):
        original_index = sequence.index(elem)
        if abs(original_index - i) > 1:
            return False
    return True


if __name__ == '__main__':
    print_flag = False
    ################### todo 注意40没有-1！！！！！
    inst_ls = [['CS', 'BNS', 6, 3, 10], ['CS', 'BNS', 8, 3, 10], ['CS', 'BNS', 10, 3, 10]]
    for ls in inst_ls:
        for i in range(10):
            C_type, B_type, group_num, block_num, miss_bay_num = ls
            inst_type = C_type + '_' + str(group_num) + '_' + B_type + '_' + str(block_num) + '_' + str(miss_bay_num)
            # file_name = '/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/' + \
            #             C_type + ' ' + B_type + '/' + inst_type + '.txt'
            file_name = 'C:\\Users\\admin\\PycharmProjects\\BayAllocation\\a_data_process\\data\\CCG\\' + C_type + '_' + str(
                group_num) + '_' + B_type + '_' + str(block_num) + '_' + \
                        str(miss_bay_num) + '_' + str(i) + '.txt'
            data = read_data(file_name)
            # ============== 生成所有可能提取序列 ================
            sequence = list(range(data.G_num))
            valid_permutations = generate_permutations(sequence, swapped=None)
            pi_num_set = [i for i in range(len(valid_permutations))]
            invalid_sequences_r = []
            # ============== 求解 ================
            if block_num == 1:
                CCG_single(data, B_type)
            else:
                CCG_multi(data, B_type,CCG_show_flag=True)
