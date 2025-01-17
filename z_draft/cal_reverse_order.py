# -*- coding: utf-8 -*-
# @Time    : 2025/1/15 1:49 PM
# @Author  : JacQ
# @File    : cal_reverse_order.py
import itertools
from itertools import permutations


def valid_permutation(sequence, perm):
    for i, elem in enumerate(perm):
        original_index = sequence.index(elem)
        if abs(original_index - i) > 1:
            return False
    return True


def generate_seq_permutations(sequence):
    all_perms = permutations(sequence)
    valid_perms = [perm for perm in all_perms if valid_permutation(sequence, perm)]
    return valid_perms


def generate_bay_permutations(sequence):
    return list(permutations(sequence))


def calculate_and_merge_inversions(bay_seq, proc_seq):
    # 创建seq1的索引字典
    index_map = {value: idx for idx, value in enumerate(bay_seq)}

    # 将seq2按照seq1的顺序映射为索引序列
    mapped_seq = [index_map[num] for num in proc_seq]

    # 计算逆序数对
    inversions = []
    for i in range(len(mapped_seq)):
        for j in range(i + 1, len(mapped_seq)):
            if mapped_seq[i] > mapped_seq[j]:
                inversions.append((mapped_seq[i], mapped_seq[j]))

    # 对逆序数对进行合并
    # 按照对的最大值排序，然后按最小值排序
    inversions = sorted(set(inversions), key=lambda x: (max(x), min(x)))

    # 合并规则：只保留最大值相同且最小值更小的对
    merged = {}
    for a, b in inversions:
        max_val = max(a, b)
        if max_val not in merged or min(a, b) < merged[max_val][0]:
            merged[max_val] = (min(a, b), max(a, b))
    inversions = sorted(set(merged.values()), key=lambda x: (min(x), max(x)))
    merged_min = {}
    for a, b in inversions:
        min_val = min(a, b)
        if min_val not in merged_min or max(a, b) > merged_min[min_val][0]:
            merged_min[min_val] = (min(a, b), max(a, b))
    result = []
    cnt = 0
    for key, pair in merged_min.items():
        # if cnt != 0:
        #     result += " + "
        cnt += 1
        for ii in range(pair[0], pair[1] + 1):
            if ii != pair[1]:
                a_a = "a" + str(bay_seq[ii])
                while True:
                    if a_a in result:
                        a_a += "'"
                    else:
                        break
                result.append(a_a)
                d_d = "d" + str(ii + 1)
                while True:
                    if d_d in result:
                        d_d += "'"
                    else:
                        break
                result.append(d_d)
            else:
                a_a = "a" + str(bay_seq[ii])
                while True:
                    if a_a in result:
                        a_a += "'"
                    else:
                        break
                result.append(a_a)

    # 输出结果
    return result


def remove_subsets(input_list):
    # 去重并初始化结果列表
    unique_list = []
    for lst in input_list:
        if lst not in unique_list:
            unique_list.append(lst)

    result = []
    for i, sublist in enumerate(unique_list):
        # 检查当前子列表是否是其他子列表的子集
        if not any(set(sublist).issubset(set(other)) for j, other in enumerate(unique_list) if i != j):
            result.append(sublist)

    return result


if __name__ == '__main__':
    # 生成0, 1, 2, 3的全排列
    nums = [0, 1, 2, 3, 4]
    all_bay_permutations = generate_bay_permutations(nums)
    all_seq_permutations = generate_seq_permutations(nums)
    res = [[0 for _ in range(len(all_seq_permutations))] for _ in range(len(all_bay_permutations))]
    for i in range(len(all_bay_permutations)):
        bay_seq = all_bay_permutations[i]
        for j in range(len(all_seq_permutations)):
            proc_seq = all_seq_permutations[j]
            res[i][j] = calculate_and_merge_inversions(bay_seq, proc_seq)

        res[i] = [sorted(sublist) for sublist in remove_subsets(res[i])]
        print(bay_seq, "\t", res[i])
    print("***************************")
    cnt = 0
    for ls in res[1:]:
        inter_set = []
        for sublss in ls:
            flag = False
            for subls in res[0]:
                if set(subls).issubset(set(sublss)):
                    flag = True
            if flag:
                continue
            else:
                inter_set.append(sublss)
        cnt += 1
        print(all_bay_permutations[cnt], inter_set)
