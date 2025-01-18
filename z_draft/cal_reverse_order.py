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
    if bay_seq == (1, 4, 0, 5, 2, 3):
        a = 1
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
                result.append(a_a + "'")
                d_d = "d" + str(ii + 1)
                while True:
                    if d_d in result:
                        d_d += "'"
                    else:
                        break
                result.append(d_d)
                result.append(d_d + "'")
            else:
                a_a = "a" + str(bay_seq[ii])
                while True:
                    if a_a in result:
                        a_a += "'"
                    else:
                        break
                result.append(a_a)
                result.append(a_a + "'")
    if proc_seq[-1] != bay_seq[-1]:
        indx = bay_seq.index(proc_seq[-1])
        lls = bay_seq[indx:]
        for q in lls:
            a_variants = [item for item in result if item.startswith("a" + str(q))]
            a_variants.sort(key=lambda x: x.count("'"), reverse=True)
            result.remove(a_variants[0])
        for ii in range(indx, len(bay_seq) - 1):
            d_variants = [item for item in result if item.startswith("d" + str(ii + 1))]
            d_variants.sort(key=lambda x: x.count("'"), reverse=True)
            result.remove(d_variants[0])
    # 输出结果
    return result


def remove_subsets(input_list):
    # 去重并初始化结果列表
    unique_list = []
    for lst in input_list:
        if set(lst) not in unique_list:
            unique_list.append(set(lst))

    result = []
    for i, sublist in enumerate(unique_list):
        # 检查当前子列表是否是其他子列表的子集
        if not any(set(sublist).issubset(set(other)) for j, other in enumerate(unique_list) if i != j):
            result.append(list(sublist))

    return result


def compare_different(ls1, ls2):
    set1 = set(ls1)
    set2 = set(ls2)
    only_in_ls1 = set1 - set2
    only_in_ls2 = set2 - set1
    # if only_in_ls1 is not None and only_in_ls2 is None:
    #     return []
    return list(only_in_ls1) + [f"-{item}" for item in only_in_ls2]


if __name__ == '__main__':
    # 生成0, 1, 2, 3的全排列
    nums = [0, 1, 2, 3,4]
    all_bay_permutations = generate_bay_permutations(nums)
    all_seq_permutations = generate_seq_permutations(nums)
    res_all = [[0 for _ in range(len(all_seq_permutations))] for _ in range(len(all_bay_permutations))]
    print(all_seq_permutations)

    print("***************************")
    for i in range(len(all_bay_permutations)):
        bay_seq = all_bay_permutations[i]
        for j in range(len(all_seq_permutations)):
            proc_seq = all_seq_permutations[j]
            res_all[i][j] = calculate_and_merge_inversions(bay_seq, proc_seq)
        res_all[i] = [sorted(sublist) for sublist in remove_subsets(res_all[i])]
        print(bay_seq, "\t", res_all[i])

    print("***************************")
    test_inst = 0
    test_inst = all_bay_permutations.index((1, 0, 2, 4, 3))
    compare_res_all = []
    for i in range(len(res_all[test_inst])):
        compare_res = []
        for ii in range(0, len(res_all[test_inst])):
            compare_res.append(compare_different(res_all[test_inst][ii], res_all[test_inst][i]))
        print(i, compare_res)
        compare_res_all.append(compare_res)

    print("***************************")
    print("test\t", all_bay_permutations[test_inst])
    cnt = 0
    candidate_perm = []
    for p_ls in res_all:
        inter_set = []
        cntt = 0
        for o_subls in res_all[test_inst]:
            flag = False
            for p_sub_ls in p_ls:
                if set(o_subls).issubset(set(p_sub_ls)):
                    flag = True
                    # break
            if not flag:
                tmp_compare_set = [str(cntt)]
                optimal_flag = True
                for p_sub_ls in p_ls:
                    tmp_res = compare_different(p_sub_ls, o_subls)
                    if tmp_res in compare_res_all[cntt]:
                        tmp_compare_set.append("equal")
                    elif all(item.startswith('-') for item in tmp_res):
                        tmp_compare_set.append("dominated")
                    else:
                        tmp_compare_set.append(tmp_res)
                        optimal_flag = False
                inter_set.append(tmp_compare_set)
                if optimal_flag:
                    candidate_perm.append(all_bay_permutations[cnt])
                # break
            cntt += 1
        if cnt< 1000:
            print(all_bay_permutations[cnt], inter_set)
        cnt += 1

    print("***************************")
    # seen = set()
    unique_list = []
    for item in candidate_perm:
        if item not in unique_list:
            unique_list.append(item)
            # seen.add(item)
    print(unique_list)

