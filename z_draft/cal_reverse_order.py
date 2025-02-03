# -*- coding: utf-8 -*-
# @Time    : 2025/1/15 1:49 PM
# @Author  : JacQ
# @File    : cal_reverse_order.py
import itertools
from collections import defaultdict
from copy import deepcopy
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

    distance_ls = ["t"]
    for i in range(len(mapped_seq)):
        distance_ls.append("a" + str(bay_seq[i]))
        distance_ls.append("t")
        if i != len(mapped_seq) - 1:
            distance_ls.append("d" + str(i))
            distance_ls.append("t")

    resultt = []
    c_i = 0
    for i in range(len(mapped_seq)):
        if c_i <= mapped_seq[i] * 4 + 1:
            for term in distance_ls[c_i:mapped_seq[i] * 4 + 2]:
                if term != "t":
                    while True:
                        if term in resultt:
                            term += "'"
                        else:
                            break
                    resultt.append(term)
            c_i = mapped_seq[i] * 4 + 2
        elif c_i > mapped_seq[i]:
            for term in reversed(distance_ls[mapped_seq[i] * 4 + 1:c_i]):
                if term != "t":
                    while True:
                        if term in resultt:
                            term += "'"
                        else:
                            break
                    resultt.append(term)
            c_i = mapped_seq[i] * 4
    for item in distance_ls:
        if item != "t":
            resultt.remove(item)

    # 计算逆序数对
    inversions = {}
    inversion_flag = [False for _ in bay_seq]
    for i in range(len(mapped_seq)):
        for j in range(i + 1, len(mapped_seq)):
            if mapped_seq[i] > mapped_seq[j] and not inversion_flag[mapped_seq[j]]:
                inversions[mapped_seq[i]] = mapped_seq[j]
        if mapped_seq[i] in inversions.keys():
            inversion_flag[inversions[mapped_seq[i]]] = True

    # # 合并规则：只保留最大值相同且最小值更小的对
    # merged = {}
    # for a, b in inversions:
    #     max_val = max(a, b)
    #     if max_val not in merged or min(a, b) < merged[max_val][0]:
    #         merged[max_val] = (min(a, b), max(a, b))
    # inversions = sorted(set(merged.values()), key=lambda x: (min(x), max(x)))
    # merged_min = {}
    # for a, b in inversions:
    #     min_val = min(a, b)
    #     if min_val not in merged_min or max(a, b) > merged_min[min_val][0]:
    #         merged_min[min_val] = (min(a, b), max(a, b))
    # result = []
    # cnt = 0
    # for pair in inversions.items():
    #     # if cnt != 0:
    #     #     result += " + "
    #     pair = [min(pair[0], pair[1]), max(pair[0], pair[1])]
    #     cnt += 1
    #     for ii in range(pair[0], pair[1] + 1):
    #         if ii != pair[1]:
    #             a_a = "a" + str(bay_seq[ii])
    #             while True:
    #                 if a_a in result:
    #                     a_a += "'"
    #                 else:
    #                     break
    #             result.append(a_a)
    #             result.append(a_a + "'")
    #             d_d = "d" + str(ii)
    #             while True:
    #                 if d_d in result:
    #                     d_d += "'"
    #                 else:
    #                     break
    #             result.append(d_d)
    #             result.append(d_d + "'")
    #         else:
    #             a_a = "a" + str(bay_seq[ii])
    #             while True:
    #                 if a_a in result:
    #                     a_a += "'"
    #                 else:
    #                     break
    #             result.append(a_a)
    #             result.append(a_a + "'")
    # if proc_seq[-1] != bay_seq[-1]:
    #     indx = bay_seq.index(proc_seq[-1])
    #     lls = bay_seq[indx:]
    #     for q in lls:
    #         a_variants = [item for item in result if item.startswith("a" + str(q))]
    #         a_variants.sort(key=lambda x: x.count("'"), reverse=True)
    #         result.remove(a_variants[0])
    #     for ii in range(indx, len(bay_seq) - 1):
    #         d_variants = [item for item in result if item.startswith("d" + str(ii))]
    #         d_variants.sort(key=lambda x: x.count("'"), reverse=True)
    #         result.remove(d_variants[0])
    # norm_result = sorted([it.replace("'", "") for it in result])
    # norm_resultt = sorted([it.replace("'", "") for it in resultt])
    # if norm_resultt != norm_result:
    #     a = 1
    # # 输出结果
    return set(resultt)


def remove_subsets(input_list):
    # 去重并初始化结果列表
    input_list_copy = [set(sub_ls) for sub_ls in input_list]
    unique_list = []
    for lst in input_list:
        if set(lst) not in unique_list:
            unique_list.append(set(lst))

    result = []
    for i, sublist in enumerate(unique_list):
        # 检查当前子列表是否是其他子列表的子集
        if not any(set(sublist).issubset(set(other)) for j, other in enumerate(unique_list) if i != j):
            result.append([input_list_copy.index(sublist), list(sublist)])
    return result


def compare_different(ls1, ls2):
    set1 = set(ls1)
    set2 = set(ls2)
    only_in_ls1 = set1 - set2
    only_in_ls2 = set2 - set1
    # if only_in_ls1 is not None and only_in_ls2 is None:
    #     return []
    return set(list(only_in_ls1) + [f"-{item}" for item in only_in_ls2])


def minus_3C(input_list):
    standard_set = set()
    for i in range(len(nums)):
        standard_set.add("a" + str(i) + "'")
        standard_set.add("a" + str(i) + "''")
        if i != len(nums) - 1:
            standard_set.add("d" + str(i) + "'")
            standard_set.add("d" + str(i) + "''")
    res = []
    for i in range(len(input_list)):
        tmp_res = list(set(input_list[i][1]) - standard_set) \
                  + ["-" + term for term in list(standard_set - set(input_list[i][1]))]
        res.append([input_list[i][0], tmp_res])
    return res


if __name__ == '__main__':
    # 生成0, 1, 2, 3的全排列
    nums = [0, 1, 2, 3, 4, 5]
    all_bay_permutations = generate_bay_permutations(nums)
    all_seq_permutations = generate_seq_permutations(nums)
    res_all = [[0 for _ in range(len(all_seq_permutations))] for _ in range(len(all_bay_permutations))]
    print(all_seq_permutations)

    print("***************************")
    for i in range(len(all_bay_permutations)):
        bay_seq = all_bay_permutations[i]
        if bay_seq == (2, 0, 1, 3, 4):
            a = 1
        for j in range(len(all_seq_permutations)):
            proc_seq = all_seq_permutations[j]
            res_all[i][j] = calculate_and_merge_inversions(bay_seq, proc_seq)
        # res_all[i] = [sublist for sublist in remove_subsets(res_all[i])]
        # res_all[i] = [sublist for sublist in minus_3C(res_all[i])]
        # print(bay_seq, "\t", res_all[i])
    grouped = defaultdict(list)
    for seq in all_bay_permutations:
        key = (seq.index(nums[-2]), seq.index(nums[-1]))  # 第3位和第4位作为键
        grouped[key].append(seq)
    test_indx = 0
    for key, group in grouped.items():
        print("***************************")
        print("test\t", group[test_indx])
        cnt = 0
        candidate_perm = []
        for bay in group:
            p_ls = res_all[all_bay_permutations.index(bay)]
            inter_set = [[] for _ in range(len(res_all[all_bay_permutations.index(group[test_indx])]))]
            optimal_flag = True

            # 从两个列表中删除相同的元素
            p_ls_copy = deepcopy(p_ls)
            o_ls_copy = deepcopy(res_all[all_bay_permutations.index(group[test_indx])])
            for i in range(len(p_ls_copy)):  p_ls_copy[i] = set(p_ls_copy[i])
            for i in range(len(o_ls_copy)): o_ls_copy[i] = set(o_ls_copy[i])
            common = []
            for i in range(len(p_ls_copy)):
                if p_ls_copy[i] in o_ls_copy:
                    common.append(p_ls_copy[i])
            for item in common:
                if item in o_ls_copy: o_ls_copy.remove(item)
                p_ls_copy.remove(item)
            for o_sub_ls in o_ls_copy:
                for p_sub_ls in p_ls_copy:
                    tmp_res = compare_different(p_sub_ls, o_sub_ls)  # p_sub_ls-o_subls
                    inv_tmp_res = compare_different(o_sub_ls, p_sub_ls)
                    tmp_res = [item.replace("'", "") for item in tmp_res]  # o_sub_ls-p_subls
                    inv_tmp_res = [item.replace("'", "") for item in inv_tmp_res]
                    # 根据前提不等式处理
                    cntt = res_all[all_bay_permutations.index(group[test_indx])].index(o_sub_ls)
                    if len(inv_tmp_res) == 0:
                        inter_set[cntt].append("equal")
                    elif all(not item.startswith('-') for item in tmp_res):
                        inter_set[cntt].append("p>o")
                        optimal_flag = False
                    elif all(not item.startswith('-') for item in inv_tmp_res):
                        inter_set[cntt].append("o>p")
                    else:
                        inter_set[cntt].append(tmp_res)  # 若满足则o dominate p, 选o
                        optimal_flag = False

            cnttt = 0
            for term1 in inter_set:
                if "p>o" in term1:
                    inter_set[cnttt] = "o dominate p"
                elif "equal" in term1 or term1 == []:
                    inter_set[cnttt] = "equal"
                elif all(tmp_term1 == "o>p" for tmp_term1 in term1):
                    inter_set[cnttt] = "p dominate o"

                cnttt += 1

            if all(tmp_term == "o dominate p" or tmp_term == "equal" or tmp_term == [] for tmp_term in inter_set):
                candidate_perm.append(group[test_indx])
            elif all(tmp_term == "p dominate o" or tmp_term == "equal" or tmp_term == [] for tmp_term in inter_set):
                candidate_perm.append("be dominated")
            elif any(tmp_term == "p dominate o" or tmp_term == [] for tmp_term in inter_set):
                candidate_perm.append("candidate")
            else:
                candidate_perm.append("-")

            print(bay, [str(term) + "/" for term in inter_set])  # inter_set[0]
            cnt += 1

        print("***************************")
        i = 0
        for bay in group:
            print(bay, "\t", candidate_perm[i])  # all_bay_permutations[i], "\t",
            i = i + 1

# if __name__ == '__main__':
#     # 生成0, 1, 2, 3的全排列
#     nums = [0, 1, 2, 3, 4, 5, 6]
#     all_bay_permutations = generate_bay_permutations(nums)
#     all_seq_permutations = generate_seq_permutations(nums)
#     res_all = [[0 for _ in range(len(all_seq_permutations))] for _ in range(len(all_bay_permutations))]
#     print(all_seq_permutations)
#
#     print("***************************")
#     for i in range(len(all_bay_permutations)):
#         bay_seq = all_bay_permutations[i]
#         if bay_seq == (2, 0, 1, 3, 4):
#             a = 1
#         for j in range(len(all_seq_permutations)):
#             proc_seq = all_seq_permutations[j]
#             res_all[i][j] = calculate_and_merge_inversions(bay_seq, proc_seq)
#         res_all[i] = [sorted(sublist) for sublist in remove_subsets(res_all[i])]
#         print(bay_seq, "\t", res_all[i])
#
#     print("***************************")
#     test_inst = 0
#     # test_inst = all_bay_permutations.index((4, 3, 2, 1, 0))
#     compare_res_all = []
#     compare_res_all_inverse = []
#     for i in range(len(res_all[test_inst])):
#         compare_res, compare_res_inverse = [], []
#         for ii in range(0, len(res_all[test_inst])):
#             compare_res.append(compare_different(res_all[test_inst][i], res_all[test_inst][ii]))
#             compare_res_inverse.append(compare_different(res_all[test_inst][ii], res_all[test_inst][i]))
#         print(i, compare_res)
#         compare_res_all.append(compare_res)
#         compare_res_all_inverse.append(compare_res_inverse)
#
#     print("***************************")
#     added_constraint = [{"a0", "-a4"}]  # {"a4", "-a2", "-d4"},{"a4", "-a2"}
#     print("added_constraint", added_constraint)
#     print("test\t", all_bay_permutations[test_inst])
#     cnt = 0
#     candidate_perm = []
#     for p_ls in res_all:
#         if all_bay_permutations[cnt] == (1, 0, 2, 3, 4):
#             a = 1
#         inter_set = [[] for _ in range(len(res_all[test_inst]))]
#
#         optimal_flag = True
#
#         # 从两个列表中删除相同的元素
#         p_ls_copy = deepcopy(p_ls)
#         o_ls_copy = deepcopy(res_all[test_inst])
#         for i in range(len(p_ls_copy)):  p_ls_copy[i] = set(p_ls_copy[i])
#         for i in range(len(o_ls_copy)): o_ls_copy[i] = set(o_ls_copy[i])
#         common = []
#         for i in range(len(p_ls_copy)):
#             if p_ls_copy[i] in o_ls_copy:
#                 common.append(p_ls_copy[i])
#         for item in common:
#             if item in o_ls_copy: o_ls_copy.remove(item)
#             p_ls_copy.remove(item)
#         for o_sub_ls in o_ls_copy:
#             for p_sub_ls in p_ls_copy:
#                 tmp_res = compare_different(p_sub_ls, o_sub_ls)  # p_sub_ls-o_subls
#                 inv_tmp_res = compare_different(o_sub_ls, p_sub_ls)
#                 tmp_res = [item.replace("'", "") for item in tmp_res]  # o_sub_ls-p_subls
#                 inv_tmp_res = [item.replace("'", "") for item in inv_tmp_res]
#                 # 根据前提不等式处理
#                 cntt = res_all[test_inst].index(o_sub_ls)
#                 for term1 in compare_res_all[cntt] + added_constraint:
#                     while True:
#                         if all(item in tmp_res for item in term1) and term1 != set():
#                             for item in term1:
#                                 tmp_res.remove(item)
#                         else:
#                             break
#                 for term1 in compare_res_all[cntt] + added_constraint:
#                     while True:
#                         if all(item in inv_tmp_res for item in term1) and term1 != set():
#                             for item in term1:
#                                 inv_tmp_res.remove(item)
#                         else:
#                             break
#                 # 记录o_sub_ls 如果是最大时，对其他序列其他的占优情况
#                 if len(inv_tmp_res) == 0:
#                     inter_set[cntt].append("equal")
#                 elif all(not item.startswith('-') for item in tmp_res):
#                     inter_set[cntt].append("p>o")
#                     optimal_flag = False
#                 elif all(not item.startswith('-') for item in inv_tmp_res):
#                     inter_set[cntt].append("o>p")
#                 else:
#                     inter_set[cntt].append(tmp_res)  # 若满足则o dominate p, 选o
#                     optimal_flag = False
#
#         cnttt = 0
#         for term1 in inter_set:
#             if "p>o" in term1:
#                 inter_set[cnttt] = "o dominate p"
#             elif "equal" in term1 or term1 == []:
#                 inter_set[cnttt] = "equal"
#             elif all(tmp_term1 == "o>p" for tmp_term1 in term1):
#                 inter_set[cnttt] = "p dominate o"
#
#             cnttt += 1
#
#         # if all(tmp_term == "o dominate p" for tmp_term in inter_set):
#         #     candidate_perm.append(all_bay_permutations[cnt])
#         # if all(tmp_term != "o dominate p" for tmp_term in inter_set):
#         #     candidate_perm.append(all_bay_permutations[cnt])
#         if all(tmp_term == "o dominate p" or tmp_term == "equal" or tmp_term == [] for tmp_term in inter_set):
#             candidate_perm.append(all_bay_permutations[test_inst])
#         elif all(tmp_term == "p dominate o" or tmp_term == "equal" or tmp_term == [] for tmp_term in inter_set):
#             candidate_perm.append("be dominated")
#         elif any(tmp_term == "p dominate o" or tmp_term == [] for tmp_term in inter_set):
#             candidate_perm.append("candidate")
#         else:
#             candidate_perm.append("-")
#
#         if cnt < 1000:
#             print(all_bay_permutations[cnt], [str(term) + "/" for term in inter_set])  # inter_set[0]
#         cnt += 1
#
#     print("***************************")
#     # seen = set()
#     # unique_list =[]
#     # for item in :
#     #     if item not in unique_list:
#     #         unique_list.append(item)
#     #         # seen.add(item)
#     # for tmp_ls in unique_list:
#     #     print(tmp_ls)
#     for i in range(len(candidate_perm)):
#         print(candidate_perm[i])  # all_bay_permutations[i], "\t",
