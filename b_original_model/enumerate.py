# -*- coding: utf-8 -*-
# @Time    : 2024/12/20 4:10 PM
# @Author  : JacQ
# @File    : enumerate.py
import itertools
from itertools import permutations, combinations
import a_data_process.config as cf
from a_data_process.read_data import read_data


def generate_g_permutations(sequence):
    def valid_g_permutation(seq, perm):
        for i, elem in enumerate(perm):
            original_index = seq.index(elem)
            if abs(original_index - i) > 1:
                return False
        return True

    all_perms = permutations(sequence)
    valid_perms = [perm for perm in all_perms if valid_g_permutation(sequence, perm)]
    return valid_perms


def is_valid_placement(plan, data):
    # constr1: 40TEU需要放在I型贝上
    forbi_bay = []
    for u, bay in plan.items():
        if u in data.U_F:
            if bay not in data.I or bay in forbi_bay:
                return False
            forbi_bay.append(bay + 2)
    return True


def generate_u_permutations(data):
    # 所有可能排列
    tmp_combinations = list(itertools.combinations(data.U_L, 2)) + list(itertools.combinations(data.U_F, 2))
    U_combinations = []
    # valid的排列
    for combi in tmp_combinations:
        if data.U_num_set[combi[0]] + data.U_num_set[combi[1]] > data.S_num * data.T_num:
            continue
        else:
            U_combinations.append(combi)
    U_combinations = U_combinations + [(u,) for u in data.U_L] + [(u,) for u in data.U_F]
    used_number = set()
    duplicate_number = set()
    for combi in U_combinations:
        if combi[0] in used_number:
            duplicate_number.add(combi[0])
        else:
            used_number.add(combi[0])
        if len(combi) == 2:
            if combi[1] in used_number:
                duplicate_number.add(combi[1])
            else:
                used_number.add(combi[1])
    not_duplicate_number = used_number - duplicate_number
    duplicate_combinations = set()
    for combi in U_combinations:
        if combi[0] in duplicate_number:
            duplicate_combinations.add(combi)
        else:
            if len(combi) == 2 and combi[1] in duplicate_number:
                duplicate_combinations.add(combi)
    unique_sequences = set()
    for r in range(int(len(duplicate_combinations) / 2), len(duplicate_combinations) + 1):
        for combo in combinations(duplicate_combinations, r):
            # 将组合中的元素合并成一个元组
            merged = tuple(set(item for sublist in combo for item in sublist))
            not_merged = [item for sublist in combo for item in sublist]
            if len(merged) == len(duplicate_number) and len(not_merged) == len(duplicate_number):
                unique_sequences.add(combo)
    unique_sequences = [list(seq) + [(i,) for i in used_number - duplicate_number] for seq in unique_sequences]
    return unique_sequences


def generate_binary_combinations(n):
    combinations = []
    for i in range(2 ** n):  # 从 0 到 2^n - 1
        # 将数字 i 转换为二进制，并填充到 n 位
        binary = format(i, f'0{n}b')
        combinations.append([int(bit) for bit in binary])
    return combinations


def calculate_extraction_time(pi_combi, plan, data):
    K_g_pos = [[[] for _ in range(data.G_num)] for _ in range(data.K_num)]
    for u, bay in plan.items():
        K_g_pos[data.J_K_dict[bay]][data.U_g_set[u]].append(u)

    g_AB_pos = [[[min(plan[u] for u in K_g_pos[k][g]), max(plan[u] for u in K_g_pos[k][g])]
                 for g in range(data.G_num)] for k in range(data.K_num)]
    # todo: 增加多场桥互相等
    AB_combis = generate_binary_combinations(data.G_num)
    best_v = float("inf")
    for AB_combi in AB_combis:
        cost = 0
        for k in range(data.K_num):
            pre_yc = data.J_K_first[k]
            for i in range(data.G_num):
                g_index = pi_combi[i]
                if g_AB_pos[k][g_index]:
                    if AB_combi[i] == 0:
                        cost += (abs(pre_yc - g_AB_pos[k][g_index][0]) + abs(
                            g_AB_pos[k][g_index][0] - g_AB_pos[k][g_index][1])) * cf.unit_move_time
                    else:
                        cost += (abs(pre_yc - g_AB_pos[k][g_index][1]) + abs(
                            g_AB_pos[k][g_index][0] - g_AB_pos[k][g_index][1])) * cf.unit_move_time
                    cost += cf.unit_process_time * sum(data.U_num_set[ii] for ii in K_g_pos[k][g_index])
        best_v = min(best_v, cost)
    return best_v


def calculate(data):
    pi_combi = list(range(data.G_num))
    pi_combinations = generate_g_permutations(pi_combi)
    global_opt_v = 0
    global_opt_plan = None
    # 枚举所有可能的优先级顺序
    for pi_combi in pi_combinations:
        local_opt_v = float("inf")
        local_opt_plan = None
        # 遍历所有可能的集装箱组合方案
        u_combinations = generate_u_permutations(data)
        for u_combi in u_combinations:
            plan = {}
            for bay_seq in permutations(data.J, len(u_combi)):
                for idx, bay in enumerate(bay_seq):
                    for u in u_combi[idx]:
                        plan[u] = bay
                if is_valid_placement(plan, data):
                    time = calculate_extraction_time(pi_combi, plan, data)
                    if time < local_opt_v:
                        local_opt_v = time
                        local_opt_plan = plan
        if local_opt_v > global_opt_v:
            global_opt_v = local_opt_v
            global_opt_plan = local_opt_plan
    print("obj", global_opt_v)


if __name__ == '__main__':
    dataa = read_data('/Users/jacq/PycharmProjects/BayAllocation/a_data_process/data/case1')
    calculate(dataa)
