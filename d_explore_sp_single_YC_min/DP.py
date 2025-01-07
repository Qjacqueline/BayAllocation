# -*- coding: utf-8 -*-
# @Time    : 2025/1/7 8:56 AM
# @Author  : JacQ
# @File    : DP.py
from itertools import permutations


def min_visit_cost(N, A, B, pt, init_pos):
    dp = [[0] * 2 for _ in range(N + 1)]
    path = [[0] * 2 for _ in range(N + 1)]  # 用于记录路径选择

    # 初始化
    dp[0][0] = abs(B[0] - init_pos) + abs(B[0] - A[0]) + pt[0]  # 访问第 0 对时的代价
    dp[0][1] = abs(A[0] - init_pos) + abs(A[0] - B[0]) + pt[0]  # 访问第 0 对时的代价

    for i in range(1, N):
        option1 = dp[i - 1][0] + abs(A[i - 1] - B[i]) + abs(B[i] - A[i]) + pt[i]
        option2 = dp[i - 1][1] + abs(B[i - 1] - B[i]) + abs(B[i] - A[i]) + pt[i]

        if option1 < option2:
            dp[i][0] = option1
            path[i][0] = 0
        else:
            dp[i][0] = option2
            path[i][0] = 1

        option3 = dp[i - 1][0] + abs(A[i - 1] - A[i]) + abs(A[i] - B[i]) + pt[i]
        option4 = dp[i - 1][1] + abs(B[i - 1] - A[i]) + abs(A[i] - B[i]) + pt[i]

        if option3 < option4:
            dp[i][1] = option3
            path[i][1] = 0
        else:
            dp[i][1] = option4
            path[i][1] = 1

    # 最终结果为访问所有点的最小代价
    min_cost = min(dp[N - 1][0], dp[N - 1][1])
    last_choice = 0 if dp[N - 1][0] < dp[N - 1][1] else 1

    # 追踪路径
    optimal_path = []
    for i in range(N - 1, -1, -1):
        if last_choice == 0:
            optimal_path.append(('A', i))
            optimal_path.append(('B', i))
        else:
            optimal_path.append(('B', i))
            optimal_path.append(('A', i))
        last_choice = path[i][last_choice]

    optimal_path.reverse()  # 反转路径

    return min_cost, optimal_path, dp


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


def find_max_permutation_cost(N, pos, pt, init_pos):
    A, B = [pair[0] for pair in pos], [pair[1] for pair in pos]

    max_cost = float('-inf')
    best_permutation = None
    best_path = None
    sequence = list(range(N))
    valid_permutations = generate_permutations(sequence, swapped=None)

    for perm in valid_permutations:
        permuted_A = [A[i] for i in perm]
        permuted_B = [B[i] for i in perm]

        # 调用min_visit_cost获取当前排列的代价
        cost, path, dp = min_visit_cost(N, permuted_A, permuted_B, pt, init_pos)

        # 调整路径索引为原始排列索引
        adjusted_path = [(label, perm[idx]) for label, idx in path]

        if cost > max_cost:  # 更新最大代价
            max_cost = cost
            best_permutation = perm
            best_path = adjusted_path

    return max_cost, best_permutation, best_path
