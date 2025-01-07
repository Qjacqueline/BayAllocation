# -*- coding: utf-8 -*-
# @Time    : 2024/11/30 12:35 PM
# @Author  : JacQ
# @File    : DP_worst case.py
# -*- coding: utf-8 -*-
# @Time    : 2024/11/25 1:41 PM
# @Author  : JacQ
# @File    : test1.py
import random

import matplotlib.pyplot as plt

from itertools import permutations

from d_explore_sp_single_YC_min.DP import min_visit_cost, find_max_permutation_cost


def plot_points(A, B):
    N = len(A)

    # 创建图形
    plt.figure(figsize=(10, 2))

    # 绘制 A 点
    plt.scatter(A, [1] * N, color='blue', label='A', s=100)

    # 绘制 B 点
    plt.scatter(B, [2] * N, color='red', label='B', s=100)

    # 添加标签
    for i in range(N):
        plt.text(A[i], 1.02, f'A{i + 1}', fontsize=12, ha='center')
        plt.text(B[i], 2.02, f'B{i + 1}', fontsize=12, ha='center')

    # 设置图形属性
    plt.axhline(1, color='gray', linewidth=0.5, linestyle='--')  # A 点水平线
    plt.axhline(2, color='gray', linewidth=0.5, linestyle='--')  # B 点水平线
    # plt.title('A 和 B 点位置')
    # plt.xlabel('坐标')
    plt.xticks(range(min(min(A), min(B)), max(max(A), max(B)) + 1))
    plt.yticks([1, 2], ['A ', 'B'])
    plt.legend()
    plt.grid()
    plt.savefig("1.png")


def generate_initial_perm(pos):
    """
    基于启发式优先级生成初始排列。
    """
    N = len(pos)
    remaining = set(range(N))
    current_position = (0, 0)  # 初始位置
    initial_perm = []

    while remaining:
        best_priority = float('-inf')
        best_idx = None

        for idx in remaining:
            a_dist = abs(current_position[0] - pos[idx][0])
            b_dist = abs(current_position[1] - pos[idx][1])
            priority = min(a_dist, b_dist)  # 使用负号表示优先级越大越好

            if priority > best_priority:
                best_priority = priority
                best_idx = idx

        # 更新当前位置为选择的点对的中心点
        current_position = (
            (pos[best_idx][0] + pos[best_idx][1]) / 2,
            (pos[best_idx][0] + pos[best_idx][1]) / 2,
        )

        # 添加选择的点对到排列
        initial_perm.append(best_idx)
        remaining.remove(best_idx)

    return initial_perm


def iterative_max_path(N, pos, initial_perm, k):
    A, B = [pair[0] for pair in pos], [pair[1] for pair in pos]
    current_perm = initial_perm
    current_A = [A[i] for i in current_perm]
    current_B = [B[i] for i in current_perm]

    # 计算初始排列的路径
    current_cost, current_path, dp = min_visit_cost(N, current_A, current_B)

    # 将路径映射回原始点对的标识
    mapped_path = [(label, current_perm[index]) for label, index in current_path]
    print(f"初始代价: {current_cost}")
    print(f"初始路径: {' -> '.join(f'{label}{idx + 1}' for label, idx in mapped_path)}")

    for step in range(k):
        best_swap_cost = current_cost
        best_perm = current_perm

        # 遍历所有swap交换
        for i in range(N - 1):
            for q in range(N - 1):
                if i == q:
                    continue
                # 生成新排列
                new_perm = current_perm[:]
                new_perm[i], new_perm[q] = new_perm[q], new_perm[i]

                # 更新对应的 A 和 B
                new_A = [A[j] for j in new_perm]
                new_B = [B[j] for j in new_perm]

                # 计算新排列的路径代价
                new_cost, _, _ = min_visit_cost(N, new_A, new_B)

                if new_cost > best_swap_cost:
                    best_swap_cost = new_cost
                    best_perm = new_perm

        # 如果找到了更优排列，则更新
        if best_swap_cost > current_cost:
            current_cost = best_swap_cost
            current_perm = best_perm
            print(f"第 {step + 1} 次迭代: 更新代价为 {current_cost}")
        else:
            print(f"第 {step + 1} 次迭代: 无更优解，停止迭代")
            break

    # 最终路径
    final_A = [A[i] for i in current_perm]
    final_B = [B[i] for i in current_perm]
    _, final_path, _ = min_visit_cost(N, final_A, final_B)

    # 将最终路径映射回原始标识
    final_mapped_path = [(label, current_perm[index]) for label, index in final_path]
    return current_cost, current_perm, final_mapped_path


def iterative_max_path_3_swap(N, pos, initial_perm, k):
    A, B = [pair[0] for pair in pos], [pair[1] for pair in pos]
    current_perm = initial_perm
    current_A = [A[i] for i in current_perm]
    current_B = [B[i] for i in current_perm]

    # 计算初始排列的路径
    current_cost, current_path, _ = min_visit_cost(N, current_A, current_B)

    # 将路径映射回原始点对的标识
    mapped_path = [(label, current_perm[index]) for label, index in current_path]
    print(f"初始代价: {current_cost}")
    print(f"初始路径: {' -> '.join(f'{label}{idx + 1}' for label, idx in mapped_path)}")

    for step in range(k):
        best_swap_cost = current_cost
        best_perm = current_perm

        # 遍历所有三点交换
        for i in range(N - 2):
            for j in range(i + 1, N - 1):
                for q in range(j + 1, N):
                    # 生成新排列
                    new_perm = current_perm[:]
                    new_perm[i], new_perm[j], new_perm[q] = (
                        new_perm[q],
                        new_perm[i],
                        new_perm[j],
                    )

                    # 更新对应的 A 和 B
                    new_A = [A[t] for t in new_perm]
                    new_B = [B[t] for t in new_perm]

                    # 计算新排列的路径代价
                    new_cost, _, _ = min_visit_cost(N, new_A, new_B)

                    if new_cost > best_swap_cost:
                        best_swap_cost = new_cost
                        best_perm = new_perm

        # 如果找到了更优排列，则更新
        if best_swap_cost > current_cost:
            current_cost = best_swap_cost
            current_perm = best_perm
            print(f"第 {step + 1} 次迭代: 更新代价为 {current_cost}")
        else:
            print(f"第 {step + 1} 次迭代: 无更优解，停止迭代")
            break

    # 最终路径
    final_A = [A[i] for i in current_perm]
    final_B = [B[i] for i in current_perm]
    _, final_path, _ = min_visit_cost(N, final_A, final_B)

    # 将最终路径映射回原始标识
    final_mapped_path = [(label, current_perm[index]) for label, index in final_path]
    return current_cost, current_perm, final_mapped_path


# Dynamic Programming Implementation
def calculate_worst_case_path_no_track(N, pos):
    # Initialize DP table and trace table
    dp = [[[0] * 2 for _ in range(3)] for _ in range(N + 1)]
    path = [[[[] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 跟踪路径

    # Init state
    ## k=1
    dp[1][1][0] = abs(pos[0][1]) + abs(pos[0][1] - pos[0][0])  # B_1 -> A_1
    dp[1][1][1] = abs(pos[0][0]) + abs(pos[0][0] - pos[0][1])  # A_1 -> B_1
    path[1][1][0] = [(1, "B1->A1")]
    path[1][1][1] = [(1, "A1->B1")]

    if N > 1:
        dp[2][0][0] = abs(pos[1][1]) + abs(pos[1][0] - pos[1][1])  # A_2 -> B_2
        dp[2][0][1] = abs(pos[1][0]) + abs(pos[1][1] - pos[1][0])  # B_2 -> A_2
        path[2][0][0] = [(1, "B2->A2")]
        path[2][0][1] = [(1, "A2->B2")]

        ## k=2
        dp[1][2][0] = min(dp[2][0][0] + abs(pos[1][0] - pos[0][1]) + abs(pos[0][0] - pos[0][1]),
                          dp[2][0][1] + abs(pos[1][1] - pos[0][1]) + abs(pos[0][0] - pos[0][1]))
        dp[1][2][1] = min(dp[2][0][0] + abs(pos[1][0] - pos[0][0]) + abs(pos[0][0] - pos[0][1]),
                          dp[2][0][1] + abs(pos[1][1] - pos[0][0]) + abs(pos[0][0] - pos[0][1]))
        dp[2][1][0] = min(dp[1][1][0] + abs(pos[0][0] - pos[1][1]) + abs(pos[1][0] - pos[1][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[1][1]) + abs(pos[1][0] - pos[1][1]))
        dp[2][1][1] = min(dp[1][1][0] + abs(pos[0][0] - pos[1][0]) + abs(pos[1][0] - pos[1][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[1][0]) + abs(pos[1][0] - pos[1][1]))  # only case 1->2
        # 计算路径的更新
        path[1][2][0] = path[2][0][0] + [(2, f"B{1}->A{1}")] if dp[2][0][0] + abs(pos[1][0] - pos[0][1]) + abs(
            pos[0][0] - pos[0][1]) < dp[2][0][1] + abs(pos[1][1] - pos[0][1]) + abs(pos[0][0] - pos[0][1]) else \
            path[2][0][1] + [(2, f"B{1}->A{1}")]
        path[1][2][1] = path[2][0][0] + [(2, f"A{1}->B{1}")] if dp[2][0][0] + abs(pos[1][0] - pos[0][0]) + abs(
            pos[0][0] - pos[0][1]) < dp[2][0][1] + abs(pos[1][1] - pos[0][0]) + abs(pos[0][0] - pos[0][1]) else \
            path[2][0][1] + [(2, f"A{1}->B{1}")]
        path[2][1][0] = path[1][1][0] + [(2, f"B{2}->A{2}")] if dp[1][1][0] + abs(pos[0][0] - pos[1][1]) + abs(
            pos[1][0] - pos[1][1]) < dp[1][1][1] + abs(pos[0][1] - pos[1][1]) + abs(pos[1][0] - pos[1][1]) else \
            path[1][1][1] + [(2, f"B{2}->A{2}")]
        path[2][1][1] = path[1][1][0] + [(2, f"A{2}->B{2}")] if dp[1][1][0] + abs(pos[0][0] - pos[1][0]) + abs(
            pos[1][0] - pos[1][1]) < dp[1][1][1] + abs(pos[0][1] - pos[1][0]) + abs(pos[1][0] - pos[1][1]) else \
            path[1][1][1] + [(2, f"A{2}->B{2}")]

        ## k=3

    if N > 2:
        dp[3][0][0] = min(dp[1][1][0] + abs(pos[0][0] - pos[2][1]) + abs(pos[2][0] - pos[2][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[2][1]) + abs(pos[2][0] - pos[2][1]))
        dp[3][0][1] = min(dp[1][1][0] + abs(pos[0][0] - pos[2][0]) + abs(pos[2][0] - pos[2][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[2][0]) + abs(pos[2][0] - pos[2][1]))
        path[3][0][0] = path[1][1][0] + [(2, f"B{3}->A{3}")] if dp[1][1][0] + abs(pos[0][0] - pos[2][1]) + abs(
            pos[2][0] - pos[2][1]) < dp[1][1][1] + abs(pos[0][1] - pos[2][1]) + abs(pos[2][0] - pos[2][1]) else \
            path[1][1][1] + [(2, f"B{3}->A{3}")]
        path[3][0][1] = path[1][1][0] + [(2, f"A{3}->B{3}")] if dp[1][1][0] + abs(pos[0][0] - pos[2][0]) + abs(
            pos[2][0] - pos[2][1]) < dp[1][1][1] + abs(pos[0][1] - pos[2][0]) + abs(pos[2][0] - pos[2][1]) else \
            path[1][1][1] + [(2, f"A{3}->B{3}")]

    # DP transitions
    for i in range(3, N + 1):
        A_imm, B_imm = pos[i - 3]
        A_im, B_im = pos[i - 2]
        A_i, B_i = pos[i - 1]

        dp[i - 1][2][0] = min(dp[i][0][0] + abs(A_i - B_im) + abs(A_im - B_im),
                              dp[i][0][1] + abs(B_i - B_im) + abs(A_im - B_im))
        dp[i - 1][2][1] = min(dp[i][0][0] + abs(A_i - A_im) + abs(A_im - B_im),
                              dp[i][0][1] + abs(B_i - A_im) + abs(A_im - B_im))

        path[i - 1][2][0] = path[i][0][0] + [(i, f"B{i - 1}->A{i - 1}")] if dp[i][0][0] + abs(A_i - B_im) + abs(
            A_im - B_im) < dp[i][0][1] + abs(B_i - B_im) + abs(A_im - B_im) else path[i][0][1] + [
            (i, f"B{i - 1}->A{i - 1}")]
        path[i - 1][2][1] = path[i][0][0] + [(i, f"A{i - 1}->B{i - 1}")] if dp[i][0][0] + abs(A_i - A_im) + abs(
            A_im - B_im) < dp[i][0][1] + abs(B_i - A_im) + abs(A_im - B_im) else path[i][0][1] + [
            (i, f"A{i - 1}->B{i - 1}")]

        dp[i][1][0] = max(min(dp[i - 2][2][0] + abs(A_imm - B_i) + abs(A_i - B_i),
                              dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i)),
                          min(dp[i - 1][1][0] + abs(A_im - B_i) + abs(A_i - B_i),
                              dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i)))
        dp[i][1][1] = max(min(dp[i - 2][2][0] + abs(A_imm - A_i) + abs(A_i - B_i),
                              dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i)),
                          min(dp[i - 1][1][0] + abs(A_im - A_i) + abs(A_i - B_i),
                              dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i)))

        if min(dp[i - 2][2][0] + abs(A_imm - B_i) + abs(A_i - B_i),
               dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i)) > \
                min(dp[i - 1][1][0] + abs(A_im - B_i) + abs(A_i - B_i),
                    dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i)):
            path[i][1][0] = path[i - 2][2][0] + [(i, f"B{i} -> A{i}")] if dp[i - 2][2][0] + abs(A_imm - B_i) + abs(
                A_i - B_i) < dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i) else path[i - 2][2][1] + [
                (i, f"B{i} -> A{i}")]
        else:
            path[i][1][0] = path[i - 1][1][0] + [(i, f"B{i} -> A{i}")] if dp[i - 1][1][0] + abs(A_im - B_i) + abs(
                A_i - B_i) < dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i) else path[i - 1][1][1] + [
                (i, f"B{i} -> A{i}")]

        if min(dp[i - 2][2][0] + abs(A_imm - A_i) + abs(A_i - B_i),
               dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i)) > \
                min(dp[i - 1][1][0] + abs(A_im - A_i) + abs(A_i - B_i),
                    dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i)):
            # 如果选择 dp[i-2][2][0] 或 dp[i-2][2][1] 更优
            path[i][1][1] = path[i - 2][2][0] + [(i, f"A{i} -> B{i}")] if dp[i - 2][2][0] + abs(A_imm - A_i) + abs(
                A_i - B_i) < dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i) else path[i - 2][2][1] + [
                (i, f"A{i} -> B{i}")]
        else:
            # 如果选择 dp[i-1][1][0] 或 dp[i-1][1][1] 更优
            path[i][1][1] = path[i - 1][1][0] + [(i, f"A{i} -> B{i}")] if dp[i - 1][1][0] + abs(A_im - A_i) + abs(
                A_i - B_i) < dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i) else path[i - 1][1][1] + [
                (i, f"A{i} -> B{i}")]

        if i < N:
            dp[i + 1][0][0] = max(min(dp[i - 1][1][0] + abs(A_im - B_i) + abs(A_i - B_i),
                                      dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i)),
                                  min(dp[i - 2][2][0] + abs(A_imm - B_i) + abs(A_i - B_i),
                                      dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i)))

            dp[i + 1][0][1] = max(min(dp[i - 1][1][0] + abs(A_im - A_i) + abs(A_i - B_i),
                                      dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i)),
                                  min(dp[i - 2][2][0] + abs(A_imm - A_i) + abs(A_i - B_i),
                                      dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i)))
            if min(dp[i - 1][1][0] + abs(A_im - B_i) + abs(A_i - B_i),
                   dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i)) < \
                    min(dp[i - 2][2][0] + abs(A_imm - B_i) + abs(A_i - B_i),
                        dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i)):
                path[i + 1][0][0] = path[i - 2][2][0] + [(i, f"B{i + 1} -> A{i + 1}")] if dp[i - 2][2][0] + abs(
                    A_imm - B_i) + abs(
                    A_i - B_i) < dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i) else path[i - 2][2][1] + [
                    (i, f"B{i + 1} -> A{i + 1}")]
            else:
                path[i + 1][0][0] = path[i - 1][1][0] + [(i, f"B{i + 1} -> A{i + 1}")] if dp[i - 1][1][0] + abs(
                    A_im - B_i) + abs(
                    A_i - B_i) < dp[i - 1][1][1] + abs(B_im - B_i) + abs(A_i - B_i) else path[i - 1][1][1] + [
                    (i, f"B{i + 1} -> A{i + 1}")]
            if min(dp[i - 1][1][0] + abs(A_im - A_i) + abs(A_i - B_i),
                   dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i)) < \
                    min(dp[i - 2][2][0] + abs(A_imm - A_i) + abs(A_i - B_i),
                        dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i)):
                # 如果从 dp[i-2][2][0] 或 dp[i-2][2][1] 更优
                path[i + 1][0][1] = path[i - 2][2][0] + [(i, f"B{i + 1} -> A{i + 1}")] if dp[i - 2][2][0] + abs(
                    A_imm - A_i) + abs(
                    A_i - B_i) < dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i) else path[i - 2][2][1] + [
                    (i, f"B{i + 1} -> A{i + 1}")]
            else:
                # 如果从 dp[i-1][1][0] 或 dp[i-1][1][1] 更优
                path[i + 1][0][1] = path[i - 1][1][0] + [(i, f"B{i + 1} -> A{i + 1}")] if dp[i - 1][1][0] + abs(
                    A_im - A_i) + abs(
                    A_i - B_i) < dp[i - 1][1][1] + abs(B_im - A_i) + abs(A_i - B_i) else path[i - 1][1][1] + [
                    (i, f"B{i + 1} -> A{i + 1}")]

    max_value = max(min(dp[N][1][0], dp[N][1][1]), min(dp[N - 1][2][0], dp[N - 1][2][1]))
    index = [dp[N][1][0], dp[N][1][1], dp[N - 1][2][0], dp[N - 1][2][1]].index(max_value)
    candi_path = [path[N][1][0], path[N][1][1], path[N - 1][2][0], path[N - 1][2][1]]
    worst_path = candi_path[index]

    return max_value, worst_path, dp, path


# Dynamic Programming Implementation
def calculate_worst_case_path_no_track_n(N, pos):
    # Initialize DP table and trace table
    dp = [[[0] * 2 for _ in range(3)] for _ in range(N + 1)]
    path = [[[[] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 跟踪路径

    # Init state
    ## k=1
    dp[1][1][0] = abs(pos[0][1]) + abs(pos[0][1] - pos[0][0])  # B_1 -> A_1
    dp[1][1][1] = abs(pos[0][0]) + abs(pos[0][0] - pos[0][1])  # A_1 -> B_1

    if N > 1:
        dp[2][0][0] = abs(pos[1][1]) + abs(pos[1][0] - pos[1][1])  # A_2 -> B_2
        dp[2][0][1] = abs(pos[1][0]) + abs(pos[1][1] - pos[1][0])  # B_2 -> A_2

        ## k=2
        dp[1][2][0] = min(dp[2][0][0] + abs(pos[1][0] - pos[0][1]) + abs(pos[0][0] - pos[0][1]),
                          dp[2][0][1] + abs(pos[1][1] - pos[0][1]) + abs(pos[0][0] - pos[0][1]))
        dp[1][2][1] = min(dp[2][0][0] + abs(pos[1][0] - pos[0][0]) + abs(pos[0][0] - pos[0][1]),
                          dp[2][0][1] + abs(pos[1][1] - pos[0][0]) + abs(pos[0][0] - pos[0][1]))
        dp[2][1][0] = min(dp[1][1][0] + abs(pos[0][0] - pos[1][1]) + abs(pos[1][0] - pos[1][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[1][1]) + abs(pos[1][0] - pos[1][1]))
        dp[2][1][1] = min(dp[1][1][0] + abs(pos[0][0] - pos[1][0]) + abs(pos[1][0] - pos[1][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[1][0]) + abs(pos[1][0] - pos[1][1]))  # only case 1->2
    if N > 2:
        dp[3][0][0] = min(dp[1][1][0] + abs(pos[0][0] - pos[2][1]) + abs(pos[2][0] - pos[2][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[2][1]) + abs(pos[2][0] - pos[2][1]))
        dp[3][0][1] = min(dp[1][1][0] + abs(pos[0][0] - pos[2][0]) + abs(pos[2][0] - pos[2][1]),
                          dp[1][1][1] + abs(pos[0][1] - pos[2][0]) + abs(pos[2][0] - pos[2][1]))
        ## k=3
        dp[2][2][0] = min(dp[3][0][0] + abs(pos[2][0] - pos[1][1]) + abs(pos[1][0] - pos[1][1]),
                          dp[3][0][1] + abs(pos[2][1] - pos[1][1]) + abs(pos[1][0] - pos[1][1]))
        dp[2][2][1] = min(dp[3][0][0] + abs(pos[2][0] - pos[1][0]) + abs(pos[1][0] - pos[1][1]),
                          dp[3][0][1] + abs(pos[2][1] - pos[1][0]) + abs(pos[1][0] - pos[1][1]))
        dp[3][1][0] = [min(dp[1][2][0] + abs(pos[0][0] - pos[2][1]) + abs(pos[2][0] - pos[2][1]),
                           dp[1][2][1] + abs(pos[0][1] - pos[2][1]) + abs(pos[2][0] - pos[2][1])),
                       min(dp[2][1][0] + abs(pos[1][0] - pos[2][1]) + abs(pos[2][0] - pos[2][1]),
                           dp[2][1][1] + abs(pos[1][1] - pos[2][1]) + abs(pos[2][0] - pos[2][1]))]
        dp[3][1][1] = [min(dp[1][2][0] + abs(pos[0][0] - pos[2][0]) + abs(pos[2][0] - pos[2][1]),
                           dp[1][2][1] + abs(pos[0][1] - pos[2][0]) + abs(pos[2][0] - pos[2][1])),
                       min(dp[2][1][0] + abs(pos[1][0] - pos[2][0]) + abs(pos[2][0] - pos[2][1]),
                           dp[2][1][1] + abs(pos[1][1] - pos[2][0]) + abs(pos[2][0] - pos[2][1]))]
    if N > 3:
        dp[4][0][0] = [min(dp[1][2][0] + abs(pos[0][0] - pos[3][1]) + abs(pos[3][0] - pos[3][1]),
                           dp[1][2][1] + abs(pos[0][1] - pos[3][1]) + abs(pos[3][0] - pos[3][1])),
                       min(dp[2][1][0] + abs(pos[1][0] - pos[3][1]) + abs(pos[3][0] - pos[3][1]),
                           dp[2][1][1] + abs(pos[1][1] - pos[3][1]) + abs(pos[3][0] - pos[3][1]))]
        dp[4][0][1] = [min(dp[1][2][0] + abs(pos[0][0] - pos[3][0]) + abs(pos[3][0] - pos[3][1]),
                           dp[1][2][1] + abs(pos[0][1] - pos[3][0]) + abs(pos[3][0] - pos[3][1])),
                       min(dp[2][1][0] + abs(pos[1][0] - pos[3][0]) + abs(pos[3][0] - pos[3][1]),
                           dp[2][1][1] + abs(pos[1][1] - pos[3][0]) + abs(pos[3][0] - pos[3][1]))]

    # DP transitions
    for i in range(4, N + 1):
        A_imm, B_imm = pos[i - 3]
        A_im, B_im = pos[i - 2]
        A_i, B_i = pos[i - 1]

        dp[i - 1][2][0] = max(min(dp[i][0][0][0] + abs(A_i - B_im),
                                  dp[i][0][1][0] + abs(B_i - B_im)) + abs(A_im - B_im),
                              min(dp[i][0][0][1] + abs(A_i - B_im),
                                  dp[i][0][1][1] + abs(B_i - B_im)) + abs(A_im - B_im))
        dp[i - 1][2][1] = max(min(dp[i][0][0][0] + abs(A_i - A_im) + abs(A_im - B_im),
                                  dp[i][0][1][0] + abs(B_i - A_im) + abs(A_im - B_im)),
                              min(dp[i][0][0][1] + abs(A_i - A_im) + abs(A_im - B_im),
                                  dp[i][0][1][1] + abs(B_i - A_im) + abs(A_im - B_im)))

        dp[i][1][0] = [min(dp[i - 2][2][0] + abs(A_imm - B_i) + abs(A_i - B_i),
                           dp[i - 2][2][1] + abs(B_imm - B_i) + abs(A_i - B_i)),
                       max(min(dp[i - 1][1][0][0] + abs(A_im - B_i) + abs(A_i - B_i),
                               dp[i - 1][1][1][0] + abs(B_im - B_i) + abs(A_i - B_i)),
                           min(dp[i - 1][1][0][1] + abs(A_im - B_i) + abs(A_i - B_i),
                               dp[i - 1][1][1][1] + abs(B_im - B_i) + abs(A_i - B_i)))]
        dp[i][1][1] = [min(dp[i - 2][2][0] + abs(A_imm - A_i) + abs(A_i - B_i),
                           dp[i - 2][2][1] + abs(B_imm - A_i) + abs(A_i - B_i)),
                       max(min(dp[i - 1][1][0][0] + abs(A_im - A_i) + abs(A_i - B_i),
                               dp[i - 1][1][1][0] + abs(B_im - A_i) + abs(A_i - B_i)),
                           min(dp[i - 1][1][0][1] + abs(A_im - A_i) + abs(A_i - B_i),
                               dp[i - 1][1][1][1] + abs(B_im - A_i) + abs(A_i - B_i)))]
        if i < N:
            A_ip, B_ip = pos[i]
            dp[i + 1][0][0] = [min(dp[i - 2][2][0] + abs(A_imm - B_ip) + abs(A_ip - B_ip),
                                   dp[i - 2][2][1] + abs(B_imm - B_ip) + abs(A_ip - B_ip)),
                               max(min(dp[i - 1][1][0][0] + abs(A_im - B_ip) + abs(A_ip - B_ip),
                                       dp[i - 1][1][1][0] + abs(B_im - B_ip) + abs(A_ip - B_ip)),
                                   min(dp[i - 1][1][0][1] + abs(A_im - B_ip) + abs(A_ip - B_ip),
                                       dp[i - 1][1][1][1] + abs(B_im - B_ip) + abs(A_ip - B_ip)))]
            dp[i + 1][0][1] = [min(dp[i - 2][2][0] + abs(A_imm - A_ip) + abs(A_ip - B_ip),
                                   dp[i - 2][2][1] + abs(B_imm - A_ip) + abs(A_ip - B_ip)),
                               max(min(dp[i - 1][1][0][0] + abs(A_im - A_ip) + abs(A_ip - B_ip),
                                       dp[i - 1][1][1][0] + abs(B_im - A_ip) + abs(A_ip - B_ip)),
                                   min(dp[i - 1][1][0][1] + abs(A_im - A_ip) + abs(A_ip - B_ip),
                                       dp[i - 1][1][1][1] + abs(B_im - A_ip) + abs(A_ip - B_ip)))]

    max_value = max(max(min(dp[N][1][0][0], dp[N][1][1][0]), min(dp[N][1][0][1], dp[N][1][1][1])),
                    min(dp[N - 1][2][0], dp[N - 1][2][1]))
    # index = [dp[N][1][0], dp[N][1][1], dp[N - 1][2][0], dp[N - 1][2][1]].index(max_value)
    # candi_path = [path[N][1][0], path[N][1][1], path[N - 1][2][0], path[N - 1][2][1]]
    # worst_path = candi_path[index]
    # \max(\min(dp[N][1][0][0], dp[N][1][1][0]), \min(dp[N][1][0][1], dp[N][1][1][1]),)
    worst_path = []

    return max_value, worst_path, dp, path


if __name__ == '__main__':

    random.seed(1)
    # 示例使用
    N, range_min, range_max, pos = 5, -10, 10, []
    for _ in range(N):
        a = random.randint(range_min, range_max)
        b = random.randint(range_min, range_max)
        pos.append([a, b])

    k = 4  # 迭代步数

    # 启发式生成初始排列
    initial_perm = [i for i in range(N)]  # 初始排列
    # initial_perm = generate_initial_perm(pos)
    print(f"生成的初始排列: {initial_perm}")

    # # 迭代优化路径
    max_cost, final_perm, final_path = iterative_max_path(N, pos, initial_perm, k)

    # 格式化输出最终路径
    formatted_path = " -> ".join(f"{label}{index + 1}" for label, index in final_path)
    print(f"最大代价: {max_cost}")
    print(f"最终排列: {final_perm}")
    print(f"最终路径: {formatted_path}")
    print(f"***********************")

    # 枚举
    max_cost, best_permutation, best_path = find_max_permutation_cost(N, pos)

    # 格式化输出最佳路径
    formatted_path = " -> ".join(f"{label}{index + 1}" for label, index in best_path)

    print(f"2最大访问代价: {max_cost}")
    print(f"2最佳排列: {best_permutation}")
    print(f"2最佳路径: {formatted_path}")

    # plot_points([pair[0] for pair in pos], [pair[1] for pair in pos])

    print(f"***********************")
    # Run the function and print the result
    result, path_taken, dp, path = calculate_worst_case_path_no_track_n(N, pos)
    print("Path:", path_taken)
    print("Worst-case shortest path:", result)
