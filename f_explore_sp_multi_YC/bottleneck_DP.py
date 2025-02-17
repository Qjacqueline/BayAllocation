# -*- coding: utf-8 -*-
# @Time    : 2024/12/4 12:08 PM
# @Author  : JacQ
# @File    : DP.py
import random
from itertools import product, permutations, combinations, combinations_with_replacement

from matplotlib import pyplot as plt


def enumerate_sequences(G, K, pos):
    single_stage = list(product([0, 1], repeat=K))
    all_sequences = list(product(single_stage, repeat=G))
    flattened_sequences = [[[value for value in stage] for stage in sequence] for sequence in all_sequences]
    best_seq, best_v, best_schedule = None, float("inf"), None
    for seq in flattened_sequences:
        seq = [list(row) for row in zip(*seq)]
        v = [0 for _ in range(K)]
        r_v = [0 for _ in range(K)]
        c_pos = [0 for _ in range(K)]
        max_v = [0 for _ in range(G)]
        schedule = [[] for _ in range(K)]
        for g in range(G):
            for k in range(K):
                if seq[k][g] == 0:
                    v[k] = v[k] + max(0, abs(c_pos[k] - pos[k][g][1]) + r_v[k]) + abs(pos[k][g][1] - pos[k][g][0])
                    schedule[k].append([v[k] - abs(pos[k][g][1] - pos[k][g][0]) - abs(c_pos[k] - pos[k][g][1]),
                                        v[k] - abs(pos[k][g][1] - pos[k][g][0]), v[k]])
                    c_pos[k] = pos[k][g][0]
                else:
                    v[k] = v[k] + max(0, abs(c_pos[k] - pos[k][g][0]) + r_v[k]) + abs(pos[k][g][1] - pos[k][g][0])
                    schedule[k].append([v[k] - abs(pos[k][g][1] - pos[k][g][0]) - abs(c_pos[k] - pos[k][g][0]),
                                        v[k] - abs(pos[k][g][1] - pos[k][g][0]), v[k]])
                    c_pos[k] = pos[k][g][1]
            max_v[g] = max(v)
            r_v = [v[k] - max_v[g] for k in range(K)]
            v = [max_v[g] for _ in range(K)]
        if best_v > max_v[-1]:
            best_v = max_v[-1]
            best_seq = [seq]
            best_schedule = [schedule]
        if best_v == max_v[-1]:
            # best_v = max_v[-1]
            best_seq.append(seq)
            best_schedule.append(schedule)
    print(best_v, best_seq[0])
    plot_gantt(best_schedule[0])
    return best_v, best_seq


def plot_gantt(data):
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan']
    data.reverse()
    k = len(data)
    for machine_idx, machine_data in enumerate(data):
        for job_idx, (s_t, a_t, f_t) in enumerate(machine_data):
            # Gray bar for idle time
            ax.barh(machine_idx, width=a_t - s_t, left=s_t, color='gray', edgecolor='black', height=0.5)
            # Colored bar for active time
            active_bar = ax.barh(machine_idx, width=f_t - a_t, left=a_t, color=colors[job_idx % len(colors)],
                                 edgecolor='black', height=0.5)
            # Add job index number to the colored bar
            ax.text(a_t + (f_t - a_t) / 2, machine_idx, str(job_idx + 1),
                    ha='center', va='center', color='white', fontsize=10, weight='bold')

    # Add labels and grid
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels([f'Machine {k - i - 1}' for i in range(len(data))])
    ax.set_xlabel('Time')
    ax.set_title('Gantt Chart')
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()


def single_DP(s, N, A, B):
    dp = [[0] * 2 for _ in range(N + 1)]
    path = [[0] * 2 for _ in range(N + 1)]  # 用于记录路径选择

    # 初始化
    if s == 1:
        dp[0][0], dp[0][1] = abs(B[0]) + abs(A[0] - B[0]), abs(A[0]) + abs(A[0] - B[0])
    else:
        dp[0][0], dp[0][1] = abs(A[0] - B[0]), abs(A[0] - B[0])

    for i in range(1, N):
        option1 = dp[i - 1][0] + abs(A[i - 1] - B[i]) + abs(B[i] - A[i])
        option2 = dp[i - 1][1] + abs(B[i - 1] - B[i]) + abs(B[i] - A[i])

        if option1 < option2:
            dp[i][0] = option1
            path[i][0] = 0
        else:
            dp[i][0] = option2
            path[i][0] = 1

        option3 = dp[i - 1][0] + abs(A[i - 1] - A[i]) + abs(A[i] - B[i])
        option4 = dp[i - 1][1] + abs(B[i - 1] - A[i]) + abs(A[i] - B[i])

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
            optimal_path.append(('A', s + i))
            optimal_path.append(('B', s + i))
        else:
            optimal_path.append(('B', s + i))
            optimal_path.append(('A', s + i))
        last_choice = path[i][last_choice]

    optimal_path.reverse()  # 反转路径

    return min_cost, optimal_path, dp


def filter_dp_dict(dp_dict):
    # 将结果存储在一个新的字典中
    filtered_dict = {}

    # 按照 machine_index 分组
    machine_paths = {}
    for key, value in dp_dict.items():
        machine_index = value[0]
        if machine_index not in machine_paths:
            machine_paths[machine_index] = []
        machine_paths[machine_index].append((key, value))

    # 对每组路径进行筛选
    for machine_index, paths in machine_paths.items():
        # 按路径长度降序排序
        paths.sort(key=lambda x: len(x[1][2]), reverse=True)

        # 用一个集合记录已保留路径的节点
        retained_paths = []
        filtered_keys = set()

        for key, value in paths:
            path_set = set(value[2])  # 将路径转化为集合

            # 检查是否与已保留路径重叠
            if all(not path_set.issubset(set(retained_path[1][2])) for retained_path in retained_paths):
                retained_paths.append((key, value))  # 保留当前路径
                filtered_keys.add(key)

        # 将筛选后的路径加入 filtered_dict
        for key in filtered_keys:
            filtered_dict[key] = dp_dict[key]

    return filtered_dict


def bottleneck_dp(G, K, pos):
    v = [[] for _ in range(int(G * (G - 1) / 2))]
    # 计算去头尾dp
    sequences = list(combinations_with_replacement(range(1, G + 1), 2))
    dp_dict = {}
    for s, e in sequences:
        if s == 1 and e == 3:
            a = 1
        max_cost = 0
        max_path = []
        for k in range(K):
            cost, path, _ = single_DP(s, e - s + 1, [pair[0] for pair in pos[k][s - 1:e]],
                                      [pair[1] for pair in pos[k][s - 1:e]])
            if cost > max_cost:
                max_machine = k
                max_cost = cost
                max_path = path
        dp_dict[s, e] = [max_machine, max_cost, max_path]
    # 去掉重合段
    filtered_dp_dict = filter_dp_dict(dp_dict)
    # 整合DP
    V = [dp_dict[1, i][1] for i in range(1, G + 1)]
    r = [dp_dict[1, i][0] for i in range(1, G + 1)]
    bt = [1 for _ in range(1, G + 1)]
    V.insert(0, 0), r.insert(0, 0), bt.insert(0, 0)

    for i in range(2, G + 1):
        for ii in range(1, i):
            # ii->i
            if V[ii] + dp_dict[ii + 1, i][1] > V[i]:
                if dp_dict[ii + 1, i][0] == r[ii]:
                    continue
                else:
                    V[i] = V[ii] + dp_dict[ii + 1, i][1]
                    r[i] = dp_dict[ii + 1, i][0]
                    bt[i] = ii
    print(V[-1])


if __name__ == '__main__':
    random.seed(1)
    # 示例使用
    G, K, range_min, range_max, pos = 6, 3, -15, 15, []  # G:color数 K:机器数
    for _ in range(K):
        tmp_pos = []
        for _ in range(G):
            a = random.randint(range_min, range_max)
            b = random.randint(range_min, range_max)
            tmp_pos.append([a, b])
        pos.append(tmp_pos)
    enumerate_sequences(G, K, pos)
    bottleneck_dp(G, K, pos)
