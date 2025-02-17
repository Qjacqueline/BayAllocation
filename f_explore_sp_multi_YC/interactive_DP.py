# -*- coding: utf-8 -*-
# @Time    : 2024/12/18 3:47 PM
# @Author  : JacQ
# @File    : CCG1.py
import random
import time
from itertools import combinations, permutations

import gurobipy as gp
from gurobipy import *
from a_data_process.read_data import read_data
import a_data_process.config as cf
from matplotlib import pyplot as plt, patches

zeta_cnt = 0
big_M = 10e4
eta_cnt = 0


def generate_permutations(sequence, swapped=None):
    if len(sequence) < 12:
        all_perms = permutations(sequence)
    else:
        all_perms = list(itertools.islice(itertools.permutations(sequence), 100000000))
    valid_perms = [perm for perm in all_perms if valid_permutation(sequence, perm)]
    return valid_perms


def valid_permutation(sequence, perm):
    for i, elem in enumerate(perm):
        original_index = sequence.index(elem)
        if abs(original_index - i) > 1:
            return False
    return True


def sub_problem_help(data, master_X):
    """
        :param data: 数据
        :param master_X: 主问题获得的解
        :return: N 几个箱组
                pos 对应A，B子箱组位置
    """
    if data.K_num == 1:
        G_u_pos = [[] for _ in range(data.G_num)]  # 每个箱组子箱组位置
        for u in range(data.U_num):
            j = master_X[u]
            G_u_pos[data.U_g_set[u]].append(j)
        pos = [[min(G_u_pos[g]) * cf.unit_move_time, max(G_u_pos[g]) * cf.unit_move_time]
               for g in range(data.G_num)]  # 每个箱组AB子箱组位置
        pt = [g_num * cf.unit_process_time for g_num in data.G_num_set]
        pt = [x for x in pt if x != 0]
    else:
        G_u_pos = [[[] for _ in range(data.G_num)] for _ in range(data.K_num)]  # 每个箱组子箱组位置
        for u in range(data.U_num):
            j = master_X[u]
            k_index = data.J_K_dict[j]
            G_u_pos[k_index][data.U_g_set[u]].append(j)

        pos = [[[min(G_u_pos[k][g]) * cf.unit_move_time if any(G_u_pos[k][g]) else None,
                 max(G_u_pos[k][g]) * cf.unit_move_time if any(G_u_pos[k][g]) else None]
                for g in range(data.G_num)] for k in range(data.K_num)]  # 每个箱组AB子箱组位置
        pt = [[sum(data.U_num_set[u] * cf.unit_process_time
                   if data.J_K_dict[master_X[u]] == k and data.U_g_set[u] == g else 0
                   for u in range(data.U_num)) for g in range(data.G_num)] for k in range(data.K_num)]
        for k in range(data.K_num):
            for i in range(data.G_num):
                if pt[k][i] == 0:
                    if i == 0:
                        pos[k][i] = [(data.J_K_first[k]) * cf.unit_move_time, (data.J_K_first[k]) * cf.unit_move_time]
                    else:
                        pos[k][i] = [None, None]
    return data.G_num, pos, pt, [(data.J_K_first[k]) * cf.unit_move_time for k in range(data.K_num)]


def sub_problem_multi(N, pos, pt, init_pos):
    # todo:顺序
    st_line = [0 for _ in range(data.G_num)]
    whole_schedule = [[] for _ in range(data.K_num)]
    touch_flag = [[False for _ in range(data.G_num)] for _ in range(data.K_num)]
    cnt = 0
    while True:
        cnt += 1
        stop_flag = [False for _ in range(data.K_num)]
        # calculate DP-T
        for k in range(data.K_num):
            schedule, min_cost = sub_problem_single_T(N=N, A=[pair[0] for pair in pos[k]],
                                                      B=[pair[1] for pair in pos[k]], pt=pt[k],
                                                      init_pos=init_pos[k], st_line=st_line, touch_flag=touch_flag[k])
            if whole_schedule[k] == schedule:
                stop_flag[k] = True
            else:
                whole_schedule[k] = schedule  # 创建一个计数器
        # terminate when no adjustment
        if all(stop_flag): break
        # adjustment
        st_line, touch_flag = [], [[False for _ in range(data.G_num)] for _ in range(data.K_num)]
        st = [[whole_schedule[k][g + 1] - whole_schedule[k][g] - pt[k][g + 1] for g in range(data.G_num - 1)] for k in
              range(data.K_num)]
        for g in range(data.G_num):
            max_l = max(whole_schedule[k][g] for k in range(data.K_num))
            for k in range(data.K_num):
                if whole_schedule[k][g] == max_l: touch_flag[k][g] = True
            st_line.append(max_l)
            if g == data.G_num - 1:
                break
            dt = [max_l - whole_schedule[k][g] for k in range(data.K_num)]
            # st = [whole_schedule[k][g + 1] - whole_schedule[k][g] - pt[k][g + 1] for k in range(data.K_num)]
            for k in range(data.K_num):
                whole_schedule[k][g + 1] = whole_schedule[k][g] + max(dt[k], st[k][g]) + pt[k][g + 1]
        a = 1
    print(cnt)
    return st_line[-1], whole_schedule


def sub_problem_single_T(N, A, B, pt, init_pos, st_line, touch_flag):
    dp = [[0] * 2 for _ in range(N + 1)]
    path = [[0] * 2 for _ in range(N + 1)]  # 用于记录路径选择

    def cal_dp_state(l, ii, ll, flag=True):
        pre_pre = A[ii - 1] if l == 0 else B[ii - 1]
        pre = A[ii] if ll == 0 else B[ii]
        suc = B[ii] if ll == 0 else A[ii]
        if ii == 0:
            return abs(pre - init_pos) + abs(pre - suc) + pt[0]
        else:
            if touch_flag[ii - 1]:
                return dp[ii - 1][l] + abs(pre_pre - pre) + abs(pre - suc) + pt[ii]
            else:
                return max(dp[ii - 1][l] + abs(pre_pre - pre), st_line[ii - 1]) + abs(pre - suc) + pt[ii]

    # 初始化
    dp[0][0] = cal_dp_state(0, 0, 1)  # abs(B[0] - init_pos) + abs(B[0] - A[0]) + pt[0]  # 访问第 0 对时的代价
    dp[0][1] = cal_dp_state(0, 0, 0)  # abs(A[0] - init_pos) + abs(A[0] - B[0]) + pt[0]  # 访问第 0 对时的代价
    path[0][0] = 0
    path[0][1] = 1
    for i in range(1, N):
        if A[i] is None:
            if touch_flag[i - 1]:
                dp[i][0] = dp[i - 1][0]
                dp[i][1] = dp[i - 1][1]
            else:
                dp[i][0] = max(dp[i - 1][0], st_line[i - 1])
                dp[i][1] = max(dp[i - 1][1], st_line[i - 1])
            path[i][0] = 0
            path[i][1] = 1
            A[i] = A[i - 1]
            B[i] = B[i - 1]
            continue
        option1 = cal_dp_state(0, i, 1)
        # max(dp[i - 1][1] + abs(B[i - 1] - B[i]), st_line[i - 1]) + abs(B[i] - A[i]) + pt[i]
        option2 = cal_dp_state(1, i, 1)
        if option1 < option2:
            dp[i][0] = option1
            path[i][0] = 0
        else:
            dp[i][0] = option2
            path[i][0] = 1

        # max(dp[i - 1][1] + abs(B[i - 1] - B[i]), st_line[i - 1]) + abs(B[i] - A[i]) + pt[i]
        option3 = cal_dp_state(0, i, 0)
        # max(dp[i - 1][1] + abs(B[i - 1] - A[i]), st_line[i - 1]) + abs(A[i] - B[i]) + pt[i]
        option4 = cal_dp_state(1, i, 0)

        if option3 < option4:
            dp[i][1] = option3
            path[i][1] = 0
        else:
            dp[i][1] = option4
            path[i][1] = 1
    min_cost = min(dp[N - 1][0], dp[N - 1][1])
    # min_cost = min(max(dp[N - 1][0], st_line[N - 1]), max(dp[N - 1][1], st_line[N - 1]))
    if dp[N - 1][0] < dp[N - 1][1]:
        last_choice = 0
    else:
        last_choice = 1
    optimal_path, schedule = [], []
    for i in range(N - 1, -1, -1):
        schedule.append(dp[i][last_choice])
        last_choice = path[i][last_choice]

    optimal_path.reverse()  # 反转路径
    schedule.reverse()

    return schedule, min_cost


def prune_bays():
    n_bay_num = len(data.U_F) * 2 + len(data.U_L)
    for k in range(data.K_num):
        # 正向
        find_flag = False
        for j in range(0, len(data.J_K[k]) - n_bay_num):
            cnt, pos = 1, data.J_K[k][j] + 2
            while True:
                if cnt == n_bay_num:
                    find_flag = True
                    break
                if pos in data.J_K[k]:
                    cnt += 1
                    pos += 2
                else:
                    break
            if find_flag:
                break
        if find_flag:
            pos_index = data.J_K[k].index(pos)
            ls = data.J_K[k][pos_index:].copy()
            for jj in ls:
                if jj in data.I:
                    data.I.remove(jj)
                data.J.remove(jj)
                data.J_K[k].remove(jj)


def generate_color_groups(g):
    cmap = plt.cm.get_cmap('tab10', g)  # 从 'tab10' 调色板中生成 g 种颜色
    colors = [cmap(i) for i in range(g)]
    return colors


def original_problem_robust(data, res=None, w_obj=None):
    """
        :param data: 数据
        :return:
    """
    s_t = time.time()

    # ============== 构造模型 ================
    big_M = 10000000
    model = Model("original problem")

    # ============== 定义变量 ================
    # X_uj: j从0开始 uN+k是虚拟job
    X = [[model.addVar(0, 1, vtype=GRB.BINARY, name='X_' + str(u) + "_" + str(j)) if j + 1 in data.J else 0
          for j in range(60 * data.K_num)] for u in range(data.U_num + 2)]

    # y^w_uuk : 每个场桥作业顺序 最后两个job分别为初始和末尾节点
    Y = [[[[model.addVar(0, 1, vtype=GRB.BINARY, name='Y_' + str(w) + "_" + str(u) + "_" + str(uu) + "_" + str(k))
            for k in range(data.K_num)] for uu in range(data.U_num + 2)] for u in range(data.U_num + 1)]
         for w in range(pi_num)]

    # z^w_uujj : 每个场桥作业顺序
    Z = [[[[[model.addVar(0, 1, vtype=GRB.BINARY,
                          name='Z_' + str(w) + "_" + str(u) + "_" + str(uu) + "_" + str(j) + "_" + str(jj))
             if j + 1 in data.J and jj + 1 in data.J else 0
             for jj in range(60 * data.K_num)] for j in range(60 * data.K_num)]
           for uu in range(data.U_num + 2)] for u in range(data.U_num + 2)] for w in range(pi_num)]

    # C^w_u : 完成时间
    C = [[model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name="C_" + str(w) + "_" + str(u))
          for u in range(data.U_num + 2)] for w in range(pi_num)]

    # Cmax^w: 最大完成时间
    C_max = [model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name="C_max_" + str(w)) for w in range(pi_num)]

    # ================ 约束 ==================
    # test
    if res is not None:
        model.addConstrs((X[u][res[u] - 1] == 1 for u in data.U), "res")

    # Con1
    model.addConstrs((C[w][u] <= C_max[w] for w in range(pi_num) for u in data.U), "1b")
    # 对于一个子箱组
    # Con2: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L + data.U_F), "1c")
    model.addConstrs((quicksum(X[u][j - 1] for j in data.I) == 1 for u in data.U_F), "1d")
    # con2: Initial position restrictions
    model.addConstrs((X[data.U_num][j - 1] == 1 for j in data.J_K_first), "1f")
    # con3:对于40ft的子箱组占了前一个后一个位置就要被虚拟子箱组占用
    model.addConstrs((X[u][j - 1] + X[uu][j - 3] <= 1 for u in data.U_F for uu in data.U for j in data.I), "1g")
    # Con4: 一个贝上放置的箱组一般不超过2个
    model.addConstrs((quicksum(X[u][j - 1] for u in data.U) <= 2 for j in data.J), "1h")
    # con5: 20和40的不能放在一个贝
    model.addConstrs((X[u][j - 1] + X[uu][j - 1] <= 1 for j in data.J for u in data.U_L for uu in data.U_F), "1i")
    # Con6: 放置集装箱数不超过贝位的最大容量
    model.addConstrs((quicksum(data.U_num_set[u] * X[u][j - 1] for u in data.U) <= data.S_num * data.T_num
                      for j in data.J), "1j")
    # Con7: x，y关系约束
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J_K[k]) + quicksum(X[uu][j - 1] for j in data.J_K[k])
                      + big_M * (1 - Y[w][u][uu][k]) >= 2 for k in range(data.K_num)
                      for u in data.U for uu in data.U for w in range(pi_num)), "1k")
    model.addConstrs((quicksum(Y[w][u][uu][k] for u in data.U + [data.U_num] if u != uu) -
                      quicksum(X[uu][j - 1] for j in data.J_K[k]) == 0 for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1l")
    model.addConstrs((quicksum(Y[w][u][uu][k] for uu in data.U + [data.U_num + 1] if u != uu) -
                      quicksum(X[u][j - 1] for j in data.J_K[k]) == 0 for u in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1m")
    model.addConstrs((quicksum(Y[w][u][data.U_num + 1][k] for u in data.U + [data.U_num]) == 1
                      for k in range(data.K_num) for w in range(pi_num)), "1l")
    model.addConstrs((quicksum(Y[w][data.U_num][uu][k] for uu in data.U + [data.U_num + 1]) == 1
                      for k in range(data.K_num) for w in range(pi_num)), "1m")
    # con8: 时序约束 pt+st
    for uu in data.U:
        for k in range(data.K_num):
            for w in range(pi_num):
                # 当 Y[w][data.U_num][uu][k] = 1 时的约束
                expr = data.U_num_set[uu] * cf.unit_process_time + \
                       gp.quicksum(Z[w][data.U_num][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                   for j in data.J_K[k] for jj in data.J_K[k])
                model.addConstr((Y[w][data.U_num][uu][k] == 1) >> (C[w][uu] >= expr), "1n")
    for u in data.U:
        for uu in data.U:
            for k in range(data.K_num):
                for w in range(pi_num):
                    # 当 Y[w][u][uu][k] = 1 时的约束
                    expr = data.U_num_set[uu] * cf.unit_process_time + \
                           gp.quicksum(Z[w][u][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                       for j in data.J_K[k] for jj in data.J_K[k])
                    model.addConstr((Y[w][u][uu][k] == 1) >> (C[w][uu] - C[w][u] >= expr), "1n")
    # Con9: z和x的关系
    model.addConstrs((2 * Z[w][u][uu][j - 1][jj - 1] <= X[u][j - 1] + X[uu][jj - 1] for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1o")
    model.addConstrs((Z[w][u][uu][j - 1][jj - 1] >= X[u][j - 1] + X[uu][jj - 1] - 1 for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1p")
    # con10:消除对称性 不同重量不等价，因为会有拼贝情况
    # model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
    #                   for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
    #                   and data.U_num_set[u] == data.U_num_set[uu]
    #                   and jj < j and u < uu), "1q")
    # Con10: 优先级完成时间约束 fixme
    model.addConstrs((C[w][u] <= C[w][uu] - data.U_num_set[uu] * cf.unit_process_time
                      for w in range(pi_num) for u in data.U for uu in data.U
                      if valid_permutations[w].index(data.U_g_set[u])
                      < valid_permutations[w].index(data.U_g_set[uu])), "1r")

    # ============== 构造目标 ================
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')  # 线性化模型变量
    model.addConstrs((obj >= C_max[w] for w in range(pi_num)), "obj")
    # model.addConstr((obj >= C_max[7] ), "obj")
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)
    model.setParam('MIPGap', 0.0001)

    # ============== 求解参数 ================
    model.Params.OutputFlag = 0
    model.Params.timelimit = 7200
    model.optimize()
    if model.status == GRB.Status.INFEASIBLE:
        print('Optimization was stopped with status %d' % model.status)
        # do IIS, find infeasible constraints
        model.computeIIS()
        model.write('a.ilp')
        for c in model.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
        return False
    else:
        # model.write('OP.sol')
        model.write('OP.lp')
        # 输出结果为
        bay_x_dict = {b: [] for b in data.J}
        res = {}
        for u in range(data.U_num):
            for b in data.J:
                if abs(X[u][b - 1].X) > 0.00001:
                    bay_x_dict[b].append(u)
                    res.setdefault(u, b)
        if print_flag:
            # 画图
            gap = 1.5
            color_groups = generate_color_groups(data.G_num)
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 1  # 长方形宽度
            space = 0  # 长方形间隔
            for b in data.J:
                x_pos = b % 60 * (bar_width + space) - 1.5  # 长方形的x坐标
                k = data.K_num - data.J_K_dict[b]
                if len(bay_x_dict[b]) == 1:
                    g = data.U_g_set[bay_x_dict[b][0]]
                    rect = patches.Rectangle((x_pos, k * gap), bar_width, 1, facecolor=color_groups[g],
                                             edgecolor='black')
                    ax.add_patch(rect)
                    ax.text(x_pos + bar_width / 2, 0.5 + k * gap, f'g{g}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                elif len(bay_x_dict[b]) == 2:
                    g1, g2 = data.U_g_set[bay_x_dict[b][0]], data.U_g_set[bay_x_dict[b][1]]
                    rect1 = patches.Rectangle((x_pos, 0.5 + k * gap), bar_width, 0.5, facecolor=color_groups[g1],
                                              edgecolor='black')
                    ax.add_patch(rect1)
                    # 下半部分
                    rect2 = patches.Rectangle((x_pos, k * gap), bar_width, 0.5, facecolor=color_groups[g2],
                                              edgecolor='black')
                    ax.add_patch(rect2)
                    # 在上半部分添加 g1 和 bay_x_dict[b][0]
                    ax.text(x_pos + bar_width / 2, 0.75 + k * gap, f'g{g1}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                    # 在下半部分添加 g2 和 bay_x_dict[b][1]
                    ax.text(x_pos + bar_width / 2, 0.25 + k * gap, f'g{g2}\nu{bay_x_dict[b][1]}',
                            ha='center', va='center', fontsize=8, color='black')
                else:
                    rect = patches.Rectangle((x_pos, k * gap), bar_width, 1, facecolor='grey', edgecolor='black')
                    ax.add_patch(rect)
            ax.set_xlim(-space, int(60 * (bar_width + space)) - space)
            ax.set_xticks([x for x in range(0, int(60 * (bar_width + space)), 2)])
            ax.set_xticklabels([str(x) for x in range(1, 60, 2)])

            ax.set_ylim(0, 3 * data.K_num)
            ax.set_yticks([])
            ax.set_xlabel('Positions', fontsize=14)
            ax.set_title('Placement of Box Groups', fontsize=16)

            # 显示图形
            plt.savefig(case + ".png")
        print("obj:", obj.X)
        # print("time:", str(time.time() - s_t))
        k = 0
        all_v = []
        for w in range(len(valid_permutations)):
            pre, v, path = data.U_num, 0, []
            for i in range(data.U_num):
                for u in data.U:
                    if abs(Y[w][pre][u][0].X - 1) <= 0.001:
                        suc = u
                        v = v + data.U_num_set[suc] * cf.unit_process_time + sum(
                            Z[w][pre][suc][j - 1][jj - 1].X * cf.unit_move_time * abs(j - jj)
                            for j in data.J_K[k] for jj in data.J_K[k])
                        pre = suc
                        path.append(v)
                        break
            all_v.append(path)

        worst_index = -1
        if w_obj is not None:
            for w in range(pi_num):
                if abs(w_obj - C_max[w].X) < 0.0001:
                    worst_index = w
                # if C_max[w].X > obj:
                #     obj = C_max[w].X
                #     worst_index = w
        # print("worst seq",valid_permutations[])
        G_u_C = [[[] for _ in range(data.G_num)] for _ in range(data.K_num)]  # 每个箱组子箱组位置
        for u in range(data.U_num):
            j = res[u]
            k_index = data.J_K_dict[j]
            G_u_C[k_index][data.U_g_set[u]].append(C[0][u].X)

        return obj.X, G_u_C, worst_index


def random_select_nums(lst, n):
    all_valid_combos = []

    def backtrack(start, current_combination):
        # 当当前组合的长度达到 n 时，将其添加到有效组合列表中
        if len(current_combination) == n:
            all_valid_combos.append(current_combination[:])
            return
        for i in range(start, len(lst)):
            # 若当前组合为空或者当前元素与组合中最后一个元素的差值大于 2
            if not current_combination or lst[i] - current_combination[-1] > 2:
                current_combination.append(lst[i])
                # 递归调用，继续寻找下一个可能的元素
                backtrack(i + 1, current_combination)
                current_combination.pop()

    backtrack(0, [])
    if all_valid_combos:
        # 若存在有效组合，随机选择一个返回
        return random.choice(all_valid_combos)
    return None


if __name__ == '__main__':
    case = 'case3m'
    random.seed(42)
    print_flag = True
    data = read_data('/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/standard/' + case)
    print(case)
    # ============== 生成所有可能提取序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = [generate_permutations(sequence, swapped=None)[0]]
    pi_num = len(valid_permutations)
    pi_num_set = [i for i in range(len(valid_permutations))]
    invalid_sequences_r = []
    # ============== 数据处理 ================
    prune_bays()
    bays = [1, 3, 11, 65, 91, 123, 71, 135, 187, 193, 213, 199, 203, 207, 183, 215, 143, 147]
    # bays=random_select_nums([i for i in data.J if i not in data.J_K_first], data.G_num)
    # random.shuffle(bays)
    res = {i: bays[i] for i in range(data.U_num)}
    N, pos, pt, init_pos = sub_problem_help(data, res)
    # ============== 求解 ================
    # OP求解
    st = time.time()
    obj1, res1, worst_index1 = original_problem_robust(data, res)
    print(obj1)
    for item in res1:
        print([max(tmp) if any(tmp) else 0 for tmp in item])
    print(time.time() - st)
    # DP求解
    stt = time.time()
    sub_results = sub_problem_multi(N, pos, pt, init_pos)
    print(sub_results[0])
    for item in sub_results[1]:
        print(item)
    print(time.time() - stt)
