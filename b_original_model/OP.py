# -*- coding: utf-8 -*-
# @Time    : 2024/12/19 4:32 PM
# @Author  : JacQ
# @File    : OP.py
import time
from itertools import permutations, islice

from gurobipy import *
from matplotlib import pyplot as plt, patches

from a_data_process.read_data import read_data
import a_data_process.config as cf


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


def original_problem_robust(data, res=None, w_obj=None):
    """
        :param data: 数据
        :return:
    """
    s_t = time.time()
    # ============== 生成所有可能序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    # valid_permutations = [ (1, 0, 3, 2)]#[(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3), (1, 0, 2, 3), (1, 0, 3, 2)]
    pi_num = len(valid_permutations)

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
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "1c")
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
    model.addConstrs((C[w][uu] + big_M * (1 - Y[w][data.U_num][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][data.U_num][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k]) for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    model.addConstrs((C[w][uu] - C[w][u] + big_M * (1 - Y[w][u][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][u][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k])
                      for u in data.U for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    # Con9: z和x的关系
    model.addConstrs((2 * Z[w][u][uu][j - 1][jj - 1] <= X[u][j - 1] + X[uu][jj - 1] for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1o")
    model.addConstrs((Z[w][u][uu][j - 1][jj - 1] >= X[u][j - 1] + X[uu][jj - 1] - 1 for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1p")
    # con10:消除对称性 不同重量不等价，因为会有拼贝情况
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")
    # Con10: 优先级完成时间约束 fixme
    model.addConstrs((C[w][u] <= C[w][uu] for w in range(pi_num) for u in data.U for uu in data.U
                      if valid_permutations[w].index(data.U_g_set[u])
                      < valid_permutations[w].index(data.U_g_set[uu])), "1r")

    # ============== 构造目标 ================
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')  # 线性化模型变量
    if res is None:
        model.addConstrs((obj >= C_max[w] for w in range(pi_num)), "obj")
        # model.addConstr((obj >= C_max[7] ), "obj")
        model.setObjective(obj, gurobipy.GRB.MINIMIZE)
    else:
        for w in range(pi_num):
            model.setObjectiveN(C_max[w], gurobipy.GRB.MINIMIZE, 1)
        a = 1

    # ============== 求解参数 ================
    model.Params.OutputFlag = 0
    model.Params.timelimit = 1800
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
        # model.write('OP.lp')
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
            color_groups = generate_color_groups(data.G_num)
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 1  # 长方形宽度
            space = 0  # 长方形间隔
            for b in data.J:
                x_pos = b * (bar_width + space) - 1.5  # 长方形的x坐标
                if len(bay_x_dict[b]) == 1:
                    g = data.U_g_set[bay_x_dict[b][0]]
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor=color_groups[g], edgecolor='black')
                    ax.add_patch(rect)
                    ax.text(x_pos + bar_width / 2, 0.5, f'g{g}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                elif len(bay_x_dict[b]) == 2:
                    g1, g2 = data.U_g_set[bay_x_dict[b][0]], data.U_g_set[bay_x_dict[b][1]]
                    rect1 = patches.Rectangle((x_pos, 0.5), bar_width, 0.5, facecolor=color_groups[g1],
                                              edgecolor='black')
                    ax.add_patch(rect1)
                    # 下半部分
                    rect2 = patches.Rectangle((x_pos, 0), bar_width, 0.5, facecolor=color_groups[g2], edgecolor='black')
                    ax.add_patch(rect2)
                    # 在上半部分添加 g1 和 bay_x_dict[b][0]
                    ax.text(x_pos + bar_width / 2, 0.75, f'g{g1}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                    # 在下半部分添加 g2 和 bay_x_dict[b][1]
                    ax.text(x_pos + bar_width / 2, 0.25, f'g{g2}\nu{bay_x_dict[b][1]}',
                            ha='center', va='center', fontsize=8, color='black')
                else:
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor='grey', edgecolor='black')
                    ax.add_patch(rect)
            ax.set_xlim(-space, int(60 * (bar_width + space)) - space)
            ax.set_xticks([x for x in range(0, int(60 * (bar_width + space)), 2)])
            ax.set_xticklabels([str(x) for x in range(1, 60, 2)])

            ax.set_ylim(0, 1.2)
            ax.set_yticks([])
            ax.set_xlabel('Positions', fontsize=14)
            ax.set_title('Placement of Box Groups', fontsize=16)

            # 显示图形
            plt.savefig(case + ".png")
        # print("obj:", obj.X)
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
        return obj.X, res, worst_index


def original_problem_robust_test_P_allocation(data):
    """
        :param data: 数据
        :return:
    """
    s_t = time.time()

    # ============== 生成所有可能序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    # valid_permutations = [ (1, 0, 3, 2)]#[(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3), (1, 0, 2, 3), (1, 0, 3, 2)]
    pi_num = len(valid_permutations)

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
    res = {}
    bay_cnt = 0
    for i in data.U:
        if i in data.U_F:
            res.setdefault(i, bay_cnt + 2)
            bay_cnt += 4
        else:
            res.setdefault(i, bay_cnt)
            bay_cnt += 2

    model.addConstrs((X[u][res[u]] == 1 for u in data.U), "res")

    # Con1
    model.addConstrs((C[w][u] <= C_max[w] for w in range(pi_num) for u in data.U), "1b")
    # 对于一个子箱组
    # Con2: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "1c")
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
    model.addConstrs((C[w][uu] + big_M * (1 - Y[w][data.U_num][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][data.U_num][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k]) for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    model.addConstrs((C[w][uu] - C[w][u] + big_M * (1 - Y[w][u][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][u][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k])
                      for u in data.U for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    # Con9: z和x的关系
    model.addConstrs((2 * Z[w][u][uu][j - 1][jj - 1] <= X[u][j - 1] + X[uu][jj - 1] for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1o")
    model.addConstrs((Z[w][u][uu][j - 1][jj - 1] >= X[u][j - 1] + X[uu][jj - 1] - 1 for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1p")
    # con10:消除对称性 不同重量不等价，因为会有拼贝情况
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")
    # Con10: 优先级完成时间约束 fixme
    model.addConstrs((C[w][u] <= C[w][uu] for w in range(pi_num) for u in data.U for uu in data.U
                      if valid_permutations[w].index(data.U_g_set[u])
                      < valid_permutations[w].index(data.U_g_set[uu])), "1r")

    # ============== 构造目标 ================
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')  # 线性化模型变量
    model.addConstrs((obj >= C_max[w] for w in range(pi_num)), "obj")
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)
    # for w in range(pi_num):
    #     model.setObjectiveN(C_max[w], gurobipy.GRB.MINIMIZE, 1)

    # ============== 求解参数 ================
    model.Params.OutputFlag = 0
    model.Params.timelimit = 1800
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
        # model.write('OP.lp')
        # 输出结果为
        bay_x_dict = {b: [] for b in data.J}
        for u in range(data.U_num):
            for b in data.J:
                if abs(X[u][b - 1].X) > 0.00001:
                    bay_x_dict[b].append(u)
        if print_flag:
            # 画图
            color_groups = generate_color_groups(data.G_num)
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 1  # 长方形宽度
            space = 0  # 长方形间隔
            for b in data.J:
                x_pos = b * (bar_width + space) - 1.5  # 长方形的x坐标
                if len(bay_x_dict[b]) == 1:
                    g = data.U_g_set[bay_x_dict[b][0]]
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor=color_groups[g], edgecolor='black')
                    ax.add_patch(rect)
                    ax.text(x_pos + bar_width / 2, 0.5, f'g{g}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                elif len(bay_x_dict[b]) == 2:
                    g1, g2 = data.U_g_set[bay_x_dict[b][0]], data.U_g_set[bay_x_dict[b][1]]
                    rect1 = patches.Rectangle((x_pos, 0.5), bar_width, 0.5, facecolor=color_groups[g1],
                                              edgecolor='black')
                    ax.add_patch(rect1)
                    # 下半部分
                    rect2 = patches.Rectangle((x_pos, 0), bar_width, 0.5, facecolor=color_groups[g2], edgecolor='black')
                    ax.add_patch(rect2)
                    # 在上半部分添加 g1 和 bay_x_dict[b][0]
                    ax.text(x_pos + bar_width / 2, 0.75, f'g{g1}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                    # 在下半部分添加 g2 和 bay_x_dict[b][1]
                    ax.text(x_pos + bar_width / 2, 0.25, f'g{g2}\nu{bay_x_dict[b][1]}',
                            ha='center', va='center', fontsize=8, color='black')
                else:
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor='grey', edgecolor='black')
                    ax.add_patch(rect)
            ax.set_xlim(-space, int(60 * (bar_width + space)) - space)
            ax.set_xticks([x for x in range(0, int(60 * (bar_width + space)), 2)])
            ax.set_xticklabels([str(x) for x in range(1, 60, 2)])

            ax.set_ylim(0, 1.2)
            ax.set_yticks([])
            ax.set_xlabel('Positions', fontsize=14)
            ax.set_title('Placement of Box Groups', fontsize=16)

            # 显示图形
            plt.savefig(case + ".png")

        print("obj:", obj.X)
        print("time:", str(time.time() - s_t))
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
        # print("worst seq",valid_permutations[])
        return obj.X


def original_problem_stochastic(data, res=None, worst_seq_idx=None):
    """
        :param worst_seq_idx:
        :param res:
        :param data: 数据
        :return:
    """
    s_t = time.time()

    # ============== 生成所有可能序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    # valid_permutations = [ (1, 0, 3, 2)]#[(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3), (1, 0, 2, 3), (1, 0, 3, 2)]
    pi_num = len(valid_permutations)

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
        # res = {0: 0, 1: 4, 2: 8, 3: 10, 4: 20, 5: 12, 6: 14, 7: 16}
        model.addConstrs((X[u][res[u] - 1] == 1 for u in data.U), "res")

    # Con1
    model.addConstrs((C[w][u] <= C_max[w] for w in range(pi_num) for u in data.U), "1b")
    # 对于一个子箱组
    # Con2: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "1c")
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
    model.addConstrs((C[w][uu] + big_M * (1 - Y[w][data.U_num][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][data.U_num][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k]) for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    model.addConstrs((C[w][uu] - C[w][u] + big_M * (1 - Y[w][u][uu][k]) >=
                      data.U_num_set[uu] * cf.unit_process_time
                      + quicksum(Z[w][u][uu][j - 1][jj - 1] * cf.unit_move_time * abs(j - jj)
                                 for j in data.J_K[k] for jj in data.J_K[k])
                      for u in data.U for uu in data.U
                      for k in range(data.K_num) for w in range(pi_num)), "1n")
    # Con9: z和x的关系
    model.addConstrs((2 * Z[w][u][uu][j - 1][jj - 1] <= X[u][j - 1] + X[uu][jj - 1] for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1o")
    model.addConstrs((Z[w][u][uu][j - 1][jj - 1] >= X[u][j - 1] + X[uu][jj - 1] - 1 for u in data.U + [data.U_num]
                      for uu in data.U for j in data.J for jj in data.J for w in range(pi_num)), "1p")
    # con10:消除对称性 不同重量不等价，因为会有拼贝情况
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")
    # Con10: 优先级完成时间约束 fixme
    model.addConstrs((C[w][u] <= C[w][uu] for w in range(pi_num) for u in data.U for uu in data.U
                      if valid_permutations[w].index(data.U_g_set[u])
                      < valid_permutations[w].index(data.U_g_set[uu])), "1r")

    # ============== 构造目标 ================
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')  # 线性化模型变量
    if worst_seq_idx is not None:
        model.addConstr(obj >= C_max[worst_seq_idx], "obj")
    else:
        model.addConstr(obj >= sum(C_max[w] for w in range(pi_num)) / pi_num, "obj")
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)
    # for w in range(pi_num):
    #     model.setObjectiveN(C_max[w], gurobipy.GRB.MINIMIZE, 1)

    # ============== 求解参数 ================
    # model.Params.OutputFlag = 0
    model.Params.OutputFlag = 0
    model.Params.timelimit = 1800
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
        # model.write('OP.lp')
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
            color_groups = generate_color_groups(data.G_num)
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 1  # 长方形宽度
            space = 0  # 长方形间隔
            for b in data.J:
                x_pos = b * (bar_width + space) - 1.5  # 长方形的x坐标
                if len(bay_x_dict[b]) == 1:
                    g = data.U_g_set[bay_x_dict[b][0]]
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor=color_groups[g], edgecolor='black')
                    ax.add_patch(rect)
                    ax.text(x_pos + bar_width / 2, 0.5, f'g{g}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                elif len(bay_x_dict[b]) == 2:
                    g1, g2 = data.U_g_set[bay_x_dict[b][0]], data.U_g_set[bay_x_dict[b][1]]
                    rect1 = patches.Rectangle((x_pos, 0.5), bar_width, 0.5, facecolor=color_groups[g1],
                                              edgecolor='black')
                    ax.add_patch(rect1)
                    # 下半部分
                    rect2 = patches.Rectangle((x_pos, 0), bar_width, 0.5, facecolor=color_groups[g2], edgecolor='black')
                    ax.add_patch(rect2)
                    # 在上半部分添加 g1 和 bay_x_dict[b][0]
                    ax.text(x_pos + bar_width / 2, 0.75, f'g{g1}\nu{bay_x_dict[b][0]}',
                            ha='center', va='center', fontsize=8, color='black')
                    # 在下半部分添加 g2 和 bay_x_dict[b][1]
                    ax.text(x_pos + bar_width / 2, 0.25, f'g{g2}\nu{bay_x_dict[b][1]}',
                            ha='center', va='center', fontsize=8, color='black')
                else:
                    rect = patches.Rectangle((x_pos, 0), bar_width, 1, facecolor='grey', edgecolor='black')
                    ax.add_patch(rect)
            ax.set_xlim(-space, int(60 * (bar_width + space)) - space)
            ax.set_xticks([x for x in range(0, int(60 * (bar_width + space)), 2)])
            ax.set_xticklabels([str(x) for x in range(1, 60, 2)])

            ax.set_ylim(0, 1.2)
            ax.set_yticks([])
            ax.set_xlabel('Positions', fontsize=14)
            ax.set_title('Placement of Box Groups', fontsize=16)

            # 显示图形
            plt.savefig(case + ".png")
        # print("obj:", obj.X)
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
        # print("worst seq",valid_permutations[])
        return obj.X, res


def generate_color_groups(g):
    cmap = plt.cm.get_cmap('tab10', g)  # 从 'tab10' 调色板中生成 g 种颜色
    colors = [cmap(i) for i in range(g)]
    return colors


def prune_bays(data):
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
            pos_index = data.J_K[k].index(pos - 2)
            ls = data.J_K[k][pos_index + 1:].copy()
            for jj in ls:
                if jj in data.I:
                    data.I.remove(jj)
                data.J.remove(jj)
                data.J_K[k].remove(jj)

        # for j in range(len(data.J_K[k]) - 1, 0, -1):
        #     S_cnt, last = 0, -1
        #     for jj in range(j - 1, 0, -1):
        #         if jj in data.I:
        #             S_cnt += 1


if __name__ == '__main__':
    print_flag = False

    for case in ['case1', 'case2', 'case3', 'case4', 'case5', 'case6', 'case11']:
        dataa = read_data('/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/standard/' + case)
        prune_bays(dataa)
        obj1, res1, worst_index1 = original_problem_robust(dataa)
        #######测试不同随机模型#########
        obj2, res2, worst_index2 = original_problem_robust(dataa, res1, obj1)
        obj3, res3 = original_problem_stochastic(dataa)
        obj4, _ = original_problem_stochastic(dataa, res3, worst_index2)

        #######测试1.5P-allocation到底有多好#########
        obj5 = original_problem_robust_test_P_allocation(dataa)

        with open("/Users/jacq/PycharmProjects/BayAllocationGit/b_original_model/output.txt", "w") as f:
            f.write("This is a test output.\n")
            f.write("Second line of output.\n")
            f.write(f"{case}\tRobust:\t{obj1}\tStochastic:\t{obj3}\tStochastic-W:\t{obj4}\t1.5P_alloc:\t{obj5}\n")
