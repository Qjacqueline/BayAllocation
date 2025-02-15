# -*- coding: utf-8 -*-
# @Time    : 2024/12/18 3:47 PM
# @Author  : JacQ
# @File    : CCG1.py
import random
import time
from itertools import combinations, permutations

import gurobipy
from gurobipy import *
from a_data_process.read_data import read_data
import a_data_process.config as cf
from d_explore_sp_single_YC_min.DP import find_max_permutation_cost
from itertools import count

zeta_cnt = 0
big_M = 10e10
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


def CCG():
    LB = 0
    UB = float("inf")
    iteration_num = 0
    start_time = time.time()
    added_cuts = []
    # prune bays
    prune_bays()

    # 求解主问题
    master_results = init_master_problem()

    while 1:
        print("=========================================iteration\t" +
              str(iteration_num) + "=========================================")

        master_results = update_master_problem_adding_cuts(master_results[2], master_results[3], added_cuts)

        if not master_results:
            print("Master problem is infeasible!")
            return
        LB = max(LB, master_results[1])  # 当前迭代主问题目标值
        print("Master obj:\t", master_results[1], "\t", master_results[0])

        # 预处理子问题
        N, pos, pt, init_pos = sub_problem_help(data, master_results[0])

        # 求解子问题
        if master_results[0] == {0: 14, 1: 10, 2: 6, 3: 16, 4: 10, 5: 6}:
            a = 1
        if data.K_num == 1:
            sub_results = sub_problem_single(N, pos, pt, init_pos)
        else:
            sub_results = sub_problem_multi(N, pos, pt, init_pos)
        # compare_results = find_max_permutation_cost(N, pos, pt, init_pos)

        UB = min(UB, sub_results[0])

        # 添加cuts
        added_cuts = generate_cuts(master_results, sub_results, added_cuts)

        # 判断迭代终止条件
        if UB - LB < 0.00001:
            end_x = master_results[0]
            print("=========================================Final results=========================================")
            print("Final solution:\t", end_x)
            print("Final UB:\t", UB)
            print("Final LB:\t", LB)
            print("Time\t", time.time() - start_time)
            break
        else:
            end_time = time.time()
            print("UB:\t", UB)
            print("LB:\t", LB)
            print("Gap:\t", round((UB - LB) * 100 / UB, 1))
            print("Time:\t", end_time - start_time)
            if end_time - start_time >= 7200:
                print("=========================================Final results=========================================")
                print("Final UB:\t", UB)
                print("Final LB:\t", LB)
                print("Gap:\t", round((UB - LB) * 100 / UB, 1))
                print("Time\t", time.time() - start_time)
                break
        iteration_num += 1

    return UB


def update_master_problem_adding_cuts(model, variables, added_cuts):
    global zeta_cnt, eta_cnt
    X, Z, Theta, theta = variables
    # ============== 添加随机lb cut, 依据序列================
    if pi_num_set and random.random() < random_seed:
        r_c = random.choice(pi_num_set)
        pi_num_set.remove(r_c)
        P = []
        for i in range(len(valid_permutations[r_c]) - 1):
            P.append([valid_permutations[r_c][i], valid_permutations[r_c][i + 1]])
        min_var1 = [[model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                  name="min_var1" + '_' + str(k) + '_' + str(g) + '_' + str(r_c))
                     for g in range(data.G_num + 1)] for k in range(data.K_num)]
        model.addConstrs((min_var1[k][g] == gurobipy.min_(Z[l][ll][k][g][gg] for l in [0, 1] for ll in [0, 1])
                          for k in range(data.K_num) for g, gg in P), "tmp")
        model.addConstrs((min_var1[k][data.G_num] == gurobipy.min_(Theta[l][k][valid_permutations[r_c][0]]
                                                                   for l in [0, 1]) for k in range(data.K_num)), "tmp")
        model.addConstrs((theta >= cf.unit_move_time * (- data.J_K_first[k])
                          + cf.unit_move_time * sum(Theta[1][k][g] - Theta[0][k][g] for g in range(data.G_num))
                          + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1]
                                for u in data.U for j in data.J_K[k])
                          + sum(min_var1[k][g] for g in range(data.G_num))
                          + cf.unit_move_time * min_var1[k][data.G_num]
                          for k in range(data.K_num)), "global lb3")

    # ============== 添加iter cuts ================
    for cut in added_cuts:
        if cut[0] == 'Feasible':
            v, XX = cut[1], cut[2]
            rhs = LinExpr(0)
            for u, j in XX:
                rhs.addTerms(1, X[u][j])
            model.addConstr(theta >= big_M * (rhs - len(XX) + 1), "fea cut")
        elif cut[0] == 'Optimal':
            v, XX = cut[1], cut[2]
            rhs = LinExpr(0)
            for u, j in XX:
                rhs.addTerms(1, X[u][j])
            model.addConstr(theta >= v * (rhs - len(XX) + 1), "opt cut")
        # elif cut[0] == 'Translation':
        #     containers = [[term[0] for term in cut[2][k]] for k in range(data.K_num)]
        #     rhs1 = LinExpr(0)
        #     for k in range(data.K_num):
        #         for pair in cut[2][k]:
        #             rhs1.addTerms(1, X[pair[0]][pair[1]])
        #     eta = model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name='eta_' + str(eta_cnt))
        #     eta_cnt += 1
        #     model.addConstr(eta >= big_M * (rhs1 - data.G_num + 1), 'eta_' + str(eta_cnt))
        #     zeta = [[0 for q in range(len(cut[3][k]))] for k in range(data.K_num)]
        #     for k in range(data.K_num):
        #         for q in range(len(cut[3][k])):
        #             compose = cut[3][k][q]
        #             rhs2 = LinExpr(0)
        #             for cnt in range(len(compose)):
        #                 rhs2.addTerms(1, X[containers[k][cnt]][compose[cnt] - 1])
        #             zeta[k][q] = model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name='zeta_' + str(zeta_cnt))
        #             model.addConstr(zeta[k][q] + eta >= big_M * (rhs2 - data.K_num + 1), 'zeta_' + str(zeta_cnt))
        #             zeta_cnt += 1

        # elif cut[0] == 'Optimal delta':
        #     rhs = LinExpr(0)
        #     model.addConstr(theta >= big_M * (rhs - len(cut[2]) + 1), "opt cut2")

    # ============== 求解参数 ================
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
        # model.write('master.sol')
        # model.write('master.lp')
        master_X = {}
        for u in range(data.U_num):
            for j in range(cf.bay_number_one_block * data.K_num):
                if j + 1 not in data.J:
                    continue
                if abs(X[u][j].X) > 0.00001:
                    master_X[u] = j
        # master_X[data.U_num + 1] = []
        # for j in range(cf.bay_number_one_block * data.K_num):
        #     if j + 1 not in data.J:
        #         continue
        #     if abs(X[data.U_num + 1][j].X) > 0.00001:
        #         master_X[data.U_num + 1].append(j)
        return master_X, theta.X, model, (X, Z, Theta, theta)


def init_master_problem():
    """
        :return: 不可行，返回 FALSE；
                 可行，返回 master_X, theta.x
                   分别表示 解、theta值
    """
    # ============== 构造模型 ================

    model = Model("Master problem")

    # ============== 定义变量 ================
    # X_uj: j从0开始 uN+k是虚拟job
    X = [[[] for _ in range(cf.bay_number_one_block * data.K_num)] for _ in range(data.U_num + 2)]
    for j in range(0, cf.bay_number_one_block * data.K_num):
        for u in range(data.U_num + 2):
            if j + 1 not in data.J:
                X[u][j] = 0
                continue
            X[u][j] = model.addVar(0, 1, vtype=GRB.BINARY, name='X_' + str(u) + "_" + str(j))

    # theta: 箱组g在箱区k上最左A和最右B的位置
    Theta = [[[model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS,
                            name='theta_' + str(l) + "_" + str(k) + "_" + str(g))
               for g in range(data.G_num)] for k in range(data.K_num)] for l in [0, 1]]
    # Z^{ll'}_{kgg'}
    Z = [[[[[0 for gg in range(data.G_num)] for g in range(data.G_num)]
           for k in range(data.K_num)] for ll in [0, 1]] for l in [0, 1]]
    ZP = [[[[[0 for gg in range(data.G_num)] for g in range(data.G_num)]
            for k in range(data.K_num)] for ll in [0, 1]] for l in [0, 1]]
    for gg in range(data.G_num):
        for g in range(data.G_num):
            for k in range(data.K_num):
                for ll in [0, 1]:
                    for l in [0, 1]:
                        if abs(gg - g) <= 1 and g != gg:
                            Z[l][ll][k][g][gg] = model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                                              name='Z_' + str(l) + '_' + str(ll) + '_' + str(
                                                                  k) + '_' + str(g) + '_' + str(gg))
                            ZP[l][ll][k][g][gg] = model.addVar(-GRB.INFINITY, GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                                               name='ZP_' + str(l) + '_' + str(ll) + '_' + str(
                                                                   k) + '_' + str(g) + '_' + str(gg))
    eta = [[model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS,
                         name='eta_' + str(k) + '_' + str(g)) for g in range(data.G_num)] for k in range(data.K_num)]

    # res = {0: 5, 1: 3, 2: 7, 3: 11, 4: 13, 5: 17}
    # model.addConstrs((X[u][res[u]-1] == 1 for u in data.U), "res")

    # ================ 约束 ==================
    # 对于一个子箱组
    # con1: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "2b")
    model.addConstrs((quicksum(X[u][j - 1] for j in data.I) == 1 for u in data.U_F), "2c")
    model.addConstrs((quicksum(X[u][j - 1] for j in set(data.J) - set(data.I)) == 0 for u in data.U_F), "2c")
    # con2: Initial position restrictions
    model.addConstrs((X[data.U_num][j - 1] == 1 for j in data.J_K_first), "2d")
    # con3:对于40ft的子箱组占了前一个后一个位置就不能被其他使用
    model.addConstrs((X[u][j - 1] + X[uu][j - 3] <= 1 for u in data.U_F for uu in data.U for j in data.I), "2e")
    # 对于一个贝
    # con4: 一个贝上放置的箱组一般不超过2个
    model.addConstrs((quicksum(X[u][j - 1] for u in data.U) <= 2 for j in data.J), "2f")
    # con5: 20和40的不能放在一个贝
    model.addConstrs((X[u][j - 1] + X[uu][j - 1] <= 1 for j in data.J for u in data.U_L for uu in data.U_F), "2g")
    # con6: 放置集装箱数不超过贝位的最大容量
    model.addConstrs((quicksum(data.U_num_set[u] * X[u][j - 1] for u in data.U) <= data.S_num * data.T_num
                      for j in data.J), "2h")
    # con7: todo + big_M * (1 - quicksum(X[u][j - 1] for j in data.J_K[k]))
    tmpp_var = [[model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="tmpp_var" + str(u) + str(k)) for u in
                 data.U] for k in range(data.K_num)]
    model.addConstrs((tmpp_var[k][u] == quicksum(j * X[u][j - 1] for j in data.J_K[k])
                      for k in range(data.K_num) for u in data.U), "tmpp_var")
    # con8: 找A_kg, B_kg的位置
    model.addConstrs((Theta[0][k][g] == gurobipy.min_(tmpp_var[k][u] for u in data.U_g[g])
                      for k in range(data.K_num) for g in range(data.G_num)), "2i")
    model.addConstrs((Theta[1][k][g] == gurobipy.max_(tmpp_var[k][u] for u in data.U_g[g])
                      for k in range(data.K_num) for g in range(data.G_num)), "2k")
    # con9:
    model.addConstrs((ZP[l][ll][k][g][gg] == (Theta[l][k][g] - Theta[1 - ll][k][gg]) * cf.unit_move_time
                      for l in [0, 1] for ll in [0, 1]
                      for g in range(data.G_num) for gg in range(data.G_num)
                      for k in range(data.K_num) if abs(gg - g) <= 1 and g != gg), "2l")
    model.addConstrs((Z[l][ll][k][g][gg] == gurobipy.abs_(ZP[l][ll][k][g][gg])
                      for l in [0, 1] for ll in [0, 1] for g in range(data.G_num)
                      for gg in range(data.G_num) for k in range(data.K_num) if abs(gg - g) <= 1 and g != gg), "2m")
    model.addConstrs((eta[k][g] == (Theta[1][k][g] - Theta[0][k][g]) * cf.unit_move_time +
                      quicksum(X[u][j - 1] * data.U_num_set[u] * cf.unit_process_time
                               for u in data.U_g[g] for j in data.J_K[k])
                      for k in range(data.K_num) for g in range(data.G_num)), "2n")

    # min
    # min_var =[]

    # ============== 构造目标 ================
    theta = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='theta')  # 线性化模型变量
    obj = LinExpr(0)
    obj.addTerms(1, theta)
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)

    # ============== 添加global cuts ================
    # extreme bay indices
    AA = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='AA' + str(k)) for k in range(data.K_num)]
    BB = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='BB' + str(k)) for k in range(data.K_num)]
    model.addConstrs((AA[k] == gurobipy.min_(Theta[1][k][g] for g in range(data.G_num))
                      for k in range(data.K_num)), "BB")
    model.addConstrs((BB[k] == gurobipy.max_(Theta[1][k][g] for g in range(data.G_num))
                      for k in range(data.K_num)), "BB")

    # 必要时间估算1
    model.addConstrs((theta >= cf.unit_move_time * (AA[k] - data.J_K_first[k])
                      + cf.unit_move_time * sum(Theta[1][k][g] - Theta[0][k][g] for g in range(data.G_num))
                      + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1]
                            for u in data.U for j in data.J_K[k]) for k in range(data.K_num)), "global lb1")
    # 必要时间估算2
    model.addConstrs((theta >= cf.unit_move_time * (BB[k] - data.J_K_first[k])
                      + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1]
                            for u in data.U for j in data.J_K[k]) for k in range(data.K_num)), "global lb2")

    # 必要时间估算3
    model.addConstrs((theta >= cf.unit_move_time * 2 * BB[k] - cf.unit_move_time * AA[k]
                      - cf.unit_move_time * data.J_K_first[k]
                      + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1]
                            for u in data.U for j in data.J_K[k]) for k in range(data.K_num)), "global lb3")

    # 消除对称性：不同重量不等价，因为会有拼贝情况
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) + big_M * (sum(X[u][jjj - 1] for jjj in data.J) - 1)
                      + big_M * (sum(X[uu][jjj - 1] for jjj in data.J) - 1)
                      for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] != data.U_num_set[uu]
                      and jj < j and u < uu), "1q")

    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if ({u, uu}.issubset(data.U_L) or {u, uu}.issubset(data.U_F))
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")

    # 连续分配
    tau = [[0 for _ in range(0, cf.bay_number_one_block * data.K_num)],
           [0 for _ in range(0, cf.bay_number_one_block * data.K_num)]]
    for j in range(0, cf.bay_number_one_block * data.K_num):
        if j + 1 in data.J:
            tau[0][j] = model.addVar(vtype=GRB.BINARY, name='tau_0_' + str(j))
            tau[1][j] = model.addVar(vtype=GRB.BINARY, name='tau_1_' + str(j))
    model.addConstrs((sum(X[u][j - 1] for u in data.U_F) == X[data.U_num + 1][j - 3] for j in data.I),
                     "continuous cut0")  # fixme m是否成立
    model.addConstrs((BB[k] - (j - 1) <= big_M * tau[0][j - 1] for k in range(data.K_num)
                      for j in set(data.I) & set(data.J_K[k])), "continuous cut1")  # fixme m是否成立
    model.addConstrs(((j - 1) - AA[k] <= big_M * tau[1][j - 1] for k in range(data.K_num)
                      for j in set(data.I) & set(data.J_K[k])), "continuous cut2")  # fixme m是否成立
    model.addConstrs((tau[0][j - 1] + tau[1][j - 1] - 1
                      <= sum(X[u][j - 1] for u in data.U + [data.U_num + 1]) for j in data.J),
                     "continuous cut3")  # fixme m是否成立
    # ============== 求解参数 ================
    model.Params.OutputFlag = 0
    model.Params.timelimit = 3600
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
        model.write('master.sol')
        model.write('master.lp')
        master_X = {}
        for u in range(data.U_num):
            for j in range(cf.bay_number_one_block * data.K_num):
                if j + 1 not in data.J:
                    continue
                if abs(X[u][j].X) > 0.00001:
                    master_X[u] = j
        if master_X == {0: 30, 1: 28, 2: 16, 3: 10, 4: 6}:
            a = 1
        return master_X, theta.X, model, (X, Z, Theta, theta)


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
            k_index = data.J_K_dict[j + 1]
            G_u_pos[k_index][data.U_g_set[u]].append(j)

        pos = [[[min(G_u_pos[k][g]) * cf.unit_move_time if any(G_u_pos[k][g]) else None,
                 max(G_u_pos[k][g]) * cf.unit_move_time if any(G_u_pos[k][g]) else None]
                for g in range(data.G_num)] for k in range(data.K_num)]  # 每个箱组AB子箱组位置
        pt = [[data.G_num_set[u] * cf.unit_process_time if data.J_K_dict[master_X[u] + 1] == k else 0
               for u in range(data.U_num)] for k in range(data.K_num)]  #
        for k in range(data.K_num):
            for i in range(data.G_num):
                if pt[k][i] == 0:
                    if i == 0:
                        pos[k][i] = [0, 0]
                    else:
                        pos[k][i] = pos[k][i - 1]
    return data.G_num, pos, pt,   [(data.J_K_first[k]) * cf.unit_move_time for k in range(data.K_num)]


def sub_problem_multi(N, pos, pt, init_pos):
    # todo:顺序
    st_line = [0 for _ in range(data.G_num)]
    whole_schedule = [[] for _ in range(data.K_num)]
    touch_flag = [[False for _ in range(data.G_num)] for _ in range(data.K_num)]
    while True:
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
        for g in range(data.G_num):
            max_l = max(whole_schedule[k][g] for k in range(data.K_num))
            for k in range(data.K_num):
                if whole_schedule[k][g] == max_l: touch_flag[k][g] = True
            st_line.append(max_l)
            if g == data.G_num - 1:
                break
            dt = [max_l - whole_schedule[k][g] for k in range(data.K_num)]
            st = [whole_schedule[k][g + 1] - whole_schedule[k][g] - pt[k][g + 1] for k in range(data.K_num)]
            for k in range(data.K_num):
                whole_schedule[k][g + 1] = whole_schedule[k][g] + max(dt[k], st[k]) + pt[k][g + 1]
        a = 1
    return st_line[-1], whole_schedule


def sub_problem_single(N, pos, pt, init_pos):
    """
        :param N: 几个箱组
        :param pos: 对应A，B子箱组位置
        :param pt: 每个箱组所需加工时间
        :return: (1/0表示是否可行, 目标值)
             不可行 (0，float("inf"))
             可行 (1，目标值)
    """

    # dp[i][k][A,B][l node]
    def f_fuc(i, k, ii, jj):
        v1 = dp[i][k][0][0] + abs(pos[ii - 1][1 - jj] - pos[i - 1][0])
        v2 = dp[i][k][1][0] + abs(pos[ii - 1][1 - jj] - pos[i - 1][1])
        v3 = dp[i][k][0][1] + abs(pos[ii - 1][1 - jj] - pos[i - 1][0])
        v4 = dp[i][k][1][1] + abs(pos[ii - 1][1 - jj] - pos[i - 1][1])
        v = max(min(v1, v2), min(v3, v4)) + abs(pos[ii - 1][0] - pos[ii - 1][1]) + pt[ii - 1]
        if min(v1, v2) > min(v3, v4):
            pre = [i, k, 0, 0] if v1 < v2 else [i, k, 1, 0]
        else:
            pre = [i, k, 0, 1] if v3 < v4 else [i, k, 1, 1]
        return v, pre

    def g_fuc(i, k, ii, jj):
        v1 = dp[i][k][0][0] + abs(pos[ii - 1][1 - jj] - pos[i - 1][0])
        v2 = dp[i][k][1][0] + abs(pos[ii - 1][1 - jj] - pos[i - 1][1])
        v = min(v1, v2) + abs(pos[ii - 1][0] - pos[ii - 1][1]) + pt[ii - 1]
        pre = [i, k, 0, 0] if v1 < v2 else [i, k, 1, 0]
        return v, pre

    # Initialize DP table and trace table
    dp = [[[[0, 0] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 3:i,k,j | 4 i,k,j,h
    path = [[[[0, 0] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 跟踪路径

    # Init state
    ## k=1
    # init_pos=0
    dp[1][1][0][0] = abs(pos[0][1] - init_pos) + abs(pos[0][1] - pos[0][0]) + pt[0]  # B_1 -> A_1
    dp[1][1][1][0] = abs(pos[0][0] - init_pos) + abs(pos[0][0] - pos[0][1]) + pt[0]  # A_1 -> B_1
    path[1][1][0][0], path[1][1][1][0] = 0, 0

    if N > 1:
        dp[2][0][0][0] = abs(pos[1][1] - init_pos) + abs(pos[1][0] - pos[1][1]) + pt[1]  # A_2 -> B_2
        dp[2][0][1][0] = abs(pos[1][0] - init_pos) + abs(pos[1][1] - pos[1][0]) + pt[1]  # B_2 -> A_2
        path[2][0][0][0], path[2][0][1][0] = 0, 0

        ## k=2
        dp[1][2][0][0], path[1][2][0][0] = g_fuc(i=2, k=0, ii=1, jj=0)
        dp[1][2][1][0], path[1][2][1][0] = g_fuc(i=2, k=0, ii=1, jj=1)
        dp[2][1][0][0], path[2][1][0][0] = g_fuc(i=1, k=1, ii=2, jj=0)
        dp[2][1][1][0], path[2][1][1][0] = g_fuc(i=1, k=1, ii=2, jj=1)

    if N > 2:
        dp[3][0][0][0], path[3][0][0][0] = g_fuc(i=1, k=1, ii=3, jj=0)
        dp[3][0][1][0], path[3][0][1][0] = g_fuc(i=1, k=1, ii=3, jj=1)
        ## k=3
        dp[2][2][0][0], path[2][2][0][0] = g_fuc(i=3, k=0, ii=2, jj=0)
        dp[2][2][1][0], path[2][2][1][0] = g_fuc(i=3, k=0, ii=2, jj=1)
        dp[3][1][0][0], path[3][1][0][0] = g_fuc(i=1, k=2, ii=3, jj=0)
        dp[3][1][1][0], path[3][1][1][0] = g_fuc(i=1, k=2, ii=3, jj=1)
        dp[3][1][0][1], path[3][1][0][1] = g_fuc(i=2, k=1, ii=3, jj=0)
        dp[3][1][1][1], path[3][1][1][1] = g_fuc(i=2, k=1, ii=3, jj=1)

    if N > 3:
        dp[4][0][0][0], path[4][0][0][0] = g_fuc(i=1, k=2, ii=4, jj=0)
        dp[4][0][1][0], path[4][0][1][0] = g_fuc(i=1, k=2, ii=4, jj=1)
        dp[4][0][0][1], path[4][0][0][1] = g_fuc(i=2, k=1, ii=4, jj=0)
        dp[4][0][1][1], path[4][0][1][1] = g_fuc(i=2, k=1, ii=4, jj=1)

    # DP transitions
    for i in range(4, N + 1):

        dp[i - 1][2][0][0], path[i - 1][2][0][0] = f_fuc(i=i, k=0, ii=i - 1, jj=0)
        dp[i - 1][2][1][0], path[i - 1][2][1][0] = f_fuc(i=i, k=0, ii=i - 1, jj=1)

        dp[i][1][0][0], path[i][1][0][0] = g_fuc(i=i - 2, k=2, ii=i, jj=0)
        dp[i][1][1][0], path[i][1][1][0] = g_fuc(i=i - 2, k=2, ii=i, jj=1)
        dp[i][1][0][1], path[i][1][0][1] = f_fuc(i=i - 1, k=1, ii=i, jj=0)
        dp[i][1][1][1], path[i][1][1][1] = f_fuc(i=i - 1, k=1, ii=i, jj=1)

        if i < N:
            dp[i + 1][0][0][0], path[i + 1][0][0][0] = f_fuc(i=i - 2, k=2, ii=i + 1, jj=0)
            dp[i + 1][0][0][1], path[i + 1][0][0][1] = f_fuc(i=i - 1, k=1, ii=i + 1, jj=0)
            dp[i + 1][0][1][0], path[i + 1][0][1][0] = f_fuc(i=i - 2, k=2, ii=i + 1, jj=1)
            dp[i + 1][0][1][1], path[i + 1][0][1][1] = f_fuc(i=i - 1, k=1, ii=i + 1, jj=1)

    v1, v2, v3, v4, v5, v6 = dp[N][1][0][0], dp[N][1][1][0], dp[N][1][0][1], dp[N][1][1][1], \
                             dp[N - 1][2][0][0], dp[N - 1][2][1][0]

    if max(min(v1, v2), min(v3, v4)) > min(v5, v6):
        if min(v1, v2) > min(v3, v4):
            path[-1][-1][-1] = [N, 1, 0, 0] if v1 < v2 else [N, 1, 1, 0]
        else:
            path[-1][-1][-1] = [N, 1, 0, 1] if v3 < v4 else [N, 1, 1, 1]
    else:
        path[-1][-1][-1] = [N - 1, 2, 0, 0] if v5 < v6 else [N - 1, 2, 1, 0]
    max_value = max(max(min(v1, v2), min(v3, v4)), min(v5, v6))

    last, worst_path = path[-1][-1][-1], [[path[-1][-1][-1][0], 'A' if path[-1][-1][-1][2] == 0 else 'B']]
    worst_path_r = [path[-1][-1][-1].__str__()]
    for i in range(N - 1):
        worst_path_r.append(path[last[0]][last[1]][last[2]][last[3]].__str__())
        last = path[last[0]][last[1]][last[2]][last[3]]
        worst_path.append([last[0], 'A' if last[2] == 0 else 'B'])
    worst_path.reverse()
    worst_path_r.reverse()
    return max_value, worst_path, dp, worst_path_r


def sub_problem_single_T(N, A, B, pt, init_pos, st_line, touch_flag):
    dp = [[0] * 2 for _ in range(N + 1)]
    path = [[0] * 2 for _ in range(N + 1)]  # 用于记录路径选择

    def cal_dp_state(l, ii, ll):
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
    for i in range(1, N):
        # max(dp[i - 1][0] + abs(A[i - 1] - B[i]), st_line[i - 1]) + abs(B[i] - A[i]) + pt[i]
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

    min_cost = min(max(dp[N - 1][0], st_line[N - 1]), max(dp[N - 1][1], st_line[N - 1]))
    last_choice = 0 if dp[N - 1][0] < dp[N - 1][1] else 1

    optimal_path, schedule = [], []
    for i in range(N - 1, -1, -1):
        if last_choice == 0:
            optimal_path.append(('A', i))
            optimal_path.append(('B', i))
        else:
            optimal_path.append(('B', i))
            optimal_path.append(('A', i))
        last_choice = path[i][last_choice]
        schedule.append(dp[i][last_choice])

    optimal_path.reverse()  # 反转路径
    schedule.reverse()

    return schedule, min_cost


def find_new_sequences(J, original_sequence, master_results):
    # Step 1: Calculate original distances
    original_distances = [original_sequence[i + 1] - original_sequence[i] for i in range(len(original_sequence) - 1)]

    # Step 2: Find all combinations of new sequences of the same length as the original sequence
    candidates = combinations(J, len(original_sequence))

    # Step 3: Filter combinations based on distance condition
    valid_sequences = []
    for candidate in candidates:
        sorted_candidate = sorted(candidate)
        new_distances = [sorted_candidate[i + 1] - sorted_candidate[i] for i in range(len(sorted_candidate) - 1)]
        if all(new_distances[i] >= original_distances[i] for i in range(len(original_distances))) and list(
                candidate) != original_sequence:
            # if sorted_candidate not in invalid_sequences_r:
            flag = True
            for item in master_results.items():
                # print(J, original_sequence, master_results)
                if item[1] in original_sequence:
                    if item[0] in data.U_F and sorted_candidate[original_sequence.index(item[1])] not in data.I:
                        flag = False
                        break
            if not flag:
                continue
            valid_sequences.append(sorted_candidate)
            invalid_sequences_r.append(sorted_candidate)
    return valid_sequences


def generate_cuts(master_results, sub_results, added_cuts):
    if sub_results is not None:
        if data.K_num==1:
            X, objVal, _, _ = master_results
            max_value, worst_path, dp, worst_path_r = sub_results
            # optimal_cuts todo: 没有区分场桥 optimal delta inactivate
            added_cuts.append(['Optimal', max_value, [[key, X[key]] for key in X.keys() if key != data.U_num + 1]])
            # if data.K_num > 1:
            #     added_cuts.append(['Optimal delta', max_value,
            #                        [[key, X[key], data.U_num_set[key] * cf.unit_process_time] for key in X.keys()]])
            # translation cut todo: add多场桥cut
            U_k = [[] for _ in range(data.K_num)]
            for pair in X.items():
                U_k[data.J_K_dict[pair[1] + 1]].append([pair[0], pair[1]])
                # if isinstance(pair[1], int):
                #     U_k[data.J_K_dict[pair[1] + 1]].append([pair[0], pair[1]])
                # else:
                #     for pos in pair[1]:
                #         if pos + 1 in data.J_K_first:
                #             continue
                #         U_k[data.J_K_dict[pos + 1]].append([pair[0], pos])
            U_k_sorted = [sorted(sublist, key=lambda x: x[1]) for sublist in U_k]  # 按照贝位排序
            Used_bays = [[pair[1] for pair in U_k_sorted[k]] for k in range(data.K_num)]
            invalid_sequences = [-1 for _ in range(data.K_num)]
            for k in range(data.K_num):
                init_b_index = data.J_K[k].index(Used_bays[k][0] + 1)
                invalid_sequences[k] = find_new_sequences(data.J_K[k][init_b_index:], Used_bays[k], master_results[0])
            if not all(num == -1 for num in invalid_sequences):
                added_cuts.append(['Translation', max_value, U_k_sorted, invalid_sequences])
            a = 1
        else:
            a=1
    else:
        # todo:generate infeasible results
        a = 1
        print("infeasible")

    return added_cuts


if __name__ == '__main__':
    case = 'case1m'

    data = read_data('/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/standard/' + case)
    print(case)
    # ============== 生成所有可能提取序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    pi_num_set = [i for i in range(len(valid_permutations))]
    invalid_sequences_r = []
    # ============== 求解 ================
    random_seed = 0.2  # 加seq cut设置系数
    CCG()
