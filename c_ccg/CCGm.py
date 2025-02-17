# -*- coding: utf-8 -*-
# @Time    : 2024/12/18 3:47 PM
# @Author  : JacQ
# @File    : CCG1.py
import random
import time
from itertools import combinations, permutations
from matplotlib import pyplot as plt, patches
import gurobipy as gp
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
        if CCG_show_flag: print("=========================================iteration\t" +
                                str(iteration_num) + "=========================================")

        master_results = update_master_problem_adding_cuts(master_results[2], master_results[3], added_cuts)

        if not master_results:
            if CCG_show_flag: print("Master problem is infeasible!")
            return
        LB = max(LB, master_results[1])  # 当前迭代主问题目标值
        if CCG_show_flag: print("Master obj:\t", master_results[1], "\t", master_results[0])

        # 求解子问题
        if master_results[0] == {0: 14, 1: 10, 2: 6, 3: 16, 4: 10, 5: 6}:
            a = 1
        if data.K_num == 1:
            sub_results = sub_problem_single(data, master_results[0])
        else:
            sub_results = sub_problem_multi(data, master_results[0])
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
            if CCG_show_flag: print("UB:\t", UB)
            if CCG_show_flag: print("LB:\t", LB)
            if CCG_show_flag: print("Gap:\t", round((UB - LB) * 100 / UB, 1))
            if CCG_show_flag: print("Time:\t", end_time - start_time)
            if end_time - start_time >= 7200:
                print("=========================================Final results=========================================")
                print("Final UB:\t", UB)
                print("Final LB:\t", LB)
                print("Gap:\t", round((UB - LB) * 100 / UB, 1))
                print("Time\t", time.time() - start_time)
                break
        iteration_num += 1

    return UB, LB, round((UB - LB) * 100 / UB, 1), time.time() - start_time, end_x


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
        elif cut[0] == 'Optimal_Delta':
            v, XX = cut[1], cut[2]
            rhs = LinExpr(0)
            pt_sum = 0
            for u, j, pt in XX:
                pt_sum += pt
                rhs.addTerms(pt, X[u][j])
            model.addConstr(theta >= v + rhs - pt_sum, "opt cut")
        elif cut[0] == 'Translation':
            a = 1
            # containers = [[term[0] for term in cut[2][k]] for k in range(data.K_num)]
            # zeta = [[0 for _ in range(len(cut[3][k]) + 1)] for k in range(data.K_num)]
            # for k in range(data.K_num):
            #     for q in range(len(cut[3][k]) + 1):
            #         if q == len(cut[3][k]):
            #             compose = [tmp[1] + 1 for tmp in cut[2][k]]
            #         else:
            #             compose = cut[3][k][q]
            #         rhs2 = LinExpr(0)
            #         for cnt in range(len(compose)):
            #             rhs2.addTerms(1, X[containers[k][cnt]][compose[cnt] - 1])
            #         zeta[k][q] = model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name='zeta_' + str(zeta_cnt))
            #         model.addConstr(zeta[k][q] >= rhs2 - len(compose) + 1, 'zeta_' + str(zeta_cnt))
            #         zeta_cnt += 1
            # rhs1, rhs2 = LinExpr(0), LinExpr(0)
            # for k in range(data.K_num):
            #     rhs1.addTerms(1, zeta[k][len(cut[3][k])])
            # eta2 = model.addVar(0, GRB.INFINITY, vtype=GRB.CONTINUOUS, name='eta2')
            # model.addConstr(eta2 >= big_M * (rhs1 - data.K_num + 1), 'eta2')
            # for k in range(data.K_num):
            #     for q in range(len(zeta[k])):
            #         rhs2.addTerms(1, zeta[k][q])
            # model.addConstr(theta + eta2 >= big_M * (rhs2 - data.K_num + 1), 'eta2')

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


def sub_problem_help(data, master_X, pi_index=0):
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
            G_u_pos[data.U_g_set[u]].append(j+1)
        pos = [[min(G_u_pos[g]) * cf.unit_move_time if g not in data.U_F else (min(G_u_pos[g]) - 2) * cf.unit_move_time,
                max(G_u_pos[g]) * cf.unit_move_time] for g in range(data.G_num)]  # 每个箱组AB子箱组位置
        pt = [g_num * cf.unit_process_time for g_num in data.G_num_set]
        pt = [x for x in pt if x != 0]
    else:
        G_u_pos = [[[] for _ in range(data.G_num)] for _ in range(data.K_num)]  # 每个箱组子箱组位置
        for u in range(data.U_num):
            j = master_X[u]
            k_index = data.J_K_dict[j + 1]
            G_u_pos[k_index][data.U_g_set[u]].append(j + 1)
        pos = [[[] for g in range(data.G_num)] for k in range(data.K_num)]
        for k in range(data.K_num):
            for g in range(data.G_num):
                if len(G_u_pos[k][g]) == 0:
                    pos[k][g] = [None, None]
                elif g not in data.U_F:
                    pos[k][g] = [min(G_u_pos[k][g]) * cf.unit_move_time, max(G_u_pos[k][g]) * cf.unit_move_time]
                else:
                    pos[k][g] = [(min(G_u_pos[k][g]) - 2) * cf.unit_move_time, max(G_u_pos[k][g]) * cf.unit_move_time]
        pt = [[data.G_num_set[u] * cf.unit_process_time if data.J_K_dict[master_X[u] + 1] == k else 0
               for u in range(data.U_num)] for k in range(data.K_num)]  #
        pi_ls = valid_permutations[pi_index]
        # 根据pi调整ls顺序
        pos = [[pos[k][index] for index in pi_ls] for k in range(data.K_num)]
        pt = [[pt[k][index] for index in pi_ls] for k in range(data.K_num)]
        for k in range(data.K_num):
            for i in range(data.G_num):
                if pt[k][i] == 0:
                    if i == 0:
                        pos[k][i] = [0, 0]
                    else:
                        pos[k][i] = pos[k][i - 1]
    return data.G_num, pos, pt, [(data.J_K_first[k]) * cf.unit_move_time for k in range(data.K_num)]


def sub_problem_multi(data, master_result):
    # todo:顺序
    max_v, max_pi, max_schedule = 0, -1, []
    for pi in range(pi_num):
        N, pos, pt, init_pos = sub_problem_help(data, master_result, pi)
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
                                                          init_pos=init_pos[k], st_line=st_line,
                                                          touch_flag=touch_flag[k])
                if whole_schedule[k] == schedule:
                    stop_flag[k] = True
                else:
                    whole_schedule[k] = schedule  # 创建一个计数器
            # terminate when no adjustment
            if all(stop_flag): break
            # adjustment
            st_line, touch_flag = [], [[False for _ in range(data.G_num)] for _ in range(data.K_num)]
            st = [[whole_schedule[k][g + 1] - whole_schedule[k][g] - pt[k][g + 1] for g in range(data.G_num - 1)] for k
                  in
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
        if st_line[-1] > max_v:
            max_v = st_line[-1]
            max_pi = pi
            max_schedule = whole_schedule
    return max_v, max_pi, max_schedule


def sub_problem_single(data, master_result):
    """
        :param N: 几个箱组
        :param pos: 对应A，B子箱组位置
        :param pt: 每个箱组所需加工时间
        :return: (1/0表示是否可行, 目标值)
             不可行 (0，float("inf"))
             可行 (1，目标值)
    """
    N, pos, pt, init_pos = sub_problem_help(data, master_result)

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
            try:
                return abs(pre - init_pos) + abs(pre - suc) + pt[0]
            except:
                a = 2
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
    if dp[N - 1][0] < dp[N - 1][1]:
        last_choice = 0
    else:
        last_choice = 1

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
        if data.K_num == 1:
            X, objVal, _, _ = master_results
            max_value, worst_path, dp, worst_path_r = sub_results
            # optimal_cuts
            added_cuts.append(['Optimal', max_value, [[key, X[key]] for key in X.keys() if key != data.U_num + 1]])

            # translation cut
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
            for k in range(data.K_num):
                init_b_index = data.J_K[k].index(Used_bays[k][0] + 1)
                invalid_sequences = find_new_sequences(data.J_K[k][init_b_index:], Used_bays[k], master_results[0])
                if invalid_sequences:
                    added_cuts.append(['Translation', max_value, [[U_k_sorted[k][i][0], seq[i]]
                                                                  for seq in invalid_sequences for i in
                                                                  range(len(seq))]])
            a = 1
        else:
            X, objVal, _, _ = master_results
            v, pi, schedule = sub_results
            # optimal_cuts
            added_cuts.append(['Optimal', v, [[key, X[key]] for key in X.keys() if key != data.U_num + 1]])
            added_cuts.append(['Optimal_Delta', v,
                               [[key, X[key], data.U_num_set[key] * cf.unit_process_time] for key in X.keys()]])
            # translation cut
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
                added_cuts.append(['Translation', v, U_k_sorted, invalid_sequences])


    else:
        # todo:generate infeasible results
        a = 1
        print("infeasible")

    return added_cuts


def original_problem_robust(data, res=None, w_obj=None):
    """
        :param data: 数据
        :return:
    """
    s_t = time.time()

    # ============== 构造模型 ================
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
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U), "1c")
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
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu]
                      and data.U_num_set[u] == data.U_num_set[uu]
                      and jj < j and u < uu), "1q")
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
        return obj.X, res, worst_index, time.time() - s_t


def generate_color_groups(g):
    cmap = plt.cm.get_cmap('tab10', g)  # 从 'tab10' 调色板中生成 g 种颜色
    colors = [cmap(i) for i in range(g)]
    return colors


if __name__ == '__main__':
    case = 'case1m'
    print_flag, CCG_show_flag = True, False
    # data = read_data('C:\\Users\\admin\\PycharmProjects\\BayAllocation\\a_data_process\\data\\standard\\' + case)
    data = read_data('/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/standard/' + case)
    print(case)
    # ============== 生成所有可能提取序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    pi_num = len(valid_permutations)
    pi_num_set = [i for i in range(len(valid_permutations))]
    invalid_sequences_r = []
    # ============== 求解 ================
    random_seed = 0.2  # 加seq cut设置系数
    res = original_problem_robust(data)
    if res is False:
        obj1, time1, res1 = -1, 7200, []
    else:
        obj1, res1, worst_index1, time1 = res

    ub, lb, gap, time2, res2 = CCG()
    print("================================================================================")
    print(f"{case}\tOP:obj\t{obj1}\ttime\t{time1:.2f}\t"
          f"CCG:ub\t{ub:.2f}\tlb\t{lb:.2f}\tgap\t{gap:.2f}\ttime\t{time2:.2f}")
    print(res1)
    print(res2)

    # with open("C:\\Users\\admin\\PycharmProjects\\BayAllocation\\b_original_model\\output_m.txt", "a") as f:
    with open("/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/output_m.txt", "a") as f:
        # f.write("This is a test output.\n")
        # f.write("Second line of output.\n")
        f.write(f"{case}\tOP:obj\t{obj1}\ttime\t{time1:.2f}\t"
                f"CCG:ub\t{ub:.2f}\tlb\t{lb:.2f}\tgap\t{gap:.2f}\ttime\t{time2:.2f}\n")
