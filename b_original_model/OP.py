# -*- coding: utf-8 -*-
# @Time    : 2024/12/19 4:32 PM
# @Author  : JacQ
# @File    : OP.py
import time
from itertools import permutations

from gurobipy import *
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


def original_problem(data):
    """
        :param data: 数据
        :return:
    """
    # ============== 生成所有可能序列 ================
    sequence = list(range(data.G_num))
    valid_permutations = generate_permutations(sequence, swapped=None)
    pi_num = len(valid_permutations)

    # ============== 构造模型 ================
    big_M = 10000
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
    # Con1
    model.addConstrs((C[w][u] <= C_max[w] for w in range(pi_num) for u in data.U), "1b")
    # 对于一个子箱组
    # Con2: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "1c")
    model.addConstrs((quicksum(X[u][j - 1] for j in data.I) == 1 for u in data.U_F_u), "1d")
    model.addConstrs((quicksum(X[u][j + 1] for j in data.I) == 1 for u in data.U_F_l), "1e")
    # con2: Initial position restrictions
    model.addConstrs((X[data.U_num][j - 1] == 1 for j in data.J_K_first), "1f")
    # con3:对于40ft的子箱组占了前一个后一个位置就要被虚拟子箱组占用
    model.addConstrs((X[u][j - 1] <= big_M * X[u + 1][j + 1] for u in data.U_F_u for j in data.I), "1g")
    # Con4: 一个贝上放置的箱组一般不超过2个
    model.addConstrs((quicksum(X[u][j - 1] for u in range(data.U_num)) <= 2 for j in data.J), "1h")
    # con5: 20和40的不能放在一个贝
    model.addConstrs((X[u][j - 1] + X[uu][j - 1] <= 1 for j in data.J for u in data.U_L for uu in data.U_F), "1i")
    # Con6: 放置集装箱数不超过贝位的最大容量
    model.addConstrs((quicksum(data.U_num_set[u] * X[u][j - 1] for u in data.U) <= data.S_num * data.T_num
                      for j in range(data.J_num)), "1j")
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
    # con10:消除对称性
    model.addConstrs((X[uu][jj - 1] <= big_M * (1 - X[u][j - 1]) for u in data.U for uu in data.U for j in data.J
                      for jj in data.J if data.U_g_set[u] == data.U_g_set[uu] and jj < j and u < uu), "1q")
    # Con10: 优先级完成时间约束 fixme
    model.addConstrs((C[w][u] <= C[w][uu] for w in range(pi_num) for u in data.U for uu in data.U
                      if valid_permutations[w].index(data.U_g_set[u])
                      < valid_permutations[w].index(data.U_g_set[uu])), "1r")

    # ============== 构造目标 ================
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')  # 线性化模型变量
    model.addConstrs((obj >= C_max[w] for w in range(pi_num)), "obj")
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)

    # ============== 求解参数 ================
    # model.Params.OutputFlag = 0
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
        model.write('a.lp')
        model.write('a.sol')
        master_X = {}
        for u in range(data.U_num):
            for j in range(60 * data.K_num):
                if j + 1 not in data.J:
                    continue
                if abs(X[u][j].X) > 0.00001:
                    master_X[u] = j
        return master_X, obj.X


if __name__ == '__main__':
    dataa = read_data('/Users/jacq/PycharmProjects/BayAllocation/a_data_process/data/case1')
    master_X, obj = original_problem(dataa)
    print(obj)
