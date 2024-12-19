# -*- coding: utf-8 -*-
# @Time    : 2024/12/17 8:08 PM
# @Author  : JacQ
# @File    : LBBD.py
import itertools
import random
import sys
import time
from math import ceil
from itertools import combinations

import gurobipy
from gurobipy import GRB, LinExpr, quicksum

import numpy as np


def sub_problem(b, D, U, U_bound, L, data, scenario):
    '''
    :return: (1/0表示是否可行, 目标值, 每个场桥在该场景下的工作时间下界 * 场景系数)
             不可行 (0，float("inf")，[])
             可行 (1/0表示是否可行，目标值，每个场桥在该场景下的工作时间下界 * 场景系数)
    '''
    # ============= 场景参数 =============
    scenarioNum = 1
    coff_scenario = data.coff_scenarios[scenario]
    priority = data.priority[scenario]
    model = gurobipy.Model("SubProblem")
    model.Params.OutputFlag = 0

    # ============= 决策变量定义 =============
    # 定义字典用来存放决策变量
    Y = {}
    C = {}
    C_d = {}
    O = {}
    O_d = {}
    C_max = {}

    # Y
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                for u1 in U_bound[k][d]:
                    for u2 in U_bound[k][d]:
                        name = 'Y_' + str(w) + '_' + str(k) + '_' + str(d) + '_' + str(u1) + '_' + str(u2)
                        Y[w, k, d, u1, u2] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)
    # C_kd
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                name = 'C_' + str(w) + '_' + str(k) + '_' + str(d)
                C[w, k, d] = model.addVar(0, 10000000, vtype=GRB.CONTINUOUS, name=name)
    # C_d
    for w in range(scenarioNum):
        for d in range(len(data.destination)):
            C_d[w, d] = model.addVar(0, 10000000, vtype=GRB.CONTINUOUS, name='C_' + str(w) + '_' + str(d))
    # C_max
    for w in range(scenarioNum):
        name = 'C_max' + '_' + str(w)
        C_max[w] = model.addVar(0, 1000000000, vtype=GRB.CONTINUOUS, name=name)

    # O_kd
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                name = 'O_' + str(w) + '_' + str(k) + '_' + str(d)
                O[w, k, d] = model.addVar(0, 10000000, vtype=GRB.CONTINUOUS, name=name)
    # O_d
    for w in range(scenarioNum):
        for d in range(len(data.destination)):
            O_d[w, d] = model.addVar(0, 10000000, vtype=GRB.CONTINUOUS, name='O_' + str(w) + '_' + str(d))

    # ============= 目标函数 =============
    obj = LinExpr(0)
    for i in range(scenarioNum):
        obj.addTerms(coff_scenario, C_max[i])
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)

    # ============= 约束 =============
    # 每个场桥区域内的每个卸货港子箱组只有一个起始子箱组和一个末尾子箱组
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                model.addConstr(quicksum(Y[w, k, d, u1, u2] for u1 in U_bound[k][d] for u2 in U_bound[k][d]) == 1)
    # 每个卸货港子箱组的最早开始提取时间和最晚结束提取时间与每个场桥上对应时间点的关系
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                model.addConstr(C[w, k, d] <= C_d[w, d])
    for w in range(scenarioNum):
        for d in range(len(data.destination)):
            model.addConstr(C_d[w, d] <= C_max[w])
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                model.addConstr(O[w, k, d] >= O_d[w, d])
    # 每个卸货港子箱组提取的时间跨度不超过𝑡_𝛿。
    for w in range(scenarioNum):
        for d in range(len(data.destination)):
            model.addConstr(C_d[w, d] - O_d[w, d] <= data.delta_time)
    # 前后提取的卸货港子箱组结束提取时间关系
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d1 in D[k]:
                for d2 in D[k]:
                    if d1 != d2:
                        for u1 in U_bound[k][d1]:
                            for u2 in U_bound[k][d2]:
                                for u3 in U_bound[k][d2]:
                                    model.addConstr(C[w, k, d2] - C[w, k, d1] + data.bigM * (
                                            4 - priority[d1][d2] - quicksum(
                                        Y[w, k, d1, u, u1] for u in U_bound[k][d1]) - quicksum(
                                        Y[w, k, d2, u2, u] for u in U_bound[k][d2]) - quicksum(
                                        Y[w, k, d2, u, u3] for u in U_bound[k][d2]))
                                                    >= data.t_move * data.disMatrix[b[u1]][b[u2]] + data.t_move * (
                                                            2 * L[k, d2] - data.disMatrix[b[u2]][b[u3]])
                                                    + quicksum((data.t_relocation * data.gamma[u] + data.t_pickup *
                                                                data.containerNum_of_subGroup[u]) for u in U[k][d2] if
                                                               u in data.subGroup1)
                                                    + quicksum((2 * data.t_relocation * data.gamma[
                                        u] + 2 * data.t_pickup * data.containerNum_of_subGroup[u]) for u in U[k][d2] if
                                                               u in data.subGroup2))
    # 每个场桥上首个提取的卸货港子箱组结束时间
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d1 in D[k]:
                if sum([priority[d1][d2] for d2 in D[k]]) == len(D[k]) - 1:
                    for u1 in U_bound[k][d1]:
                        model.addConstr(
                            O[w, k, d1] + data.bigM * (1 - quicksum(Y[w, k, d1, u1, u2] for u2 in U_bound[k][d1]))
                            >= data.t_move * data.disMatrix[data.bayNum_of_Block[k][0]][b[u1]])
    # 场桥𝑘上每个卸货港提取的开始和结束时间之差至少为提取该港子箱组的总翻倒、提取时间以及移动时间之和
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d in D[k]:
                for u1 in U_bound[k][d]:
                    for u2 in U_bound[k][d]:
                        model.addConstr(
                            C[w, k, d] - data.t_move * (2 * L[k, d] - data.disMatrix[b[u1]][b[u2]]) + data.bigM * (
                                    1 - Y[w, k, d, u1, u2])
                            - quicksum(
                                (data.t_relocation * data.gamma[u] + data.t_pickup * data.containerNum_of_subGroup[u])
                                for u in U[k][d] if u in data.subGroup1)
                            - quicksum((2 * data.t_relocation * data.gamma[u] + 2 * data.t_pickup *
                                        data.containerNum_of_subGroup[u]) for u in U[k][d] if u in data.subGroup2)
                            >= O[w, k, d])
    # 若卸货港𝑑在卸货港𝑑′前提取，那么卸货港𝑑子箱组的最晚结束提取时间不超过每个场桥上卸货港𝑑′的起始子箱组结束时间
    for w in range(scenarioNum):
        for k in range(data.craneNum):
            for d1 in range(len(data.destination)):
                for d2 in D[k]:
                    for u1 in U_bound[k][d2]:
                        if (d1 != d2) and (priority[d1][d2] == 1):
                            if u1 in data.subGroup1:
                                model.addConstr(C_d[w, d1] - data.bigM * (
                                        1 - quicksum(Y[w, k, d2, u1, u2] for u2 in U_bound[k][d2])) <=
                                                O[w, k, d2] + (data.t_relocation * data.gamma[u1] + data.t_pickup *
                                                               data.containerNum_of_subGroup[u1]))
                            if u1 in data.subGroup2:
                                model.addConstr(C_d[w, d1] - data.bigM * (
                                        1 - quicksum(Y[w, k, d2, u1, u2] for u2 in U_bound[k][d2])) <=
                                                O[w, k, d2] + (
                                                        2 * data.t_relocation * data.gamma[u1] + 2 * data.t_pickup *
                                                        data.containerNum_of_subGroup[u1]))
    # ============= 模型求解 & 结果打印 =============
    model.setParam(GRB.Param.InfUnbdInfo, 1)
    model.optimize()

    # ============= 返回目标值 & 是否可行 =============
    if model.status == GRB.INFEASIBLE:
        return 0, float("inf"), []
    else:
        print("w=", scenario)
        print(C_d)
        print(C)
        print(model.objVal)
        K = []
        for k in range(data.craneNum):
            temp = [C[key].x for key in C.keys() if (key[0] == 0) and (key[1] == k)]
            if not temp:
                K.append(0)
                continue
            cum_time = 0
            for d in D[k]:
                for u in U[k][d]:
                    if u in data.subGroup1:
                        cum_time += data.t_relocation * data.gamma[u] + data.t_pickup * data.containerNum_of_subGroup[u]
                    if u in data.subGroup2:
                        cum_time += 2 * data.t_relocation * data.gamma[u] + 2 * data.t_pickup * \
                                    data.containerNum_of_subGroup[u]
            for d in D[k]:
                temp = 0
                for u1 in U_bound[k][d]:
                    for u2 in U_bound[k][d]:
                        temp = max(data.disMatrix[b[u1]][b[u2]], temp)
                cum_time += temp * data.t_move
            for d1 in D[k]:
                for d2 in D[k]:
                    if (d1 != d2) and (sum([priority[d1][d] for d in D[k]]) - sum(
                            [priority[d2][d] for d in D[k]]) == 2):  # d1在d2前面
                        temp = 100000000
                        for u1 in U_bound[k][d1]:
                            for u2 in U_bound[k][d2]:
                                temp = min(data.disMatrix[b[u1]][b[u2]], temp)
                        cum_time += temp * data.t_move
            for d in D[k]:
                if sum([priority[d][d2] for d2 in D[k]]) == len(D[k]) - 1:
                    temp = 100000000
                    for u1 in U_bound[k][d]:
                        temp = min(data.disMatrix[data.bayNum_of_Block[k][0]][b[u1]], temp)
                    cum_time += temp * data.t_move
                    for d1 in D[k]:
                        if (d1 != d) and (
                                sum([priority[d][d2] for d2 in D[k]]) - sum([priority[d1][d2] for d2 in D[k]]) == 2):
                            if max([b[u] for u in U[k][d1]]) <= min([b[u] for u in U[k][d]]):
                                cum_time += temp * data.t_move
            for d1 in D[k]:
                for d2 in D[k]:
                    for d3 in D[k]:
                        if (d1 != d2) and (d2 != d3) and (d1 != d3) and sum([priority[d1][d] for d in D[k]]) - sum(
                                [priority[d2][d] for d in D[k]]) == 2 \
                                and sum([priority[d2][d] for d in D[k]]) - sum([priority[d3][d] for d in D[k]]) == 2:
                            if (max([b[u] for u in U[k][d1]]) - min([b[u] for u in U[k][d2]]) <= 0 and max(
                                    [b[u] for u in U[k][d3]]) - min([b[u] for u in U[k][d2]]) <= 0) \
                                    | (min([b[u] for u in U[k][d1]]) - max([b[u] for u in U[k][d2]]) >= 0 and min(
                                [b[u] for u in U[k][d3]]) - max([b[u] for u in U[k][d2]]) >= 0):
                                temp = 0
                                for u1 in U_bound[k][d2]:
                                    for u2 in U_bound[k][d2]:
                                        temp = max(data.disMatrix[b[u1]][b[u2]], temp)
                                cum_time += temp * data.t_move
            K.append(cum_time * coff_scenario)
        return 1, model.objVal, K


def master_problem(data, added_cuts, LB):
    '''
    :param data: 数据
    :param added_cuts: 所有cut的集合
    :param LB: 当前迭代的下界LB
    :return: 不可行，返回 FALSE；
             可行，返回 master_X, theta.x, model.objVal, master_obj；分别表示解、theta值、主问题目标值、主问题目标值-theta值
    '''
    model = gurobipy.Model("MasterProblem")
    X = {}
    for u in range(data.subGroupNum):
        for j in range(data.bayNum):
            name = 'X_' + str(u) + '_' + str(j)
            X[u, j] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, u放在j上则为1
    for k in range(data.craneNum):
        j = data.bayNum_of_Block[k][0]
        X[data.subGroupNum, j] = model.addVar(1, 1, vtype=GRB.BINARY,
                                              name='X_' + str(data.subGroupNum) + '_' + str(j))  # 0-1变量, u放在j上则为1
    theta = model.addVar(0, float("inf"), vtype=GRB.CONTINUOUS, name='theta')

    # ============= 目标函数 =============
    obj = LinExpr(0)
    obj.addTerms(1, theta)
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)
    model.Params.OutputFlag = 0

    # ============= 约束 =============
    # 一个贝位最多混贝两个子箱组
    for j in range(data.bayNum):
        lhs = LinExpr(0)
        for u in range(data.subGroupNum):
            lhs.addTerms(1, X[u, j])
        model.addConstr(lhs <= 2)
    # 20TEU和40TEU不能混贝
    for j in range(data.bayNum):
        for u1 in data.subGroup1:
            for u2 in data.subGroup2:
                model.addConstr(X[u1, j] + X[u2, j] <= 1)
    # 总箱数不能超过贝位限制
    for j in range(data.bayNum):
        lhs = LinExpr(0)
        for u in range(data.subGroupNum):
            lhs.addTerms(ceil(data.containerNum_of_subGroup[u] / data.t), X[u, j])
        model.addConstr(lhs <= data.s)
    # 一个子箱组只能放在一个贝位上
    for u in range(data.subGroupNum):
        lhs = LinExpr(0)
        for j in range(data.bayNum):
            lhs.addTerms(1, X[u, j])
        model.addConstr(lhs == 1)
    # 40TEU的两个子箱组必须放在相邻两个贝位都可以堆放的贝位上
    for u in data.subGroup2:
        lhs = LinExpr(0)
        for j in data.F:
            lhs.addTerms(1, X[u, j])
        model.addConstr(lhs == 1)
    # 40TEU所在贝位的后一贝位不放子箱组
    for j in data.F:
        for u in data.subGroup2:
            lhs = LinExpr(0)
            lhs.addTerms(2, X[u, j])
            for u_temp in range(data.subGroupNum):
                lhs.addTerms(1, X[u_temp, j + 1])
            model.addConstr(lhs - 2 <= 0)

    # ============= 初始有效不等式 =============
    # 消除对称性
    for l in range(len(data.subGroup_with_same_destination)):
        temp1 = [u for u in data.subGroup_with_same_destination[l] if u in data.subGroup1]
        if len(temp1) > 1:
            for u1 in temp1:
                for u2 in temp1:
                    if (u1 < u2) & (data.containerNum_of_subGroup[u1] == data.containerNum_of_subGroup[u2]):
                        for j2 in range(data.bayNum):
                            for j1 in range(j2 + 1, data.bayNum):
                                model.addConstr(X[u2, j2] - (1 - X[u1, j1]) <= 0)
        temp2 = [u for u in data.subGroup_with_same_destination[l] if u in data.subGroup2]
        if len(temp2) > 1:
            for u1 in temp2:
                for u2 in temp2:
                    if (u1 < u2) & (data.containerNum_of_subGroup[u1] == data.containerNum_of_subGroup[u2]):
                        for j2 in range(data.bayNum):
                            for j1 in range(j2 + 1, data.bayNum):
                                model.addConstr(X[u2, j2] - (1 - X[u1, j1]) <= 0)
    # 初始下界
    model.addConstr(data.craneNum * theta >= quicksum(
        (data.t_relocation * data.gamma[u] + data.t_pickup * data.containerNum_of_subGroup[u]) for u in
        list(range(data.subGroupNum)) if u in data.subGroup1)
                    + quicksum(
        (2 * data.t_relocation * data.gamma[u] + 2 * data.t_pickup * data.containerNum_of_subGroup[u]) for u in
        list(range(data.subGroupNum)) if u in data.subGroup2))
    # 可行性不等式
    L = {}
    for k in range(data.craneNum):
        for d in range(len(data.destination)):
            L[k, d] = model.addVar(0, 1000000000, vtype=GRB.CONTINUOUS, name='L_' + str(k) + '_' + str(d))
            for u1 in data.subGroup_with_same_destination[d]:
                for u2 in data.subGroup_with_same_destination[d]:
                    for j1 in data.bayNum_of_Block[k]:
                        for j2 in data.bayNum_of_Block[k]:
                            if u1 != u2: model.addConstr(L[k, d] >= data.t_move * data.disMatrix[j1][j2] - data.bigM * (
                                    2 - X[u1, j1] - X[u2, j2]))
            model.addConstr(L[k, d] + quicksum(data.t_pickup * data.containerNum_of_subGroup[u] * X[u, j] for u in
                                               data.subGroup_with_same_destination[d] if u in data.subGroup1 for j in
                                               data.bayNum_of_Block[k])
                            + quicksum(2 * data.t_pickup * data.containerNum_of_subGroup[u] * X[u, j]
                                       for u in data.subGroup_with_same_destination[d] if u in data.subGroup2
                                       for j in data.bayNum_of_Block[k])
                            + quicksum(data.t_relocation * data.gamma[u] * X[u, j]
                                       for u in data.subGroup_with_same_destination[d] if u in data.subGroup1
                                       for j in data.bayNum_of_Block[k])
                            + quicksum(2 * data.t_relocation * data.gamma[u] * X[u, j]
                                       for u in data.subGroup_with_same_destination[d] if u in data.subGroup2
                                       for j in data.bayNum_of_Block[k]) <= data.delta_time)

    # ============= 加cut =============
    for item in added_cuts:
        if len(item) == 2:
            if item[1] < LB:
                continue
            keys = [key for key in item[0].keys() if (item[0][key] == 1) and (key[0] < data.subGroupNum)]
            rhs = LinExpr(0)
            for key in keys:
                rhs.addTerms(item[1], X[key])
            model.addConstr(theta + item[1] * len(keys) - item[1] >= rhs)
        else:
            if item[2] < LB:
                continue
            if item[3] == data.craneNum:
                model.addConstr(theta >= item[2] *
                                (quicksum((quicksum(X[u, j] for u in item[0][k][d] for j in item[1][k][d]) - len(
                                    item[0][k][d])) for k in range(data.craneNum) for d in range(len(data.destination)))
                                 + quicksum((quicksum(X[u, j] for u in item[0][k][d]
                                                      for j in item[4][k, d][0]) - len(item[4][k, d][0]))
                                            for k in range(data.craneNum) for d in range(len(data.destination)))
                                 + quicksum((quicksum(
                                            X[u, j] for u in [u for u in item[0][k][d] if u in data.subGroup1] for j in
                                            item[4][k, d][1]) - len(item[4][k, d][1])) for k in range(data.craneNum) for
                                            d in range(len(data.destination)))
                                 + quicksum((quicksum(
                                            X[u, j] for u in [u for u in item[0][k][d] if u in data.subGroup2] for j in
                                            item[4][k, d][2]) - len(item[4][k, d][2])) for k in range(data.craneNum) for
                                            d in range(len(data.destination))) + 1))
            else:
                model.addConstr(theta >=
                                item[2] * (quicksum(quicksum(X[u, j] for u in item[0][item[3]][d]
                                                             for j in item[1][item[3]][d]) -
                                                    len(item[0][item[3]][d]) +
                                                    quicksum(X[u, j] for u in item[0][item[3]][d]
                                                             for j in item[4][item[3], d][0])
                                                    - len(item[4][item[3], d][0])
                                                    + quicksum(X[u, j] for u in item[0][item[3]][d]
                                                               if u in data.subGroup1
                                                               for j in item[4][item[3], d][1])
                                                    - len(item[4][item[3], d][1])
                                                    + quicksum(X[u, j] for u in item[0][item[3]][d]
                                                               if u in data.subGroup2
                                                               for j in item[4][item[3], d][2])
                                                    - len(item[4][item[3], d][2])
                                                    for d in range(len(data.destination))) + 1))

    # ============= 主问题的解 =============
    model.optimize()
    # print("\n\n-----optimal value-----")
    # print(model.ObjVal)
    # print(theta.VarName + ' = ', theta.x)
    master_X = {}
    master_obj = 0
    for u in range(data.subGroupNum):
        for j in range(data.bayNum):
            # print(X[u, j].VarName + ' = ', X[u, j].x)
            master_X[u, j] = 0 if (abs(X[u, j].x) < 0.00001) else 1
    if model.status == GRB.INFEASIBLE:
        return False
    else:
        return master_X, theta.x, model.objVal, master_obj


def added_cuts_permutation(master_X, added_cuts, sub_obj_sum, real_k):
    longest_row = np.max([len(row) for row in data.bayNum_of_Block])
    keys = [key for key in master_X.keys() if (master_X[key] == 1) & (key[0] < data.subGroupNum)]
    # 记录每种类型箱组的：位置
    table = {}
    for key in keys:
        u = key[0]
        for d in range(len(data.destination)):
            if u in data.subGroup_with_same_destination[d]:
                if u in data.subGroup1:
                    if (d, '20', data.containerNum_of_subGroup[u]) in table.keys():
                        table[d, '20', data.containerNum_of_subGroup[u]].append(key[1])
                    else:
                        table[d, '20', data.containerNum_of_subGroup[u]] = [key[1]]
                else:
                    if (d, '40', data.containerNum_of_subGroup[u]) in table.keys():
                        table[d, '40', data.containerNum_of_subGroup[u]].append(key[1])
                    else:
                        table[d, '40', data.containerNum_of_subGroup[u]] = [key[1]]
                break
            else:
                continue

    types = list(table.keys())

    def generate_combinations(master_X, table, current_combination, result, index):
        if index == len(table):
            result.append(current_combination.copy())
            return

        key = types[index]
        candidates = [u for u in data.subGroup_with_same_destination[key[0]]
                      if (u in data.subGroup1 if key[1] == '20' else u in data.subGroup2)
                      and (data.containerNum_of_subGroup[u] == key[2])]
        if len(candidates) - len(table[key]) + 1 <= 0:
            generate_combinations(master_X, table, current_combination, result, index + 1)
        for i in range(len(candidates) - len(table[key]) + 1):  # 0,1,2
            for j in range(len(table[key])):
                current_combination[candidates[i + j], table[key][j]] = 1
            generate_combinations(master_X, table, current_combination, result, index + 1)
            for j in range(len(table[key])):
                current_combination[candidates[i + j], table[key][j]] = 0

    def replace_combinations(master_X, table):
        result = []
        current_combination = {}
        generate_combinations(master_X, table, current_combination, result, 0)
        return result

    result = replace_combinations(master_X, table)
    for i in range(len(result)):
        keys = [key for key in result[i].keys() if (result[i][key] == 1) & (key[0] < data.subGroupNum)]
        X = {}
        for key in keys:
            X[key] = 1
        U = U_dict(data, X)
        added_cuts = inner_change(X, added_cuts, sub_obj_sum, U, real_k)
        # added_cuts.append([X, sub_obj_sum])
    return added_cuts


def inner_change(master_X, added_cuts, sub_obj_sum, U, i):
    """
    同港互换割：在added_cuts集合中添加同港互换割所需集合
    :param master_X: 完整/部分 贝位分配方案
    :param added_cuts: cut集合
    :param sub_obj_sum: 目标值
    :param U: 各场桥上各卸货港的子箱组集合
    :param i: 完整贝位分配方案 —— 总场桥数；部分贝位分配方案 —— 场桥标号
    """
    if i < data.craneNum:
        J = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
        b = {}
        bound_j = {}
        for key in master_X.keys():
            b[key[0]] = key[1]
        for d in range(len(data.destination)):
            bound_j[i, d] = [[], [], []]
            for u in U[i][d]:
                if u in data.subGroup1 and b[u] not in J[i][d]:
                    J[i][d].append(b[u])
                elif u in data.subGroup2:
                    if b[u] not in J[i][d]:
                        J[i][d].append(b[u])
                    if b[u] + 1 not in J[i][d]:
                        J[i][d].append(b[u] + 1)
            bound_j[i, d][0] = [j for j in J[i][d] if j - 1 not in J[i][d]]
            temp = [j for j in J[i][d] if j + 1 not in J[i][d]]
            for j in temp:
                if j not in [b[u] for u in U[i][d]]:
                    bound_j[i, d][2].append(j - 1)
                else:
                    bound_j[i, d][1].append(j)
        if [U, J, sub_obj_sum, i, bound_j] not in added_cuts:
            added_cuts.append([U, J, sub_obj_sum, i, bound_j])
        # bound_j[i, d][0]表示场桥i上卸货港d的子箱组所占用的连续贝位的端点（包括40TEU子箱组的后一贝位）
        # bound_j[i, d][1]表示场桥i上卸货港d的子箱组连续段端点中，20英尺箱所占的贝位集合
        # bound_j[i, d][2]表示场桥i上卸货港d的子箱组连续段端点中，40英尺箱所占的贝位集合
    else:
        J = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
        b = {}
        bound_j = {}
        keys = [key for key in master_X.keys() if (master_X[key] == 1) & (key[0] < data.subGroupNum)]
        for key in keys:
            b[key[0]] = key[1]
        for k in range(data.craneNum):
            for d in range(len(data.destination)):
                bound_j[k, d] = [[], [], []]
                for u in U[k][d]:
                    if u in data.subGroup1 and b[u] not in J[k][d]:
                        J[k][d].append(b[u])
                    elif u in data.subGroup2:
                        if b[u] not in J[k][d]:
                            J[k][d].append(b[u])
                        if b[u] + 1 not in J[k][d]:
                            J[k][d].append(b[u] + 1)
                bound_j[k, d][0] = [j for j in J[k][d] if j - 1 not in J[k][d]]
                temp = [j for j in J[k][d] if j + 1 not in J[k][d]]
                for j in temp:
                    if j not in [b[u] for u in U[k][d]]:
                        bound_j[k, d][2].append(j - 1)
                    else:
                        bound_j[k, d][1].append(j)
        if [U, J, sub_obj_sum, data.craneNum, bound_j] not in added_cuts:
            added_cuts.append([U, J, sub_obj_sum, data.craneNum, bound_j])
    return added_cuts


def sub_help(master_X, data):
    '''
    :return: b 每个子箱组的贝位号,
             D 每个场桥上箱子的卸货港集合,
             U 记录每个场桥上每个卸货港的子箱组,
             U_bound 记录每个场桥上每个卸货港的连续贝位段的端点子箱组,
             L 计算每个场桥上每个卸货港的子箱组最远距离
    '''
    b = [-1] * data.subGroupNum
    keys = [key for key in master_X.keys() if (master_X[key] == 1) & (key[0] < data.subGroupNum)]
    for key in keys:
        b[key[0]] = key[1]
    # D_k
    D = []
    for k in range(data.craneNum):
        tmp = []
        for d in range(len(data.destination)):
            for u in data.destination[d]:
                for sub_u in data.subGroup_dict[u]:
                    if b[sub_u] in data.bayNum_of_Block[k]:
                        tmp.append(d)
                        break
                if d in tmp:
                    break
        D.append(tmp)
    # U_dk
    U = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
    U_bound = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
    b_U = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
    for d in range(len(data.destination)):
        for u in data.destination[d]:
            for sub_u in data.subGroup_dict[u]:
                for k in range(data.craneNum):
                    if b[sub_u] in data.bayNum_of_Block[k]:
                        U[k][d].append(sub_u)
                        b_U[k][d].append(b[sub_u])
                        break
    for k in range(data.craneNum):
        for d in D[k]:
            for u in U[k][d]:
                if (b[u] == data.bayNum_of_Block[k][0]) or (b[u] == data.bayNum_of_Block[k][-1]) or (
                        b[u] - 1 not in b_U[k][d]) or (b[u] + 1 not in b_U[k][d]):
                    U_bound[k][d].append(u)
    # L_kd
    L = {}
    for k in range(data.craneNum):
        for d in D[k]:
            tmp = [b[u] for u in U[k][d]]
            L[k, d] = max(tmp) - min(tmp)
    return b, D, U, U_bound, L


def U_dict(data, dominated_master_X):
    U = [[[] for _ in range(len(data.destination))] for _ in range(data.craneNum)]
    for key in dominated_master_X.keys():
        for d in range(len(data.destination)):
            for u in data.destination[d]:
                for k in range(data.craneNum):
                    if key[0] in data.subGroup_dict[u] and key[1] in data.bayNum_of_Block[k]:
                        U[k][d].append(key[0])
                        continue
                continue
            continue
        continue
    return U


def crane_symmetry(data, dominated_master_X, added_cuts, sub_obj_sum_k, allCrane):
    '''
    只针对各场桥工作区域的空贝位对称的情况
    :return:
    '''
    perm = itertools.permutations(list(range(data.craneNum)))
    longest_row = np.max([len(row) for row in data.bayNum_of_Block])
    for seq in perm:
        X = {}
        for key in dominated_master_X.keys():
            X[key[0], data.bayNum_of_Block[seq[key[1] // longest_row]][
                key[1] % longest_row]] = 1
            real_k = seq[key[1] // longest_row]
        if allCrane:
            added_cuts = added_cuts_permutation(X, added_cuts, sub_obj_sum_k, data.craneNum)
        else:
            added_cuts = added_cuts_permutation(X, added_cuts, sub_obj_sum_k, real_k)
    return added_cuts


def allocated_bay_set(data, master_result, crane):
    '''
    统计被分配到的贝位集合
    :return:
    '''
    dominated_master_X = {}
    max_j = 0
    j_set = []
    for key in master_result.keys():
        if (key[1] in data.bayNum_of_Block[crane]) and (master_result[key] == 1) and (key[0] < data.subGroupNum):
            max_j = max(max_j, key[1])
            dominated_master_X[key] = master_result[key]
            if key[1] not in j_set:
                j_set.append(key[1])
    j_set.sort()
    return j_set, max_j, dominated_master_X


def benders_MIP(data):
    LB = 0
    UB = float("inf")
    added_cuts = []
    master_results_history = []
    iteration_num = 0
    start_time = time.time()

    while 1:
        # step 2
        master_results = master_problem(data, added_cuts, LB)
        if not master_results:
            print("Master problem is infeasible!")
            return
        LB = max(LB, master_results[2])
        # print("master_x:", master_results[0])
        print("master_obj", master_results[2])
        print("LB:", LB)
        b, D, U, U_bound, L = sub_help(master_results[0], data)  # 主问题求解后、子问题求解前，计算一些辅助集合和参数
        iteration_num += 1

        # step 3
        sub_obj_sum = 0
        flag = 0
        sub_cranes_sum = []
        for w in range(data.scenarioNum):
            sub_status, sub_obj, sub_cranes = sub_problem(b, D, U, U_bound, L, data, w)
            if sub_status == 0:  # 子问题不可行
                # TODO: 可以根据不可行原因提出割
                added_cuts.append([master_results[0], 10000000000])
                break
            else:  # 子问题有最优解
                sub_obj_sum += sub_obj
                sub_cranes_sum.append(sub_cranes)
                if w == data.scenarioNum - 1:
                    flag = 1
        if flag == 1:  # 每个场景下的子问题都有解
            # =================== 单场桥割、同港互换割、平移对称割的交叉使用 ===================
            # 针对单场桥的部分分配
            for i in range(data.craneNum):
                sub_obj_sum_k = 0
                for w in range(data.scenarioNum):
                    sub_obj_sum_k += sub_cranes_sum[w][i]
                j_set, max_j, dominated_master_X = allocated_bay_set(data, master_results[0], i)
                # 对于各场桥区域贝位不对称的算例，仅在该场桥上添加同港互换割
                U = U_dict(data, dominated_master_X)
                added_cuts = inner_change(dominated_master_X, added_cuts, sub_obj_sum_k, U, i)
                # 对于各场桥区域贝位对称的算例，考虑消除对称性，并均使用同港互换割
                # TODO：与初始有效不等式中的消除对称性有部分重复作用，需更谨慎设计
                added_cuts = crane_symmetry(data, dominated_master_X, added_cuts, sub_obj_sum_k, False)

                # 平移对称割 + 同港互换割 交叉使用
                for g in range(len(j_set)):
                    same_j_set = j_set[:g]
                    change_j_set = j_set[g:]
                    for j in range(1, data.bayNum_of_Block[i][-1] - max_j):
                        new_dominated_master_X = {}
                        for key in dominated_master_X.keys():
                            if key[1] in same_j_set:
                                new_dominated_master_X[key[0], key[1]] = 1
                            if key[1] in change_j_set:
                                new_dominated_master_X[key[0], key[1] + j] = 1
                        U = U_dict(data, new_dominated_master_X)
                        added_cuts = inner_change(new_dominated_master_X, added_cuts, sub_obj_sum_k + j, U, i)
                        added_cuts = crane_symmetry(data, new_dominated_master_X, added_cuts, sub_obj_sum_k + j, False)

            # 针对完整解：平移对称割 + 同港互换割
            res = [[] for _ in range(data.craneNum)]
            for i in range(data.craneNum):
                j_set, max_j, dominated_master_X = allocated_bay_set(data, master_results[0], i)
                res[i].append(dominated_master_X)
                for g in range(len(j_set)):
                    same_j_set = j_set[:g]
                    change_j_set = j_set[g:]
                    for j in range(1, data.bayNum_of_Block[i][-1] - max_j):
                        new_dominated_master_X = {}
                        for key in dominated_master_X.keys():
                            if key[1] in same_j_set:
                                new_dominated_master_X[key[0], key[1]] = 1
                            if key[1] in change_j_set:
                                new_dominated_master_X[key[0], key[1] + j] = 1
                        res[i].append(new_dominated_master_X)
            # res维度 = 场桥个数 = 循环次数
            for i in res[0]:
                for j in res[1]:
                    result_dict = {}
                    result_dict.update(i)
                    result_dict.update(j)
                    U = U_dict(data, result_dict)
                    added_cuts = inner_change(result_dict, added_cuts, sub_obj_sum, U, data.craneNum)
                    added_cuts = crane_symmetry(data, result_dict, added_cuts, sub_obj_sum, True)
            UB = min(UB, sub_obj_sum + master_results[3])
            print("sub_status:", sub_status)
            print("sub_obj_sum:", sub_obj_sum)
            print("UB:", UB)

        if UB - LB < 0.00001:
            end_x = master_results[0]
            break
        else:
            print("=========================================我是分割线=========================================")
            print(UB)
            print(LB)
            print("iteration_num:", iteration_num)
            end_time = time.time()
            if end_time - start_time >= 7200:
                break

    print(end_x)
    print(UB)
    print("iteration_num:", iteration_num)
    return UB


class Data:
    # 算例基础数据
    containerNum_of_Group1 = [50] * 3  # 20TEU每个箱组的箱数
    containerNum_of_Group2 = [25] * 3  # 40TEU每个箱组的箱数
    craneNum = 2
    bayNum_of_Block = [list(range(6)), list(range(6, 12))]
    bayIndex = list(range(6)) + list(range(10000000, 10000006))
    gamma = [0.0002] * 9
    destination = [
        [0, 3],
        [1, 4],
        [2, 5]
    ]  # 每个卸货港的箱组

    s = 5
    t = 5
    w = 2
    t_move = 5.89 / 2.17  # 单位：秒
    t_relocation = 2.44 * s / 1.8 + 20  # 单位：秒
    t_pickup = 30  # 单位：秒
    bigM = 1000000
    delta_time = 1506
    gamma_dict = {}

    # 算例生成
    subGroupNum = None  # 子箱组数
    subGroup1 = []  # 20TEU子箱组
    subGroup2 = []  # 40TEU子箱组
    F = []
    subGroup_with_same_destination = []  # [[same destination subgroup indexes],[]]
    priority = None
    containerNum_of_subGroup = []
    coff_scenarios = None
    scenarioNum = None
    disMatrix = None
    bayNum = None
    subGroup_dict = {}  # group:

    # 生成每个目的港的子箱组
    def generate_subGroup_with_same_destination(self):
        for d in range(len(self.destination)):
            subGroup_of_d = []
            for group in self.destination[d]:
                subGroup_of_d.extend(self.subGroup_dict[group])
            self.subGroup_with_same_destination.append(subGroup_of_d)

    # 根据stacking结果选择gamma取值，目前直接在算例基础数据中直接给定
    def generate_gamma(self):
        for i in self.containerNum_of_subGroup:
            self.gamma.append(self.gamma_dict[i, ceil(i / self.t), self.t, self.w])

    # 输入箱组分子箱组
    def generate_subGroup(self):
        subGroupIndex = -1
        for j in range(len(self.containerNum_of_Group1)):
            group = self.containerNum_of_Group1[j]
            self.subGroup_dict[j] = []
            if group <= self.s * self.t:
                subGroupIndex += 1
                self.subGroup1.append(subGroupIndex)
                self.subGroup_dict[j].append(subGroupIndex)
                self.containerNum_of_subGroup.append(group)
            else:
                for i in range(ceil(group / (self.s * self.t))):
                    subGroupIndex += 1
                    self.subGroup_dict[j].append(subGroupIndex)
                    self.subGroup1.append(subGroupIndex)
                    if i == ceil(group / (self.s * self.t)) - 1:
                        self.containerNum_of_subGroup.append(
                            group % (self.s * self.t) if group % (self.s * self.t) != 0 else (self.s * self.t))
                    else:
                        self.containerNum_of_subGroup.append(self.s * self.t)
        for j in range(len(self.containerNum_of_Group2)):
            group = self.containerNum_of_Group2[j]
            j = j + len(self.containerNum_of_Group1)
            self.subGroup_dict[j] = []
            if group <= self.s * self.t:
                subGroupIndex += 1
                self.subGroup2.append(subGroupIndex)
                self.subGroup_dict[j].append(subGroupIndex)
                self.containerNum_of_subGroup.append(group)
            else:
                for i in range(ceil(group / (self.s * self.t))):
                    subGroupIndex += 1
                    self.subGroup_dict[j].append(subGroupIndex)
                    self.subGroup2.append(subGroupIndex)
                    if i == ceil(group / (self.s * self.t)) - 1:
                        self.containerNum_of_subGroup.append(
                            group % (self.s * self.t) if group % (self.s * self.t) != 0 else (self.s * self.t))
                    else:
                        self.containerNum_of_subGroup.append(self.s * self.t)
        self.subGroupNum = len(self.subGroup1) + len(self.subGroup2)

    # 输入贝位标号算距离矩阵
    def generate_disMatrix(self):
        self.bayNum = len(self.bayIndex)
        self.disMatrix = []
        for i in range(self.bayNum):
            distance = [0] * self.bayNum
            for j in range(self.bayNum):
                distance[j] = abs(self.bayIndex[i] - self.bayIndex[j])
            self.disMatrix.append(distance)

    # 生成各子箱组提取顺序
    def generate_priority(self):
        if len(self.destination) <= 3:
            perm = itertools.permutations(list(range(len(self.destination))))  # res = []
            self.priority = []
            for s in perm:
                p_matrix = [[-1] * len(self.destination) for _ in range(len(self.destination))]
                for j in range(len(s) - 1, -1, -1):
                    if j < len(s) - 1:
                        p_matrix[s[j]] = p_matrix[s[j + 1]].copy()
                        p_matrix[s[j]][s[j + 1]] = 1
                    p_matrix[s[j]][s[j]] = 0
                self.priority.append(p_matrix)
            print(self.priority)
            self.coff_scenarios = [1 / len(self.priority)] * len(self.priority)
            self.scenarioNum = len(self.priority)
        else:
            # 随机抽三个为一组，进行划分
            perm = self.extract_groups(list(range(len(self.destination))))
            self.priority = []
            for s in perm:
                p_matrix = [[-1] * len(self.destination) for _ in range(len(self.destination))]
                for j in range(len(s) - 1, -1, -1):
                    if j < len(s) - 1:
                        p_matrix[s[j]] = p_matrix[s[j + 1]].copy()
                        p_matrix[s[j]][s[j + 1]] = 1
                    p_matrix[s[j]][s[j]] = 0
                self.priority.append(p_matrix)
            print(self.priority)
            self.coff_scenarios = [1 / len(self.priority)] * len(self.priority)
            self.scenarioNum = len(self.priority)

    def extract_groups(self, destination):
        result = []
        for combo in combinations(destination, 3):
            remaining = [num for num in destination if num not in combo]
            if len(remaining) <= 3:
                result.append(list(combo) + remaining)
                continue
            for l in self.extract_groups(remaining):
                result.append(list(combo) + l)
        return result

    def generate_F(self):
        for i in range(len(self.bayIndex)):
            if self.bayIndex[i] + 1 in self.bayIndex:
                self.F.append(i)


if __name__ == "__main__":
    data = Data()
    data.generate_subGroup()  # 20TEU和40TEU分开存放
    data.generate_disMatrix()  # 计算距离
    data.generate_F()  # 标记连续贝位
    data.generate_priority()
    data.generate_subGroup_with_same_destination()
    start_time = time.time()

    benders_MIP(data)
    end_time = time.time()
    run_time = end_time - start_time
    print("运行时间：", run_time, "秒")
