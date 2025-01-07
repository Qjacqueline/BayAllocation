# -*- coding: utf-8 -*-
# @Time    : 2024/12/20 8:41 AM
# @Author  : JacQ
# @File    : OP_hepu.py
import itertools
import time
from math import ceil
from random import random

import gurobipy
import numpy as np
from gurobipy import GRB, LinExpr, quicksum


def permute(nums):
    res = []
    def backtrack(nums, tmp):
        if not nums:
            res.append(tmp)
            return
        for i in range(len(nums)):
            backtrack(nums[:i] + nums[i + 1:], tmp + [nums[i]])

    backtrack(nums, [])
    return res


class Data:
    # 算例基础数据
    containerNum_of_Group1 = [50] * 3  # 20TEU每个箱组的箱数
    containerNum_of_Group2 = [25] * 3  # 40TEU每个箱组的箱数
    craneNum = 2
    bayNum_of_Block = [list(range(6)), list(range(6, 12))]
    bayIndex = list(range(6)) + list(range(1000000, 1000006))
    gamma = [0.0002] * 18
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
    delta_time = 1510
    gamma_dict = {}

    # 算例生成
    subGroupNum = None  # 子箱组数
    subGroup1 = []  # 20TEU子箱组
    subGroup2 = []  # 40TEU子箱组
    F = []  # 后一连续贝位可放置的贝位下标集合

    priority = None
    containerNum_of_subGroup = []
    coff_scenarios = None
    scenarioNum = None
    disMatrix = None
    bayNum = None
    subGroup_dict = {}
    subGroup_with_same_destination = []

    def generate_subGroup_with_same_destination(self):
        for d in range(len(self.destination)):
            subGroup_of_d = []
            for group in self.destination[d]:
                subGroup_of_d.extend(self.subGroup_dict[group])
            self.subGroup_with_same_destination.append(subGroup_of_d)

    # 根据stacking结果选择gamma取值，目前直接在算例基础数据中直接给定
    def generate_gamma(self):
        for i in self.containerNum_of_subGroup:
            self.gamma.append(self.gamma_dict[i, ceil(i/self.t), self.t, self.w])

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
                for i in range(ceil(group/(self.s * self.t))):
                    subGroupIndex += 1
                    self.subGroup_dict[j].append(subGroupIndex)
                    self.subGroup1.append(subGroupIndex)
                    if i == ceil(group/(self.s * self.t)) - 1:
                        self.containerNum_of_subGroup.append(group % (self.s * self.t) if group % (self.s * self.t) != 0 else (self.s * self.t))
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
                for i in range(ceil(group/(self.s * self.t))):
                    subGroupIndex += 1
                    self.subGroup_dict[j].append(subGroupIndex)
                    self.subGroup2.append(subGroupIndex)
                    if i == ceil(group/(self.s * self.t)) - 1:
                        self.containerNum_of_subGroup.append(group % (self.s * self.t) if group % (self.s * self.t) != 0 else (self.s * self.t))
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
        for combo in itertools.combinations(destination, 3):
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


if __name__ == '__main__':
    data = Data()
    data.generate_subGroup()
    data.generate_priority()
    data.generate_disMatrix()
    data.generate_F()
    data.generate_subGroup_with_same_destination()
    start_time = time.time()
    # ============= 创建模型 & 参数定义 =============
    model = gurobipy.Model('贝位分配模型')  # 创建一个模型

    # ============= 决策变量定义 =============
    # 定义字典用来存放决策变量
    X = {}
    Y = {}
    C = {}
    C_max = {}
    Z = {}
    A = {}
    B = {}

    # A, B：每个场景下、每个卸货港子箱组的开始提取时间和结束提取时间
    for w in range(data.scenarioNum):
        for d in range(len(data.destination)):
            A[w, d] = model.addVar(0, float("inf"), vtype=GRB.CONTINUOUS, name='A_' + str(w) + '_' + str(d))
            B[w, d] = model.addVar(0, float("inf"), vtype=GRB.CONTINUOUS, name='B_' + str(w) + '_' + str(d))
    # X
    for u in range(data.subGroupNum):
        for j in range(data.bayNum):
            name = 'X_' + str(u) + '_' + str(j)
            X[u, j] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, u放在j上则为1

    for k in range(data.craneNum):
        j = data.bayNum_of_Block[k][0]
        X[data.subGroupNum, j] = model.addVar(1, 1, vtype=GRB.BINARY, name='X_' + str(data.subGroupNum) + '_' + str(j))  # 0-1变量, u放在j上则为1
    # S
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            for u1 in range(data.subGroupNum):
                for u2 in range(data.subGroupNum):
                    if u1 != u2:
                        name = 'Y_' + str(w) + '_' + str(k) + '_' + str(u1) + '_' + str(u2)
                        Y[w, k, u1, u2] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, 场景w下，场桥k上u1和u2连续提取则为1
    # 初始节点
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            Y[w, k, data.subGroupNum, data.subGroupNum + 1] = model.addVar(0, 1, vtype=GRB.BINARY, name='S_' + str(w) + '_' + str(k) + '_' + str(data.subGroupNum) + '_' + str(data.subGroupNum + 1))
            for u in range(data.subGroupNum):
                name = 'Y_' + str(w) + '_' + str(k) + '_' + str(data.subGroupNum) + '_' + str(u)
                Y[w, k, data.subGroupNum, u] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, 场景w下，场桥k上u1和u2连续提取则为1
    # 末尾节点
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            Y[w, k, data.subGroupNum + 1, data.subGroupNum] = model.addVar(0, 1, vtype=GRB.BINARY, name='S_' + str(w) + '_' + str(k) + '_' + str(data.subGroupNum + 1) + '_' + str(data.subGroupNum))
            for u in range(data.subGroupNum):
                name = 'Y_' + str(w) + '_' + str(k) + '_' + str(u) + '_' + str(data.subGroupNum+1)
                Y[w, k, u, data.subGroupNum + 1] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, 场景w下，场桥k上u1和u2连续提取则为1
    # C
    for w in range(data.scenarioNum):
        C[w, data.subGroupNum] = model.addVar(0, 0, vtype=GRB.CONTINUOUS, name='C_' + str(w) + '_' + str(data.subGroupNum))  # 连续变量, 场景w下，u的提取结束时间
        for u in range(data.subGroupNum):
            name = 'C_' + str(w) + '_' + str(u)
            C[w, u] = model.addVar(0, 10000000, vtype=GRB.CONTINUOUS, name=name)  # 连续变量, 场景w下，u的提取结束时间
    for w in range(data.scenarioNum):
        name = 'C_max' + '_' + str(w)
        C_max[w] = model.addVar(0, 1000000000, vtype=GRB.CONTINUOUS, name=name)  # 连续变量, 场景w下的makespan

    # Z
    for j1 in range(data.bayNum):
        for j2 in range(data.bayNum):
            for u1 in range(data.subGroupNum):
                for u2 in range(data.subGroupNum):
                    if u1 < u2:
                        name = 'Z_' + str(u1) + '_' + str(u2) + '_' + str(j1) + '_' + str(j2)
                        Z[u1, u2, j1, j2] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, u1和u2同时放在j上则为1
    for k in range(data.craneNum):
        for j in data.bayNum_of_Block[k]:
            for u in range(data.subGroupNum):
                name = 'Z_' + str(data.subGroupNum) + '_' + str(u) + '_' + str(data.bayNum_of_Block[k][0]) + '_' + str(j)
                Z[data.subGroupNum, u, data.bayNum_of_Block[k][0], j] = model.addVar(0, 1, vtype=GRB.BINARY, name=name)  # 0-1变量, u1和u2同时放在j上则为1
    # ============= 目标函数 =============
    obj = LinExpr(0)
    for i in range(data.scenarioNum):
        obj.addTerms(data.coff_scenarios[i], C_max[i])
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)

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
            lhs.addTerms(ceil(data.containerNum_of_subGroup[u]/data.t), X[u, j])
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
    # makespan >= 每个u的提取结束时间
    for scenario in range(data.scenarioNum):
        for u in range(data.subGroupNum):
            model.addConstr(C[scenario, u] <= C_max[scenario])
    # S == 1, 则u1和u2必在k上
    for scenario in range(data.scenarioNum):
        for k in range(data.craneNum):
            for u1 in range(data.subGroupNum):
                for u2 in range(data.subGroupNum):
                    if u1 != u2:
                        lhs = LinExpr(0)
                        for j in data.bayNum_of_Block[k]:
                            lhs.addTerms(1, X[u1, j])
                            lhs.addTerms(1, X[u2, j])
                        # lhs.addTerms(data.bigM)
                        lhs.addTerms(- data.bigM, Y[scenario, k, u1, u2])
                        model.addConstr(lhs + data.bigM >= 2)
    # 前任务只有一个连续
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            for u2 in range(data.subGroupNum):
                rhs = LinExpr(0)
                for j in data.bayNum_of_Block[k]:
                    rhs.addTerms(1, X[u2, j])
                lhs = LinExpr(0)
                for u1 in range(data.subGroupNum+1):
                    if u1 != u2:
                        lhs.addTerms(1, Y[w, k, u1, u2])
                model.addConstr(lhs == rhs)
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            lhs = LinExpr(0)
            for u1 in range(data.subGroupNum + 1):
                lhs.addTerms(1, Y[w, k, u1, data.subGroupNum + 1])
            model.addConstr(lhs == 1)
    # 后任务只有一个
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            for u1 in range(data.subGroupNum):
                rhs = LinExpr(0)
                for j in data.bayNum_of_Block[k]:
                    rhs.addTerms(1, X[u1, j])
                lhs = LinExpr(0)
                for u2 in range(data.subGroupNum):
                    if u1 != u2:
                        lhs.addTerms(1, Y[w, k, u1, u2])
                lhs.addTerms(1, Y[w, k, u1, data.subGroupNum + 1])
                model.addConstr(lhs == rhs)
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            lhs = LinExpr(0)
            for u2 in range(data.subGroupNum):
                lhs.addTerms(1, Y[w, k, data.subGroupNum, u2])
            lhs.addTerms(1, Y[w, k, data.subGroupNum, data.subGroupNum + 1])
            model.addConstr(lhs == 1)
    # 前后集装箱提取时间关系
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            for u1 in range(data.subGroupNum):
                for u2 in range(data.subGroupNum):
                    if u1 != u2:
                        rhs = LinExpr(0)
                        for j1 in data.bayNum_of_Block[k]:
                            for j2 in data.bayNum_of_Block[k]:
                                rhs.addTerms(data.t_move * data.disMatrix[j1][j2], Z[min(u1, u2), max(u1, u2), j1, j2])
                        if u2 in data.subGroup1:
                            model.addConstr(C[w, u2] - C[w, u1] + data.bigM * (1 - Y[w, k, u1, u2]) >= rhs + data.t_relocation * data.gamma[u2] + data.containerNum_of_subGroup[u2] * data.t_pickup)
                        if u2 in data.subGroup2:
                            model.addConstr(
                                C[w, u2] - C[w, u1] + data.bigM * (1 - Y[w, k, u1, u2]) >= rhs + 2 * data.t_relocation *
                                data.gamma[u2] + 2 * data.containerNum_of_subGroup[u2] * data.t_pickup)
    for w in range(data.scenarioNum):
        for k in range(data.craneNum):
            u1 = data.subGroupNum
            for u2 in range(data.subGroupNum):
                if u1 != u2:
                    rhs = LinExpr(0)
                    j1 = data.bayNum_of_Block[k][0]
                    for j2 in data.bayNum_of_Block[k]:
                        rhs.addTerms(data.t_move * data.disMatrix[j1][j2], X[u2, j2])
                    if u2 in data.subGroup1:
                        model.addConstr(C[w, u2] - C[w, u1] + data.bigM * (1 - Y[w, k, u1, u2]) >= rhs + data.t_relocation * data.gamma[u2] + data.containerNum_of_subGroup[u2] * data.t_pickup)
                    if u2 in data.subGroup2:
                        model.addConstr(
                            C[w, u2] - C[w, u1] + data.bigM * (1 - Y[w, k, u1, u2]) >= rhs + 2 * data.t_relocation *
                            data.gamma[u2] + 2 * data.containerNum_of_subGroup[u2] * data.t_pickup)
    # 提取优先级与完成提取的时间关系
    sub_u = [[] for _ in range(len(data.destination))]
    for d in range(len(data.destination)):
        for u in data.destination[d]:
            sub_u[d].extend(data.subGroup_dict[u])
    # A、B相关约束
    for scenario in range(data.scenarioNum):
        for d in range(len(data.destination)):
            for u in sub_u[d]:
                if u in data.subGroup1:
                    model.addConstr(A[scenario, d] <= C[scenario, u] - data.t_pickup * data.containerNum_of_subGroup[u] - data.t_relocation * data.gamma[u])
                if u in data.subGroup2:
                    model.addConstr(A[scenario, d] <= C[scenario, u] - 2 * data.t_pickup * data.containerNum_of_subGroup[u] - 2 * data.t_relocation * data.gamma[u])
                model.addConstr(B[scenario, d] >= C[scenario, u])
    for scenario in range(data.scenarioNum):
        for d in range(len(data.destination)):
                model.addConstr(B[scenario, d] - A[scenario, d] <= data.delta_time)
    for scenario in range(data.scenarioNum):
        for d1 in range(len(data.destination)):
            for d2 in range(len(data.destination)):
                if data.priority[scenario][d1][d2] == 1:
                    for u1 in sub_u[d1]:
                        for u2 in sub_u[d2]:
                            model.addConstr(C[scenario, u2] >= C[scenario, u1])

    # 40TEU所在贝位的后一贝位不放子箱组
    for j in data.F:
        for u in data.subGroup2:
            lhs = LinExpr(0)
            lhs.addTerms(data.bigM, X[u, j])
            for u_temp in range(data.subGroupNum):
                lhs.addTerms(1, X[u_temp, j+1])
            model.addConstr(lhs - data.bigM <= 0)
    # Z的约束
    for j1 in range(data.bayNum):
        for j2 in range(data.bayNum):
            for u1 in range(data.subGroupNum):
                for u2 in range(data.subGroupNum):
                    if u1 < u2:
                        model.addConstr(2 * Z[u1, u2, j1, j2] <= X[u1, j1] + X[u2, j2])
                        model.addConstr(Z[u1, u2, j1, j2] >= X[u1, j1] + X[u2, j2] - 1)

    # ============= 模型求解 & 结果打印 =============
    model.optimize()
    print("\n\n-----optimal value-----")
    print(model.ObjVal)

    end_time = time.time()
    run_time = end_time - start_time
    print("运行时间：", run_time, "秒")

    for key in X.keys():
        print(X[key].VarName + ' = ', X[key].x)
    for key in Y.keys():
        print(Y[key].VarName + ' = ', Y[key].x)
    for key in C.keys():
        print(C[key].VarName + ' = ', C[key].x)