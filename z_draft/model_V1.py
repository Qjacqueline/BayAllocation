# -*- coding: utf-8 -*-
# @Time    : 2022/6/11 6:35 PM
# @Author  : JacQ
# @File    : model_V1.py

from gurobipy import *

from a_data_process.read_data import Data


def model(data: Data):
    # ============== 构造模型 ================
    big_M = 10000
    model_V1 = Model("bay_allocation")

    # ============== 定义变量 ================
    # X_ij
    X = [[[] for _ in range(data.U_num)] for _ in range(data.J_num)]
    for i in range(data.J_num):
        for j in range(data.U_num):
            name = 'X_' + str(i) + "_" + str(j)
            X[i][j] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name)

    # Y_jn
    Y = [[[] for _ in range(data.N_num)] for _ in range(data.J_num)]
    for j in range(data.U_num):
        for n in range(data.N_num):
            name = 'Y_' + str(j) + "_" + str(n)
            Y[j][n] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name)

    # Z_jj'n
    Z = [[[[] for _ in range(data.N_num)] for _ in range(data.U_num)] for _ in range(data.U_num)]
    for j in range(data.U_num):
        for jj in range(data.U_num):
            for n in range(data.N_num):
                name = 'Z_' + str(j) + "_" + str(jj) + "_" + str(n)
                Z[j][jj][n] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name)
    # FTjn
    FT = [[[] for _ in range(data.N_num)] for _ in range(data.U_num)]
    for j in range(data.U_num):
        for n in range(data.N_num):
            name = 'FT_' + str(j) + "_" + str(n)
            FT[j][n] = model_V1.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name=name)

    # ============== 构造目标 ================
    q = model_V1.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='q')  # 线性化模型变量
    # model.setObjective(q, GRB.MINIMIZE)
    model_V1.addConstrs((q >= FT[j][n] for j in range(data.U_num) for n in range(data.N_num)), "obj1")

    obj2 = LinExpr(0)
    M_1 = [[[] for _ in range(data.U_num)] for _ in range(data.U_num)]
    for j in range(data.U_num):
        for jj in range(data.U_num):
            name = 'M1_' + str(j) + "_" + str(jj)
            M_1[j][jj] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name)

    for i in range(data.J_num):
        for j in range(data.U_num):
            for jj in range(data.U_num):
                model_V1.addConstr(M_1[j][jj] <= X[i][j], 'M1_1')
                model_V1.addConstr(M_1[j][jj] <= X[i][jj], 'M1_2')
                model_V1.addConstr(M_1[j][jj] >= X[i][j] + X[i][jj] - 1, 'M1_3')
                obj2.addTerms(data.gamma_uu[j][jj], M_1[j][jj])

    # model.setObjective(obj2, GRB.MINIMIZE)
    model_V1.setObjective(q + obj2, GRB.MINIMIZE)

    # ============== 构造约束 ================
    # 对于一个贝
    # con1: 一个贝上放置的箱组一般不超过2个
    model_V1.addConstrs((sum(X[i][j] for j in range(data.U_num)) <= 2 for i in range(data.J_num)), "con1")
    # con2: 放置集装箱数不超过贝位的最大容量
    model_V1.addConstrs(
        (sum(data.U_L_num_set[j] * X[i][j] for j in range(data.U_num)) <= data.S_num * data.T_num for i in range(data.J_num)),
        "con2")

    # 对于一个子箱组
    # con3: 每个子箱组必须放到一个贝上
    model_V1.addConstrs((sum(X[i][j] for i in range(data.J_num)) == 1 for j in range(data.U_num)), "con3")
    # con4: 每个子箱组必须由一个岸桥服务
    model_V1.addConstrs((sum(Y[j][n] for n in range(data.N_num)) == 1 for j in range(data.U_num)), "con4")

    # 对于40ft的子箱组 TODO

    # 对于岸桥n
    # con5: 如果j,j'是相邻任务，则完成j后马上完成j'
    unit_mt = 5
    M_2 = [[[[[] for _ in range(data.U_num)] for _ in range(data.U_num)] for _ in range(data.J_num)] for _ in range(data.J_num)]
    for i in range(data.J_num):
        for ii in range(data.J_num):
            for j in range(data.U_num):
                for jj in range(data.U_num):
                    name = 'M2_' + str(i) + "_" + str(ii) + str(j) + "_" + str(jj)
                    M_2[i][ii][j][jj] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name)

    # mt_jj'
    mt = [[[] for _ in range(data.U_num)] for _ in range(data.U_num)]
    for j in range(data.U_num):
        for jj in range(data.U_num):
            name = 'mt_' + str(j) + "_" + str(jj)
            mt[j][jj] = model_V1.addVar(lb=0, ub=big_M, vtype=GRB.CONTINUOUS, name=name)
            expr = LinExpr(0)
            for i in range(data.J_num):
                for ii in range(data.J_num):
                    model_V1.addConstr(M_2[i][ii][j][jj] <= X[i][j], 'M2_1')
                    model_V1.addConstr(M_2[i][ii][j][jj] <= X[ii][jj], 'M2_2')
                    model_V1.addConstr(M_2[i][ii][j][jj] >= X[i][j] + X[ii][jj] - 1, 'M2_3')
                    expr.addTerms(5 * data.d_jj[i][ii], M_2[i][ii][j][jj])
            model_V1.addConstr(mt[j][jj] == expr, 'define mt_jj')
            expr.clear()

    for j in range(data.U_num):
        for jj in range(data.U_num):
            for n in range(data.N_num):
                model_V1.addConstr(Z[j][jj][n] * (FT[j][n] + mt[j][jj] + data.t_u[jj] - FT[jj][n]) == 0, 'con5')

    # con6: 提取需满足优先级
    M_3 = [[[] for _ in range(data.U_num)] for _ in range(data.U_num)]
    M_4 = [[[] for _ in range(data.U_num)] for _ in range(data.U_num)]

    for j in range(data.U_num):
        for jj in range(data.U_num):
            name1 = 'M3_' + str(j) + "_" + str(jj)
            M_3[j][jj] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name1)
            name2 = 'M4_' + str(j) + "_" + str(jj)
            M_4[j][jj] = model_V1.addVar(0, 1, vtype=GRB.BINARY, name=name2)
            model_V1.addConstr(data.p_u[jj] - data.p_u[j] <= big_M * (1 - M_3[j][jj]), 'con6_1')  # FIXME
            model_V1.addConstrs((FT[jj][n] - FT[j][n] <= big_M * (1 - M_4[j][jj]) for n in range(data.N_num)), 'con6_2')
            model_V1.addConstr(M_4[j][jj] >= M_3[j][jj], 'con6_3')

    # con7: 保证安全距离 如果两贝之间距离小于安全距离，则操作时间不能重合
    # model_V1.addConstrs((X[i][j] * X[ii][jj] * (data.F - data.d_ii[i][ii])
    #                      <= big_M * (FT[j][n] - data.t_j[j] - FT[jj][nn])
    #                      * (FT[jj][nn] - data.t_j[jj] - FT[j][n])
    #                      for i in range(data.I) for ii in range(data.I)
    #                      for j in range(data.J) for jj in range(data.J)
    #                      for n in range(data.N) for nn in range(data.N)), 'con7')

    # 变量间关系
    # con8: Xij和Yjn
    for j in range(len(data.J)):
        ac_index = data.J[j]
    cor_k = math.floor(ac_index / 120)  # 对应所在的箱区
    model_V1.addConstrs(((X[i][j] - 1) * (Y[j][cor_k] + Y[j][cor_k + 1] - 1) == 0
                         for i in range(data.J_num) for j in range(data.U_num)), 'con8')
    # con9: Yjn和Zjj'n
    model_V1.addConstrs((Z[j][jj][n] <= big_M * Y[j][n]
                         for j in range(data.U_num) for jj in range(data.U_num)
                         for n in range(data.N_num)), 'con9')
    model_V1.addConstrs((Z[j][jj][n] <= big_M * Y[jj][n]
                         for j in range(data.U_num) for jj in range(data.U_num)
                         for n in range(data.N_num)), 'con10')

    model_V1.write('a.lp')
    model_V1.Params.timelimit = 3600
    model_V1.optimize()
    if model_V1.status == GRB.Status.INFEASIBLE:
        print('Optimization was stopped with status %d' % model_V1.status)
        # do IIS, find infeasible constraints
        model_V1.computeIIS()
        for c in model_V1.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
