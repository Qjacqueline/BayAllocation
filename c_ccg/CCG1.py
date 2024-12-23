# -*- coding: utf-8 -*-
# @Time    : 2024/12/18 3:47 PM
# @Author  : JacQ
# @File    : CCG1.py
import time
from gurobipy import *
from a_data_process.read_data import read_data
import a_data_process.config as cf


def CCG(data):
    LB = 0
    UB = float("inf")
    iteration_num = 0
    start_time = time.time()
    added_cuts = []

    while 1:
        print("=========================================iteration\t" +
              str(iteration_num) + "=========================================")

        # 求解主问题
        master_results = master_problem(data, added_cuts)
        if not master_results:
            print("Master problem is infeasible!")
            return
        LB = max(LB, master_results[1])  # 当前迭代主问题目标值
        print("Master obj:\t", master_results[1])

        # 预处理子问题
        N, pos, pt = sub_problem_help(data, master_results[0])

        # 求解子问题
        sub_results = sub_problem(N, pos, pt)
        UB = min(UB, sub_results[0])

        # 添加cuts
        added_cuts = generate_cuts(master_results, sub_results, added_cuts)

        # 判断迭代终止条件
        if UB - LB < 0.00001:
            end_x = master_results[0]
            break
        else:
            end_time = time.time()
            print("UB:\t", UB)
            print("LB:\t", LB)
            print("Time:\t", end_time - start_time)
            if end_time - start_time >= 7200:
                break
        iteration_num += 1

    print("=========================================Final results=========================================")
    print("Final solution:\t", end_x)
    print("Final UB:\t", UB)
    print("Final LB:\t", LB)
    return UB


def master_problem(data, added_cuts):
    """
        :param data: 数据
        :param added_cuts: 所有cut的集合
        :return: 不可行，返回 FALSE；
                 可行，返回 master_X, theta.x
                   分别表示 解、theta值
    """
    # ============== 构造模型 ================
    big_M = 10000
    model = Model("Master problem")

    # ============== 定义变量 ================
    # X_uj: j从0开始 uN+k是虚拟job
    X = [[[] for _ in range(60 * data.K_num)] for _ in range(data.U_num + 2)]
    for j in range(0, 60 * data.K_num):
        for u in range(data.U_num + 2):
            if j + 1 not in data.J:
                X[u][j] = 0
                continue
            X[u][j] = model.addVar(0, 1, vtype=GRB.BINARY, name='X_' + str(u) + "_" + str(j))

    # ================ 约束 ==================
    # 对于一个子箱组
    # con1: Each sub-container group must be placed on one bay
    model.addConstrs((quicksum(X[u][j - 1] for j in data.J) == 1 for u in data.U_L), "2b")
    model.addConstrs((quicksum(X[u][j - 1] for j in data.I) == 1 for u in data.U_F), "2c")
    # con2: Initial position restrictions
    model.addConstrs((X[data.U_num][j - 1] == 1 for j in data.J_K_first), "2d")
    # con3:对于40ft的子箱组占了前一个后一个位置就不能被其他使用
    model.addConstrs((X[u][j - 1] <= big_M * (1 - X[uu][j + 1]) for u in data.U_F for uu in data.U
                      for j in data.I), "2e")
    # 对于一个贝
    # con4: 一个贝上放置的箱组一般不超过2个
    model.addConstrs((quicksum(X[u][j - 1] for u in data.U) <= 2 for j in data.J), "2f")
    # con5: 20和40的不能放在一个贝
    model.addConstrs((X[u][j - 1] + X[uu][j - 1] <= 1 for j in data.J for u in data.U_L for uu in data.U_F), "2g")
    # con6: 放置集装箱数不超过贝位的最大容量
    model.addConstrs((quicksum(data.U_num_set[u] * X[u][j - 1] for u in data.U) <= data.S_num * data.T_num
                      for j in data.J), "2h")

    # ============== 构造目标 ================
    theta = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='theta')  # 线性化模型变量
    obj = LinExpr(0)
    obj.addTerms(1, theta)
    model.setObjective(obj, gurobipy.GRB.MINIMIZE)

    # ============== 添加global cuts ================
    # extreme bay indices
    A = [[model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='A' + str(k) + "_" + str(g))
          for g in range(data.G_num)] for k in range(data.K_num)]
    B = [[model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='B' + str(k) + "_" + str(g))
          for g in range(data.G_num)] for k in range(data.K_num)]
    AA = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='AA' + str(k)) for k in range(data.K_num)]
    BB = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='BB' + str(k)) for k in range(data.K_num)]

    model.addConstrs((AA[k] <= j * X[u][j - 1] for k in range(data.K_num) for j in data.J_K[k] for u in data.U), "AA")
    model.addConstrs((BB[k] >= j * X[u][j - 1] for k in range(data.K_num) for j in data.J_K[k] for u in data.U), "BB")
    model.addConstrs((theta + cf.unit_move_time * AA[k] >= cf.unit_move_time * BB[k]
                      + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1] for u in data.U for j in data.J_K[k])
                      for k in range(data.K_num)), "global lb")

    model.addConstrs((A[k][g] <= j * X[u][j - 1] for k in range(data.K_num) for j in data.J_K[k]
                      for g in range(data.G_num) for u in data.U_g[g]), "A")
    model.addConstrs((B[k][g] >= j * X[u][j - 1] for k in range(data.K_num) for j in data.J_K[k]
                      for g in range(data.G_num) for u in data.U_g[g]), "B")
    model.addConstrs((theta + sum(cf.unit_move_time * A[k][g] for g in range(data.G_num)) >=
                      sum(cf.unit_move_time * B[k][g] for g in range(data.G_num))
                      + sum(cf.unit_process_time * data.U_num_set[u] * X[u][j - 1] for u in data.U for j in data.J_K[k])
                      for k in range(data.K_num)), "global lb")

    # ============== 添加iter cuts ================
    for cut in added_cuts:
        if cut[0] == 'Feasible':
            v, XX = cut[1], cut[2]
            rhs = LinExpr(0)
            for u, j in XX:
                rhs.addTerms(1, X[u][j])
            model.addConstr(theta >= big_M * (rhs - len(XX) + 1), "fea cut")
        elif cut[0] == 'Optimal' or cut[0] == 'Translation':
            v, XX = cut[1], cut[2]
            rhs = LinExpr(0)
            for u, j in XX:
                rhs.addTerms(1, X[u][j])
            model.addConstr(theta >= v * (rhs - len(XX) + 1), "opt cut")

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
        master_X = {}
        for u in range(data.U_num):
            for j in range(60 * data.K_num):
                if j + 1 not in data.J:
                    continue
                if abs(X[u][j].X) > 0.00001:
                    master_X[u] = j
        return master_X, theta.X


def sub_problem_help(data, master_X):
    """
        :param data: 数据
        :param master_X: 主问题获得的解
        :return: N 几个箱组
                pos 对应A，B子箱组位置
    """
    G_u_pos = [[] for _ in range(data.G_num)]  # 每个箱组子箱组位置
    for u in range(data.U_num):
        j = master_X[u]
        G_u_pos[data.U_g_set[u]].append(j)
    pos = [[min(G_u_pos[g]) * cf.unit_move_time, max(G_u_pos[g]) * cf.unit_move_time]
           for g in range(data.D_num * 2)]  # 每个箱组AB子箱组位置
    pt = [g_num * cf.unit_process_time for g_num in data.G_num_set]
    return data.D_num * 2, pos, pt


def sub_problem(N, pos, pt):
    """
        :param N: 几个箱组
        :param pos: 对应A，B子箱组位置
        :param pt: 每个箱组所需加工时间
        :return: (1/0表示是否可行, 目标值)
             不可行 (0，float("inf"))
             可行 (1，目标值)
    """

    def f_fuc(i, k, ii, jj):
        v = max(min(dp[i][k][j][h] + abs(pos[ii - 1][1 - jj] - pos[i - 1][j])
                    for j in [0, 1]) for h in [0, 1]) + abs(pos[ii - 1][0] - pos[ii - 1][1]) + pt[ii - 1]
        return v

    def g_fuc(i, k, ii, jj):
        v = min(dp[i][k][j][0] + abs(pos[ii - 1][1 - jj] - pos[i - 1][j])
                for j in [0, 1]) + abs(pos[ii - 1][0] - pos[ii - 1][1]) + pt[ii - 1]
        return v

    # Initialize DP table and trace table
    dp = [[[[0, 0] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 3:i,k,j | 4 i,k,j,h
    path = [[[[] for _ in range(2)] for _ in range(3)] for _ in range(N + 1)]  # 跟踪路径

    # Init state
    ## k=1
    dp[1][1][0][0] = abs(pos[0][1]) + abs(pos[0][1] - pos[0][0]) + pt[0]  # B_1 -> A_1
    dp[1][1][1][0] = abs(pos[0][0]) + abs(pos[0][0] - pos[0][1]) + pt[1]  # A_1 -> B_1

    if N > 1:
        dp[2][0][0][0] = abs(pos[1][1]) + abs(pos[1][0] - pos[1][1]) + pt[1]  # A_2 -> B_2
        dp[2][0][1][0] = abs(pos[1][0]) + abs(pos[1][1] - pos[1][0]) + pt[1]  # B_2 -> A_2

        ## k=2
        dp[1][2][0][0] = g_fuc(i=2, k=0, ii=1, jj=0)
        dp[1][2][1][0] = g_fuc(i=2, k=0, ii=1, jj=1)
        dp[2][1][0][0] = g_fuc(i=1, k=1, ii=2, jj=0)
        dp[2][1][1][0] = g_fuc(i=1, k=1, ii=2, jj=1)

    if N > 2:
        dp[3][0][0][0] = g_fuc(i=1, k=1, ii=3, jj=0)
        dp[3][0][1][0] = g_fuc(i=1, k=1, ii=3, jj=1)
        ## k=3
        dp[2][2][0][0] = g_fuc(i=3, k=0, ii=2, jj=0)
        dp[2][2][1][0] = g_fuc(i=3, k=0, ii=2, jj=1)
        dp[3][1][0][0] = f_fuc(i=1, k=2, ii=3, jj=0)
        dp[3][1][1][0] = f_fuc(i=1, k=2, ii=3, jj=1)

    if N > 3:
        dp[4][0][0] = [g_fuc(i=1, k=2, ii=4, jj=0), g_fuc(i=2, k=1, ii=4, jj=0)]
        dp[4][0][1] = [g_fuc(i=1, k=2, ii=4, jj=1), g_fuc(i=2, k=1, ii=4, jj=1)]

    # DP transitions
    for i in range(4, N + 1):

        dp[i - 1][2][0][0] = f_fuc(i=i, k=0, ii=i - 1, jj=0)
        dp[i - 1][2][1][0] = f_fuc(i=i, k=0, ii=i - 1, jj=1)

        dp[i][1][0] = [g_fuc(i=i - 2, k=2, ii=i, jj=0), f_fuc(i=i - 1, k=1, ii=i, jj=0)]
        dp[i][1][1] = [g_fuc(i=i - 2, k=2, ii=i, jj=1), f_fuc(i=i - 1, k=1, ii=i, jj=1)]

        if i < N:
            dp[i + 1][0][0] = [f_fuc(i=i - 2, k=2, ii=i + 1, jj=0), f_fuc(i=i - 1, k=1, ii=i + 1, jj=0)]
            dp[i + 1][0][1] = [f_fuc(i=i - 2, k=2, ii=i + 1, jj=1), f_fuc(i=i - 1, k=1, ii=i + 1, jj=1)]

    max_value = max(max(min(dp[N][1][0][0], dp[N][1][1][0]), min(dp[N][1][0][1], dp[N][1][1][1])),
                    min(dp[N - 1][2][0][0], dp[N - 1][2][1][0]))
    # index = [dp[N][1][0], dp[N][1][1], dp[N - 1][2][0], dp[N - 1][2][1]].index(max_value)
    # candi_path = [path[N][1][0], path[N][1][1], path[N - 1][2][0], path[N - 1][2][1]]
    # worst_path = candi_path[index]
    # \max(\min(dp[N][1][0][0], dp[N][1][1][0]), \min(dp[N][1][0][1], dp[N][1][1][1]),)
    worst_path = []

    return max_value, worst_path, dp, path


def generate_cuts(master_results, sub_results, added_cuts):
    if sub_results is not None:
        X, objVal = master_results
        max_value, worst_path, dp, path = sub_results
        # optimal_cuts
        added_cuts.append(['Optimal', max_value, [[key, X[key]] for key in X.keys()]])
        # translation cut todo: add多场桥cut
        # U_k = [[] for _ in range(data.K_num)]
        # for pair in X.items():
        #     U_k[data.J_K_dict[pair[1] + 1]].append([pair[0], pair[1]])
        # added_cuts.append(['Translation', max_value, [[key, X[key]] for key in X.keys()]])

    else:
        # todo:generate infeasible results
        a = 1
    return added_cuts


if __name__ == '__main__':
    data = read_data('/Users/jacq/PycharmProjects/BayAllocation/a_data_process/data/case1')
    CCG(data)
