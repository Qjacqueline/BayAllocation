# -*- coding: utf-8 -*-
# @Time    : 2024/11/11 7:06 PM
# @Author  : JacQ
# @File    : SP_YC.py

import gurobipy as gp
from gurobipy import GRB


def SP_gurobi():
    # Initialize Gurobi model
    model = gp.Model("SPR")

    # Constants
    G_num = 5  # Replace this with the actual size of G (number of layers)
    N_num = 3  # Number of nodes, including 0, 1, and 2

    c = {}
    for p in range(G_num + 1):
        for l in range(0, N_num):
            for ll in range(0, N_num):
                if l == 1 and ll == 1:
                    c[p, l, ll] = 4
                elif ll == 0:
                    c[p, l, ll] = 0
                else:
                    c[p, l, ll] = 3

    lambd = [[[model.addVar(vtype=GRB.CONTINUOUS, name=f'lambd_{p}_{l}_{ll}')
               for ll in range(N_num)] for l in range(N_num)] for p in range(G_num + 1)]
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='obj')

    model.addConstr(obj == gp.quicksum(c[p, l, ll] * lambd[p][l][ll]
                                       for p in range(G_num + 1)
                                       for l in range(0, N_num)
                                       for ll in range(0, N_num)), name="objective")
    model.setObjective(obj, GRB.MINIMIZE)

    model.addConstr(lambd[0][0][1] + lambd[0][0][2] == 1, name="4b")
    for l in range(1, N_num):
        model.addConstr(lambd[0][0][l] - (lambd[1][l][1] + lambd[1][l][2]) == 0, name=f"4c_l{l}")
    for p in range(2, G_num):
        for l in range(1, N_num):
            model.addConstr((lambd[p - 1][1][l] + lambd[p - 1][2][l]) -
                            (lambd[p][l][1] + lambd[p][l][2]) == 0, name=f"4d_p{p}_l{l}")
    for l in range(1, N_num):
        model.addConstr((lambd[G_num - 1][1][l] + lambd[G_num - 1][2][l]) -
                        lambd[G_num][l][0] == 0, name=f"4e_l{l}")
    model.addConstr(lambd[G_num][1][0] + lambd[G_num][2][0] == 1, name="4f")

    model.optimize()

    # Output solution
    if model.status == GRB.OPTIMAL:
        print("Optimal solution found:")
        for p in range(G_num + 1):
            for l in range(N_num):
                for ll in range(N_num):
                    if lambd[p][l][ll].x > 0.5:
                        print(f"lambd[{p}][{l}][{ll}] = 1")
        print(f"Objective value: {obj.x}")
    else:
        print("No optimal solution found.")

    # Initialize Gurobi model
    model_d = gp.Model("SPRD")

    # Variables
    alpha = model_d.addVar(vtype=GRB.CONTINUOUS, name="alpha")
    beta = model_d.addVar(vtype=GRB.CONTINUOUS, name="beta")
    pi = [[model_d.addVar(vtype=GRB.CONTINUOUS, name=f"pi_{p}_{l}") for l in range(0, 3)] for p in range(G_num + 1)]

    # Objective function
    model_d.setObjective(alpha + beta, GRB.MAXIMIZE)

    # Constraints
    # Constraint (5b)
    for l in range(1, 3):
        model_d.addConstr(alpha + pi[1][l] <= c[0, 0, l], name=f"5b_l{l}")

    # Constraint (5c)
    for p in range(1, G_num):
        for l in range(1, 3):
            for ll in range(1, 3):
                model_d.addConstr(pi[p + 1][ll] - pi[p][l] <= c[p, l, ll],
                                  name=f"5c_p{p}_l{l}_lprime{ll}")

    # Constraint (5d)
    for l in range(1, 3):
        model_d.addConstr(beta - pi[G_num][l] <= c[G_num, l, 0], name=f"5d_l{l}")

    # Optimize the model
    model_d.optimize()

    # Output solution
    if model_d.status == GRB.OPTIMAL:
        print("Optimal solution found:")
        print(f"alpha = {alpha.x}")
        print(f"beta = {beta.x}")
        for p in range(G_num + 1):
            for l in range(1, 3):
                print(f"pi[{p}][{l}] = {pi[p][l - 1].x}")
        print(f"Objective value: {model_d.objVal}")
    else:
        print("Model is infeasible. Computing IIS...")
        model_d.computeIIS()
        model_d.write("model_d.ilp")  # Save the IIS to a file for review
        print("Irreducible Infeasible Set (IIS) written to 'model_d.ilp'")
