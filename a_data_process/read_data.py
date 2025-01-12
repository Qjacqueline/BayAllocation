# -*- coding: utf-8 -*-
# @Time    : 2022/6/11 11:30 AM
# @Author  : JacQ
# @File    : read_data.py
import math
from a_data_process import config as cf


class Data:
    S_num = 0  # 堆垛数
    T_num = 0  # 层高
    K_num = 0  # 箱区数
    D_num = 0  # 卸货港数
    G_num = 0  # 箱组数需要
    J_num = 0  # 空贝位数
    U_num = 0  # 子箱组数
    U_L_num = 0  # 20TEU子箱组数
    U_F_num = 0  # 40TEU子箱组数，含double虚拟

    G_num_set = []  # 箱组内集装箱个数 [14，13，12，12]
    G_sub_num_set = []  # 每个箱组的子箱组数group_sub_g_num
    # G_sub_index_set = []  # 每个箱组内子箱组标号集合[]
    U_L_num_set = []  # 20TEU子箱组集合记录对应数量             [12, 2, 12]
    U_L_g_set = []  # 20TEU子箱组集合记录对应箱组                 [0, 0, 2]
    U_F_num_set = []  # 40TEU子箱组集合记录对应数量             [12, 1, 12]
    U_F_g_set = []  # 40TEU子箱组集合记录对应箱组                 [1, 1, 3]
    U_num_set = []  # 所有子箱组集合记录对应数量
    U_g_set = []  # 所有子箱组集合记录对应箱组

    J = []  # 总贝集合 [3,5,11,...]
    I = []  # 存在后继贝位的贝位编号 [3]
    J_K_first = []  # 每个箱区第一个空贝位置
    J_K = []  # 每个block内bay标号
    J_K_dict = {}  # 每个block内bay标号

    U = []  # 子箱组标号集合
    U_L = []  # 20TEU子箱组标号集合
    U_F = []  # 40TEU子箱组标号集合
    U_g = []  # g箱组子箱组集合
    p_u = []  # 子箱组的优先级
    # t_u = []  # 子箱组操作所需时间
    gamma_uu = []  # 箱组间混贝成本

    def init_process(self):
        # todo 处理个数为0的情况
        u = 0
        for q in range(len(self.G_num_set)):
            c_num = self.G_num_set[q]
            if c_num == 0:
                continue
            if q % 2 == 0:
                tmp = math.ceil(c_num / (self.S_num * self.T_num))
                self.G_sub_num_set.append(tmp)
                for qq in range(tmp - 1):
                    self.U_L_num_set.append(self.S_num * self.T_num)
                    self.U_num_set.append(self.S_num * self.T_num)
                    self.U_L_g_set.append(q)
                    self.U_g_set.append(q)
                    self.U_L.append(u)
                    u = u + 1
                if c_num == 0:
                    self.U_L_num_set.append(0)
                    self.U_num_set.append(0)
                else:
                    self.U_L_num_set.append(c_num - (self.S_num * self.T_num) * (tmp - 1))
                    self.U_num_set.append(c_num - (self.S_num * self.T_num) * (tmp - 1))
                self.U_L_g_set.append(q)
                self.U_g_set.append(q)
                self.U_L.append(u)
                u = u + 1
            else:
                tmp = math.ceil(c_num / (self.S_num * self.T_num))
                self.G_sub_num_set.append(tmp)
                for qq in range(tmp - 1):
                    self.U_F_num_set.append(self.S_num * self.T_num)
                    self.U_num_set.append(self.S_num * self.T_num)
                    self.U_F_g_set.append(q)
                    self.U_g_set.append(q)
                    self.U_F.append(u)
                    u = u + 1
                if c_num == 0:
                    self.U_F_num_set.append(0)
                    self.U_num_set.append(0)
                else:
                    self.U_F_num_set.append(int(c_num - (self.S_num * self.T_num) * (tmp - 1)))
                    self.U_num_set.append(int(c_num - (self.S_num * self.T_num) * (tmp - 1)))
                self.U_F_g_set.append(q)
                self.U_g_set.append(q)
                self.U_F.append(u)
                u = u + 1
        self.G_num = len(self.U_g_set)
        self.U_L_num = len(self.U_L_num_set)
        self.U_F_num = len(self.U_F_num_set)
        self.U_num = self.U_L_num + self.U_F_num
        self.U = [i for i in range(self.U_num)]  # 子箱组标号集合
        self.U_g = [[] for g in range(self.G_num)]
        for i in range(self.U_num):
            self.U_g[self.U_g_set[i]].append(i)


def read_data(path: str):
    f = open(path, 'r')
    lines = f.readlines()
    data = Data()
    # read the info
    data.S_num, data.T_num = list(map(int, lines[0].split(" ")))
    data.K_num = int(lines[1][0])
    data.J = []
    for k in range(data.K_num):
        bays = list(map(int, lines[2 + k].split(" ")))
        bays = [j + cf.bay_number_one_block * k for j in bays]
        data.J_K_first.append(bays[0])
        data.J.extend(bays)
        data.J_K.append(bays)
        for b in bays:
            data.J_K_dict[b] = k
        for j in bays:
            if j + 2 in bays:
                data.I.append(j)
    data.G_num_set = []
    data.D_num = int(lines[data.K_num + 2][0])
    for k in range(data.D_num):
        groups = list(map(int, lines[data.K_num + k + 3].split(" ")))
        data.G_num_set.extend(groups)
    data.init_process()
    return data


if __name__ == '__main__':
    data = read_data('/Users/jacq/PycharmProjects/BayAllocation/a_data_process/data/case4')
