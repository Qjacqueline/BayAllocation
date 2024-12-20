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

    G_num_set = []  # 箱组内集装箱个数 [20,40,20,40]
    G_sub_num_set = []  # 每个箱组的子箱组数group_sub_g_num
    U_L_num_set = []  # 20TEU子箱组集合记录对应数量 [12,12,12,3...]
    U_L_g_set = []  # 20TEU子箱组集合记录对应箱组 [1,1,1,1,...]
    U_F_num_set = []  # 40TEU子箱组集合记录对应数量 [12,12,12,3...]
    U_F_g_set = []  # 40TEU子箱组集合记录对应箱组 [1,1,1,1,...]
    U_num_set = []  # 所有子箱组集合记录对应数量 [12,12,12,3...]
    U_g_set = []  # 所有子箱组集合记录对应箱组 [1,1,1,1,...]

    J = []  # 总贝集合 [3,5,11,...]
    I = []  # 存在后继贝位的贝位编号 [3]
    J_K_first = []  # 每个箱区第一个空贝位置
    J_K = []  # 每个block内bay标号

    U = []  # 子箱组标号集合
    U_L = []  # 20TEU子箱组标号集合
    U_F_u = []  # 40TEU子箱组前标号集合
    U_F_l = []  # 40TEU子箱组后标号集合
    U_F = []  # 40TEU子箱组标号集合
    p_u = []  # 子箱组的优先级
    t_u = []  # 子箱组操作所需时间
    gamma_uu = []  # 箱组间混贝成本

    def init_process(self):
        for q in range(len(self.G_num_set)):
            if q % 2 == 0:
                c_num = self.G_num_set[q]
                tmp = math.ceil(c_num / (self.S_num * self.T_num))
                self.G_sub_num_set.append(tmp)
                for qq in range(tmp - 1):
                    self.U_L_num_set.append(self.S_num * self.T_num)
                    self.U_L_g_set.append(q)
                self.U_L_num_set.append(c_num - (self.S_num * self.T_num) * (tmp - 1))
                self.U_L_g_set.append(q)
            else:
                c_num = self.G_num_set[q]
                tmp = math.ceil(c_num / (self.S_num * self.T_num))
                self.G_sub_num_set.append(tmp * 2)
                for qq in range(tmp - 1):
                    self.U_F_num_set.append(self.S_num * self.T_num)
                    self.U_F_g_set.append(q)
                    self.U_F_num_set.append(self.S_num * self.T_num)
                    self.U_F_g_set.append(q)
                self.U_F_num_set.append(int(c_num - (self.S_num * self.T_num) * (tmp - 1)))
                self.U_F_g_set.append(q)
                self.U_F_num_set.append(int(c_num - (self.S_num * self.T_num) * (tmp - 1)))
                self.U_F_g_set.append(q)
        self.U_L_num = len(self.U_L_num_set)
        self.U_F_num = len(self.U_F_num_set)
        self.U_num = self.U_L_num + self.U_F_num
        self.t_u = [c_num * cf.unit_process_time for c_num in self.U_L_num_set]
        self.gamma_uu = [[1 for _ in range(self.U_num)] for _ in range(self.U_num)]
        self.U_num_set.extend(self.U_L_num_set)
        self.U_num_set.extend(self.U_F_num_set)
        self.U_g_set.extend(self.U_L_g_set)
        self.U_g_set.extend(self.U_F_g_set)
        self.U = [i for i in range(self.U_num)]  # 子箱组标号集合
        self.U_L = [i for i in range(self.U_L_num)]  # 20TEU子箱组标号集合
        self.U_F_u = [i for i in range(self.U_L_num, self.U_num, 2)]  # 40TEU子箱组前标号集合
        self.U_F_l = [i for i in range(self.U_L_num + 1, self.U_num, 2)]  # 40TEU子箱组后标号集合
        self.U_F = [i for i in range(self.U_L_num, self.U_num)]


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
        bays = [j + 60 * k for j in bays]
        data.J_K_first.append(bays[0])
        data.J.extend(bays)
        data.J_K.append(bays)
        for j in bays:
            if j + 2 in bays:
                data.I.append(j)
    data.G_num_set = []
    data.D_num = int(lines[data.K_num + 2][0])
    data.G_num = data.D_num * 2
    for k in range(data.D_num):
        groups = list(map(int, lines[data.K_num + k + 3].split(" ")))
        data.G_num_set.extend(groups)
    data.init_process()
    return data


if __name__ == '__main__':
    data = read_data('data/case1')
