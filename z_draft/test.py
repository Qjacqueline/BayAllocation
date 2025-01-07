# -*- coding: utf-8 -*-
# @Time    : 2024/12/27 12:55 PM
# @Author  : JacQ
# @File    : test.py
# import matplotlib
from matplotlib import pyplot as plt

x=[1, 2, 3, 4]
#y坐标轴上点的数值
y=[1, 4, 9, 16]
#第2步：使用plot绘制线条第1个参数是x的坐标值，第2个参数是y的坐标值
plt.plot(x,y)
#第3步：显示图形
plt.show()