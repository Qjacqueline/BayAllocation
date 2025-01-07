# -*- coding: utf-8 -*-
# @Time    : 2022/6/10 15:30 AM
# @Author  : main.py

from z_draft.model_V1 import model
from a_data_process.read_data import read_data

if __name__ == '__main__':

    # input data
    path = 'data_process/data/case1'
    data = read_data(path)

    # Create a new model
    m = model(data)
