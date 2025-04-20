import random


def generate_instance(container_setting, bay_setting, num_container_groups, num_block, num_occupied_bays, filename):
    """
    生成算例文件
    :param container_setting: 容器设置，'CS' 或 'CNS'
    :param bay_setting: 舱位设置，'BS' 或 'BNS'
    :param num_container_groups: 容器组数量
    :param num_occupied_bays: 被占用的舱位数（仅当 bay_setting 为 'BNS' 时有效）
    :param filename: 输出文件名
    """

    if container_setting == 'CS':
        container_groups = [12] * num_container_groups
    elif container_setting == 'CNS':
        container_groups = [random.randint(1, 5) * 12 for _ in range(num_container_groups)]
    else:
        raise ValueError("Invalid container setting. Must be 'CS' or 'CNS'.")

    bays = list(range(1, 61, 2))
    if bay_setting == 'BS':
        available_bays = [[bay for bay in bays] for _ in range(num_block)]
    elif bay_setting == 'BNS':
        if num_occupied_bays > len(bays):
            raise ValueError("Number of occupied bays cannot exceed total bays.")
        occupied_bays = [random.sample(bays, num_occupied_bays) for _ in range(num_block)]
        available_bays = [[bay for bay in bays if bay not in occupied_bays[k]] for k in range(num_block)]
    else:
        raise ValueError("Invalid bay setting. Must be 'BS' or 'BNS'.")

    # 写入文件
    with open(filename, 'w') as file:
        file.write("3 4\n")
        file.write(str(num_block) + "\n")
        for k in range(num_block):
            file.write(" ".join(map(str, available_bays[k])) + "\n")
        file.write(str(int(num_container_groups / 2)) + "\n")
        cnt = 0
        for _ in range(int(num_container_groups / 2)):
            file.write(str(container_groups[cnt]) + " ")
            file.write(str(container_groups[cnt + 1]) + "\n")
            cnt += 2


if __name__ == '__main__':
    # 示例：生成一个 CS-6 BNS-5 的算例
    # for
    C_type = 'CNS'
    B_type = 'BNS'  # 'BS'
    group_num = 6
    block_num = 2
    miss_bay_num = 5
    inst_ls = [[2, 1, 10], [4, 1, 10]]
    for ls in inst_ls:
        group_num, block_num, miss_bay_num = ls
        file_name = '/Users/jacq/PycharmProjects/BayAllocationGit/a_data_process/data/' + C_type + ' ' + B_type + '/' + \
                    C_type + '_' + str(group_num) + '_' + B_type + '_' + str(block_num) + '_' + str(
            miss_bay_num) + '.txt'
        generate_instance(C_type, B_type, group_num, block_num, miss_bay_num, file_name)
