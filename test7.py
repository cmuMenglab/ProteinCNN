import numpy as np
import pandas as pd

# 将输入的氨基酸序列转化为相对应的理化属性
def ammino_order_to_physi_property(ammino_order_alist):
    file_path = 'abc.xltx'
    df_data = pd.DataFrame(pd.read_excel(file_path))
    physi_property = []
    for item in ammino_order_alist:
        physi_property.append(df_data.loc[item, :].values)
    physi_property = np.array(physi_property)
    return physi_property


# 实现通过lag和j求第j个理化属性
def AC_from_lag_and_j(P, j, L, lag):
    weight1 = (1.0 / (L - lag))
    weight2 = (1.0 / L)

    # Sigma_Pij
    Sigma_Pij = 0
    for i in range(L):
        Sigma_Pij += P[i][j]
    # Sigma_AC
    Sigma_AC = 0
    for i in range(L - lag):
        Sigma_AC += (P[i][j] - weight2 * Sigma_Pij) * (P[i + lag][j] - weight2 * Sigma_Pij)
    return weight1 * Sigma_AC



AC_list = []
# 求shape为lag x 7 的特征向量
def AC_eigenvector(P, L, lag):

    for j in range(7):
        AC_list.append(AC_from_lag_and_j(P, j, L, lag))
    ac_eigenvector = np.array(AC_list)
    return ac_eigenvector

if __name__ == '__main__':
    # ammino_order = 'QRRQRRGEERKAPENQEEEEERAELNQSEEPEAGESSTGGP'
    ammino_order = input('请输入氨基酸序列：')
    ammino_order_alist = list(ammino_order)
    L = len(ammino_order_alist)
    lg = int(input('请输入lag的值：'))
    lags = [i + 1 for i in range(lg)]

    # j = int(input("请输入j的值："))
    physi_property = ammino_order_to_physi_property(ammino_order_alist)    # 获取序列的理化属性值
    # print(physi_property)
    for lag in lags:
        ac_eigenvector_from_lag = AC_eigenvector(physi_property, L, lag)
    print(ac_eigenvector_from_lag)





# 导出结果
data = pd.DataFrame(ac_eigenvector_from_lag,)
# print(data)

# writer = pd.ExcelWriter('output.xlsx')       #导出到excel
# data.to_excel(writer, float_format = '%.5f')
# writer.save()

data.to_csv("125.csv", header = None, index=None)   #导出到csv
