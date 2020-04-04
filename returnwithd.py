import pandas as pd
import numpy as np
from datetime import datetime


# read data
origin_d = pd.read_csv('dailyclean.csv')
origin_d['PERMNO'] = origin_d['PERMNO'].apply(str)
origin_d['date'] = pd.to_datetime(origin_d['date'], format='%Y%m%d')
origin_d['RET'].replace('C', 0.0, inplace=True)
origin_d['RET'].replace('B', 0.0, inplace=True)
origin_d['RET'] = origin_d['RET'].apply(float)
origin_d['DLRET'].replace('A', 0.0, inplace=True)
origin_d['DLRET'].replace('S', 0.0, inplace=True)
origin_d['DLRET'].replace('T', 0.0, inplace=True)
origin_d['DLRET'].replace('P', 0.0, inplace=True)
origin_d['DLRET'] = origin_d['DLRET'].apply(float)
origin_d['RET'].replace(np.nan, 0.0, inplace=True)
origin_d['DLRET'].replace(np.nan, 0.0, inplace=True)
origin_d['retur'] = (1 + origin_d.RET) * (1 + origin_d.DLRET) - 1
origin_d['price'] = np.abs(origin_d.PRC) / origin_d.CFACPR
origin_d['cap'] = origin_d.price * (origin_d.SHROUT * origin_d.CFACSHR) * 1000
ignore = origin_d[origin_d['DLRET'] == -1.0]


data = pd.read_csv('sp500_8618.csv')
dateindex = data['caldt']
print(dateindex)
data['caldt'] = pd.to_datetime(data['caldt'], format='%Y%m%d')

col = origin_d['PERMNO'].unique().tolist()
date_port = pd.date_range(start="19870228", end="20181231", freq="M")
date_port = pd.to_datetime(date_port).strftime('%Y%m%d')


# build return and cap table
df_price = pd.DataFrame(index = data['caldt'])
df_return = pd.DataFrame(index = data['caldt'])
df_cap = pd.DataFrame(index = data['caldt'])

for a in range(0, len(col)):
    b = col[a]
    index_col = origin_d[(origin_d.PERMNO == b)].index.tolist()
    start = min(index_col)
    end = max(index_col)
    temdata = origin_d[start:end+1]
    temdata.set_index("date", inplace=True)
    df_return[b] = temdata.retur
    df_cap[b] = temdata.cap
    df_price[b] = temdata.price
a = 0
df_return.to_csv("df_return.csv", index = True)
df_cap.to_csv("df_cap.csv", index = True)
df_price.to_csv("df_price.csv", index = True)


# delete stocks in portfolio
df_portf = pd.DataFrame()

for a in range(0,date_port.__len__()):
    df_ret_test = pd.DataFrame(df_return)
    b = date_port[a]
    year = b[0:4]
    month = int(b[4:6])
    # delete t-13 price = nan or 0
    if month < 2:
        t_13 = str(int(year) - 2) + '-' + str(month + 11)
    else:
        t_13 = str(int(year) - 1) + '-' + str(month - 1)
    tem_df_price = df_price[t_13]
    nap = tem_df_price.isnull().all()
    nap = nap.to_frame()
    nap_ind = nap[(nap[0] == True)].index.tolist()
    print(str(b)+':', len(nap_ind))
    df_ret_test.drop(columns=nap_ind, axis=1, inplace=True)
    # delete t-2 return = nan
    if month < 3:
        t_2 = str(int(year) - 1) + '-' + str(month + 10)
    else:
        t_2 = str(int(year)) + '-' + str(month - 2)
    tem_df_return = df_ret_test[t_2]
    nar = tem_df_return.isnull().all()
    nar = nar.to_frame()
    nar_ind = nar[(nar[0] == True)].index.tolist()
    print(str(b) + ':', len(nar_ind))
    df_ret_test.drop(columns=nar_ind, axis=1, inplace=True)
    testlist = (pd.DataFrame(df_ret_test.columns)).T
    df_portf = df_portf.append(testlist.iloc[0], ignore_index=True)
a = 0
df_portf.index = date_port
df_portf.to_csv("df_portf.csv", index = True)


# check the weight
nawtest = pd.DataFrame()

for a in range(0,date_port.__len__()):
    df_wei_test = pd.DataFrame(df_cap)
    b = date_port[a]
    year = b[0:4]
    month = int(b[4:6])
    # delete t-1 weight = nan
    if month < 2:
        t_1 = str(int(year) - 1) + '-' + str(month + 11)
    else:
        t_1 = str(int(year)) + '-' + str(month - 1)
    tem_df_cap = df_cap[t_1]
    naw = tem_df_cap.isnull().all()
    naw = naw.to_frame()
    naw_ind = (pd.DataFrame(naw[(naw[0] == False)].index)).T
    nawtest = nawtest.append(naw_ind, ignore_index=True)
a = 0
nawtest.index = date_port
nawtest.to_csv("nawtest.csv", index = True)


# build portfolio
df_port = pd.DataFrame()
for a in range(0,date_port.__len__()):
    b = date_port[a]
    year = b[0:4]
    month = int(b[4:6])
    list1 = df_portf.iloc[a].tolist()
    list2 = nawtest.iloc[a].tolist()
    list12 = [l for l in list1 if l in list2]
    clist = [x for x in list12 if str(x) != 'nan']
    print(len(clist))
    port_list = (pd.DataFrame(clist)).T
    df_port = df_port.append(port_list, ignore_index=True)
a = 0
df_port.index = date_port
df_port.to_csv("df_port.csv", index = True)



#calculate cumulative return for each month
df_cum = pd.DataFrame()
for a in range(0,date_port.__len__()):
    b = date_port[a]
    year = b[0:4]
    month = int(b[4:6])
    start_port = str(int(year) - 1) + '-' + str(month)
    if month < 3:
        end_port = str(int(year) - 1) + '-' + str(month + 10)
    else:
        end_port = str(year) + '-' + str(month - 2)
    tem_df_cure = df_return[start_port : end_port]
    cum_return = np.cumprod(1 + tem_df_cure, axis=0)[-1:] - 1
    df_cum = df_cum.append(cum_return, ignore_index=True)
a = 0
df_cum.index = date_port
df_cum.to_csv("df_cum.csv", index = True)



# build 10 decile portfolio
date = dateindex[274:8318].tolist()

names = locals()
for i in range(1, 11):
    names['port' + str(i)] = pd.DataFrame()
i = 0
for a in range(0, date_port.__len__()):
    b = date_port[a]
    dftest = df_port.iloc[a].tolist()
    dftest = [x for x in dftest if str(x) != 'nan']
    ptest = pd.DataFrame(df_cum, columns=dftest)
    port = pd.DataFrame(ptest.iloc[a])
    sort_port = port.sort_values(axis=0, by=b, ascending=True)
    for i in range(1, 11):
        c = round(len(sort_port) * (i - 1) / 10)
        d = round(len(sort_port) * i / 10)
        porttest = sort_port[c:d].T
        portl = porttest.columns.values.tolist()
        names['port' + str(i)] = names['port' + str(i)].append((pd.DataFrame(portl)).T)
a = 0
i = 0
for i in range(1, 11):
    names['port' + str(i)].index = date_port
i = 0

############################################################################ trial 1
# calculate cap-weighted return
#for i in range(1, 11):
#    g = 0
#    e = 0
#    names['ret_port' + str(i)] = pd.DataFrame()
#    month_b = int(date_port[g][0:6])
#    f = str(date[e])
#   month_e = int(f[0:6])
#    while g < 383:
#        while month_b == month_e:
#            col_list = names['port' + str(i)].iloc[g].tolist()
#            col_list = [x for x in col_list if str(x) != 'nan']
#            tem_ret = pd.DataFrame(df_return[f: f], columns=col_list)
#            tem_cap = pd.DataFrame(df_cap[f: f], columns=col_list)
#            tem_weight = tem_cap.apply(lambda x: x / x.sum(), axis = 1)
#            tem_wei_ret = tem_ret * tem_weight
#            wei_ret = tem_wei_ret.sum(axis=1)
#            names['ret_port' + str(i)] = names['ret_port' + str(i)].append(pd.DataFrame(wei_ret))
#            if e < 8043:
#                e = e + 1
#                f = str(date[e])
#                month_e = int(f[0:6])
#            else:
#                e = 0
#                f = str(date[e])
#                month_e = int(f[0:6])
#        if g < 382:
#            g = g + 1
#            month_b = int(date_port[g][0:6])
#        else:
#            g = 383
#    names['ret_port' + str(i)].index = data['caldt'][274:8318]
#i = 0

# calculate cumulative return
#for i in range(1, 11):
#    sta = '19870202'
#    names['cumret_port' + str(i)] = pd.DataFrame()
#    for a in range(0, len(date)):
#        b = str(date[a])
#        tem_cum = names['ret_port' + str(i)][sta: b]
#        cum_re = np.cumprod(1 + tem_cum, axis=0)[-1:]
#        names['cumret_port' + str(i)] = names['cumret_port' + str(i)].append(pd.DataFrame(cum_re))
#    names['cumret_port' + str(i)].index = data['caldt'][274:8318]
#    a = 0
#i = 0


############################################################################ trial 2
date_port1 = pd.date_range(start="19870131", end="20181130", freq="M")
date_port1 = pd.to_datetime(date_port1).strftime('%Y-%m-%d')
# calculate cap-weighted return
for i in range(1, 11):
    g = 0
    e = 0
    names['ret_port' + str(i)] = pd.DataFrame()
    month_b = int(date_port[g][0:6])
    f = str(date[e])
    month_e = int(f[0:6])
    while g < 383:
        while month_b == month_e:
            col_list = names['port' + str(i)].iloc[g].tolist()
            col_list = [x for x in col_list if str(x) != 'nan']
            tem_ret = pd.DataFrame(df_return[f: f], columns=col_list)
            tem_ret = np.array(tem_ret)
            h = str(date_port1[g][0:7])
            tem_cap = pd.DataFrame(df_cap[h], columns=col_list)
            tem_weight = tem_cap[-1:].apply(lambda x: x / x.sum(), axis = 1)
            tem_weight = np.array(tem_weight)
            tem_wei_ret = pd.DataFrame(np.multiply(tem_ret, tem_weight)).sum(axis=1)
            names['ret_port' + str(i)] = names['ret_port' + str(i)].append(pd.DataFrame(tem_wei_ret))
            if e < 8043:
                e = e + 1
                f = str(date[e])
                month_e = int(f[0:6])
            else:
                e = 0
                f = str(date[e])
                month_e = int(f[0:6])
        if g < 382:
            g = g + 1
            month_b = int(date_port[g][0:6])
        else:
            g = 383
    names['ret_port' + str(i)].index = data['caldt'][274:8318]
i = 0

# calculate cumulative return
for i in range(1, 11):
    sta = '19870202'
    names['cumret_port' + str(i)] = pd.DataFrame()
    for a in range(0, len(date)):
        b = str(date[a])
        tem_cum = names['ret_port' + str(i)][sta: b]
        cum_re = np.cumprod(1 + tem_cum, axis=0)[-1:]
        names['cumret_port' + str(i)] = names['cumret_port' + str(i)].append(pd.DataFrame(cum_re))
    names['cumret_port' + str(i)].index = data['caldt'][274:8318]
    a = 0
i = 0


# connect CRSP with TAQ and Option Metrics
listper = []
for a in range(0,date_port.__len__()):
    orlist = df_port.iloc[a].tolist()
    clist = [x for x in orlist if str(x) != 'nan']
    listper.extend(clist)
    listper = list(set(listper))
a = 0
file = pd.DataFrame(listper)
file.columns = ['PERMNO']
file.to_csv('file.txt', sep='\t', index=False)

# find the stocks that appear most times in port10
listper10 = []
for a in range(0,date_port.__len__()):
    orlist = port10.iloc[a].tolist()
    clist = [x for x in orlist if str(x) != 'nan']
    listper10.extend(clist)
a = 0
result10 = pd.value_counts(listper10)