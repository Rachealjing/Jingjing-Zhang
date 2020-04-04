import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta

# read data
clean_0012 = pd.read_csv('clean_0012.csv')
clean_0012['PERMNO'] = clean_0012['PERMNO'].apply(str)
clean_0012['DATE'] = pd.to_datetime(clean_0012['DATE'])
clean_0012['NAMEENDT'] = pd.to_datetime(clean_0012['NAMEENDT'])

datedf = pd.DataFrame()
RV_list = pd.DataFrame()
list1 = [np.nan]

k = 0
#path = r'D:\jing\intraday\newintraday_31.csv'
path = r'E:\thesis\data\merge\intraday_50.csv'
for df in pd.read_csv(path, chunksize=1000000):
    k = k + 1
    print(df)
    df['DATE'] = pd.to_datetime(df['DATE'], format='%Y%m%d')
    df['TIME'] = pd.to_timedelta(df['TIME'])
    df['DATETIME'] = df['DATE'] + df['TIME']
    df['DATE'] = df['DATE'].apply(lambda x: datetime.date(x))
    datedf = datedf.append(pd.DataFrame(df['DATE']), ignore_index=True)
    df['DATE'] = df['DATE'].apply(str)
    datelist = df['DATE'].unique().tolist()
    df.set_index(df['DATETIME'], inplace=True)
    df.drop(columns=['SIZE', 'G127', 'CORR', 'COND', 'EX', 'TSEQ'], inplace=True)
    if k == 1:
        temdate = datelist[-1]
        tem = df[temdate]
    else:
        temdate_n = datelist[-1]
        df = pd.concat([tem, df], axis=0)
        if datelist[0] != temdate:
            datelist.insert(0, temdate)
    for b in range(0, len(datelist)-1):
        df_clean = pd.DataFrame()
        listime = []
        c = datelist[b]
        current = df[df['DATE'] == c]
        dt1 = datetime(int(c[0:4]), int(c[5:7]), int(c[8:10]), 9, 30, 0)
        dt2 = datetime(int(c[0:4]), int(c[5:7]), int(c[8:10]), 16, 0, 0)
        current = current[dt1:dt2]
        dt = dt1
        while dt <= dt2:
            listime.append(dt)
            dtt = dt + timedelta(minutes=5)
            valid = current[dt:dtt]
            if valid.empty == True:
                RV_5m = pd.DataFrame(list1, columns=['PRICE'])
                df_clean = df_clean.append(RV_5m, ignore_index=True)
                dt = dt + timedelta(minutes=5)
            else:
                RV_5m = (pd.DataFrame(valid.mean())).T
                df_clean = df_clean.append(RV_5m, ignore_index=True)
                dt = dt + timedelta(minutes=5)
        df_clean.index = listime
        df_clean = df_clean.fillna(axis=0, method='ffill')
        df_clean = df_clean.fillna(axis=0, method='bfill')
        dfmean = []
        for a in range(0, 3):
            listmean = pd.DataFrame()
            dtm = dt1 + a * timedelta(minutes=5)
            while dtm < dt2:
                listmean = listmean.append(df_clean[dtm:dtm], ignore_index=False)
                dtm = dtm + timedelta(minutes=15)
            listmean['return'] = listmean['PRICE'].pct_change()
            listmean['R_square'] = listmean['return'].map(lambda x: x ** 2)
            RV15Min = listmean['R_square'].sum()
            dfmean.append(RV15Min)
        RV = (pd.DataFrame(dfmean)).mean()
        RV_list = RV_list.append(RV, ignore_index=True)
    if k > 1:
        temdate = temdate_n
        tem = df[temdate_n]
    if temdate == '2012-12-31':
        df_clean = pd.DataFrame()
        listime = []
        c = datelist[-1]
        current = df[df['DATE'] == c]
        dt1 = datetime(int(c[0:4]), int(c[5:7]), int(c[8:10]), 9, 30, 0)
        dt2 = datetime(int(c[0:4]), int(c[5:7]), int(c[8:10]), 16, 0, 0)
        current = current[dt1:dt2]
        dt = dt1
        while dt <= dt2:
            listime.append(dt)
            dtt = dt + timedelta(minutes=5)
            valid = current[dt:dtt]
            if valid.empty == True:
                RV_5m = pd.DataFrame(list1, columns=['PRICE'])
                df_clean = df_clean.append(RV_5m, ignore_index=True)
                dt = dt + timedelta(minutes=5)
            else:
                RV_5m = (pd.DataFrame(valid.mean())).T
                df_clean = df_clean.append(RV_5m, ignore_index=True)
                dt = dt + timedelta(minutes=5)
        df_clean.index = listime
        df_clean = df_clean.fillna(axis=0, method='ffill')
        df_clean = df_clean.fillna(axis=0, method='bfill')
        dfmean = []
        for a in range(0, 3):
            listmean = pd.DataFrame()
            dtm = dt1 + a * timedelta(minutes=5)
            while dtm < dt2:
                listmean = listmean.append(df_clean[dtm:dtm], ignore_index=False)
                dtm = dtm + timedelta(minutes=15)
            listmean['return'] = listmean['PRICE'].pct_change()
            listmean['R_square'] = listmean['return'].map(lambda x: x ** 2)
            RV15Min = listmean['R_square'].sum()
            dfmean.append(RV15Min)
        RV = (pd.DataFrame(dfmean)).mean()
        RV_list = RV_list.append(RV, ignore_index=True)
datep = datedf['DATE'].unique().tolist()
tickerlist = df['SYMBOL'].unique().tolist()
perline = clean_0012[clean_0012['TICKER'] == tickerlist[0]]
per = perline['PERMNO'].unique().tolist()
RV_list.index = datep
RV_list.columns = [per[0]]

# read df
df_0012_RVd = pd.read_csv("df_0012_RVd.csv")
df_0012_RVd['caldt'] = pd.to_datetime(df_0012_RVd['caldt'])
df_0012_RVd['caldt'] = df_0012_RVd['caldt'].apply(lambda x: datetime.date(x))
df_0012_RVd.set_index("caldt", inplace=True)

# write df
#df_0012_RVd[per[0]] = RV_list[per[0]]
#df_0012_RVd.to_csv("df_0012_RVd.csv", index = True)



