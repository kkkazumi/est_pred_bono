import numpy as np
import pandas as pd

#user_list={'B9001','B9002','B9003','B9004','B9005'}
user_list={'B9001','B9002','B9003'}
cond_list={"ML","CN"}

ROWS=14

def divide_data(data):
    new_data=np.copy(data)
    new_data[:,8]=data[:,8]/20.0
    return new_data

def extract_first_pause(df,target_label):
    return df[df[target_label]==0]

def drop_data(df,drop_list):
    return df.drop(drop_list,axis=1)

def switch_data(df,new_columns):
    return df.reindex(columns=new_columns)

def get_factor_data(userid,condition):
    filename = "./data/"+condition+"/"+userid+"_sudata_log.csv"
    data = np.loadtxt(filename,delimiter=",")
    new_data = divide_data(data)
    columns=['oqtype','topic','pausenum','val1','val2','val3','pause_dur','val4','turn','day','hour','min','sec','msec']
    df = pd.DataFrame(data=new_data,columns=columns,dtype='float')
    df['turn num'] = data[:,8]+1
    target_label='pausenum'
    drop_list=['pausenum','val1','val2','val3','val4','day','hour','min','sec','msec']
    new_columns=['pause_dur','oqtype','topic','turn','turn num']
    _new_df = drop_data(extract_first_pause(df,target_label),drop_list)
    new_df = switch_data(_new_df,new_columns)
    new=data[:,8] + 1
    return new_df

def get_signal_data(userid,condition):
    filename = "./data/sound_emo/"+userid+"_"+condition+"_emotion.csv"
    data=pd.read_csv(filename)
    drop_list = ['calm','energy','sorrow']
    new_df=drop_data(data,drop_list)
    _data=new_df.groupby('turn num').mean()
    _data.reset_index(inplace=True)
    return _data

for userid in user_list:
    for condition in cond_list:
        factor=get_factor_data(userid,condition)
        signal=get_signal_data(userid,condition)
        TURN_NUM = 20
        new_factor= np.zeros((1,5))
        new_signal= np.zeros((1,3))
        count = 0
        for t in range(TURN_NUM):
            if(any(factor['turn num']==t)):
                if(any(signal['turn num']==t)):
                    fval=factor[factor['turn num']==t].values
                    sval=signal[signal['turn num']==t].values/100.0
                    if(count==0):
                        new_factor[0,:]=fval[0]
                        new_signal[0,:]=sval[0]
                        count+=1
                    else:
                        new_factor=np.append(new_factor,fval,axis=0)
                        new_signal=np.append(new_signal,sval,axis=0)
        fsave_arr=new_factor[:,:4]
        ssave_arr=new_signal[:,1:]

        ffilename = "./data/sudata/"+userid+"_"+condition+"_factor.csv"
        sfilename = "./data/sudata/"+userid+"_"+condition+"_signal.csv"
        np.savetxt(ffilename,fsave_arr,delimiter=",",fmt='%1.3f')
        np.savetxt(sfilename,ssave_arr,delimiter=",",fmt='%1.3f')
