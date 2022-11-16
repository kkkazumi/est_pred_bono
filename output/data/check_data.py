import numpy as np
from matplotlib import pyplot as plt

userlist = ["B9001","B9002","B9003"]
conlist = ["CN","ML"]

for username in userlist:
    for condition in conlist:
        filename="./"+username+"_"+condition+"/mental_last.csv"
        data=np.loadtxt(filename,delimiter=",")
        xdata_file="../../bono_datamaker/data/sudata/"+username+"_"+condition+"_factor.csv"
        _x = np.loadtxt(xdata_file,delimiter=",")
        x=_x[:,-1]*20+1
        plt.scatter(x,data,label=condition)
        plt.legend()
        print(condition,sum(data))
    plt.xlabel("the number of turns")
    plt.xlim(0,20)
    plt.xticks(range(21))
    plt.ylabel("estimated engagement of user")
    plt.show()
    
