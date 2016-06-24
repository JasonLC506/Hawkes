import math
THRESHOLD = [0.01 for i in range(10)]


def Hawkes_predict(series, tr, tp):
    if len(series)<= tp+1:
        print "too short series"
        return [-1,-1],[-1,-1]
    paras = Hawkes_fit_EM(series[:(tr+1)])
    # according to TingTing He et al. paper
    predict = [0 for n in range(2)]
    for n in range(2):
        predict[n] = 1 - pow((tp-tr),-paras[n]) + paras[2*(n+1)]*(1-math.exp(-pow(tp-tr,2)/2/pow(paras[n*2+6],2)))\
        + paras[5-n*2]*(1-math.exp(-pow(tp-tr,2)/2/pow(paras[9-n*2],2)))
        predict[n] += series[tr][n]
    print series[tp][0], predict[0]
    print series[tp][1], predict[1]
    abs_err = [abs(predict[i]-series[tp][i]) for i in range(2)]
    rel_err = [0 for i in range(2)]
    for i in range(2):
        if series[tp][i]>0:
            rel_err[i] = abs_err[i]/series[tp][i]
        else:
            rel_err[i] = abs_err[i]
    return abs_err, rel_err

def Hawkes_fit_EM(series):

    # para_list
    # alpha1, alpha2, yita11, yita12, yita22, yita21,
    # sigma11, sigma12, sigma22, sigma21
    # derived para
    # lamda(t,i,t,j), prob(t,i,t,j)
    # series_para
    # length N,

    N = len(series)
    if N < 12:
        print "too short series"
        return [0 for i in range(6)] +[1.0 for i in range(6,10)]
    paras = [1.0 for i in range(10)]
    improve = [0.0 for i in range(10)]
    lamda = [[[] for i in range(N)] for j in range(6)] # 0 for 11, 1 for 12, 2 for 22, 3 for 21, 4 for 1, 5 for 2
    prob = [[[] for i in range(N)] for j in range(6)]

    paras = paras_init(paras)
    improve = [THRESHOLD[i]+0.1 for i in range(len(THRESHOLD)) ]
    while not converge(improve, THRESHOLD):
        paras_old = paras
        lamda = lamda_update(paras, lamda)
        prob  = prob_update(lamda, series, prob)
        paras_new = paras_update(prob, series, paras)
        improve = improve_cal(paras_old, paras_new)
        paras = paras_new

    return paras

def paras_init(paras):

    paras[0] = 0.8 # alpha1
    paras[1] = 0.8 # alpha2
    paras[2] = 0.8 # yita11
    paras[3] = 0.8 # yita12
    paras[4] = 0.8 # yita22
    paras[5] = 0.8 # yita21
    paras[6] = 1.0 # sigma11
    paras[7] = 1.0 # sigma12
    paras[8] = 1.0 # sigma22
    paras[9] = 1.0 # sigma21
    return paras


def lamda_update(paras, lamda):

    for n in range(4):
        for i in range(len(lamda[n])):
            j = i+1
            lamda[n][i] = paras[n+2]*j/paras[n+6]/paras[n+6]*math.exp(0-(j*j)/2/paras[n+6]/paras[n+6])
    for n in [4,5]:
        for i in range(len(lamda[n])):
            lamda[n][i] = paras[n-4]/pow(j,(paras[n-4]+1))
    return lamda


def prob_update(lamda, series, prob):
    """
    altitude at each time index is included in prob, so that later sum will be over time index
    """
    for i in range(len(series)):
        lamda_sum = [0 for n in range(2)]
        lamda_sum[0] = sum([series[j][0]*lamda[0][i-j] + series[j][1]*lamda[3][i-j] for j in range(i)]) + lamda[4][i]
        lamda_sum[1] = sum([series[j][1]*lamda[2][i-j] + series[j][0]*lamda[1][i-j] for j in range(i)]) + lamda[5][i]
        for n in range(4):
            if n ==0 or n == 3:
                tag = 0
            else:
                tag = 1
            prob[n][i] = [series[j][tag]*lamda[n][i-j]/lamda_sum[tag] for j in range(i)]
        for n in [4,5]:
            prob[n][i] = lamda[n][i]/lamda_sum[n-4]
    return prob


def paras_update(prob, series, paras):

    paras_new = [0 for i in range(len(paras))]
    for n in [0,1]:
        denominator = sum([math.log(i+1)*prob[n+4][i] for i in range(len(prob[n+4]))])
        if denominator == 0:
            paras_new[n] = 0
        else:
            paras_new[n] = sum(prob[n+4])/denominator
    for n in range(2,6):
        try:
            series_sum = [sum(series[:][i]) for i in [0,1]]
        except IndexError, e:
            print series
        if series_sum[(n-2)/2] == 0:
            paras_new[n] = 0
        else:
            paras_new[n] = sum([sum(prob[n-2][i]) for i in range(len(prob[n-2]))])/series_sum[(n-2)/2]
    for n in range(6,10):
        denominator = sum([sum(prob[n-6][i]) for i in range(len(prob[n-6]))])
        if denominator == 0:
            paras_new[n] = paras[n]
        else:
            paras_new[n] = sum([sum([prob[n-6][i][j]*math.pow(i-j,2) for j in range(i)])for i in range(len(prob[n-6]))])/denominator
    return paras_new


def improve_cal(paras_old, paras_new):

    improve = [0 for i in range(len(paras_old))]
    for i in range(len(paras_old)):
        if paras_old[i] > 0:
            improve[i] = (paras_new[i]-paras_old[i])/paras_old[i]
        else:
            improve[i] = paras_new[i] - paras_old[i]
    return improve


def converge(improve, THRESHOLD):
    for i in range(len(improve)):
        if(improve[i]>THRESHOLD[i]):
            return 0
    return 1

def series_prep():
    pass
    
if __name__ == "__main__":
    Hawkes_predict(series, tr, tp)