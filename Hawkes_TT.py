import math
import random

THRESHOLD = [0.01 for i in range(10)]
MAX_ITERATION = 1000
T0 = 1.5 # time of first stamp (simple test shows different T0>1.0 influence predicted value only a little)
EPSILON = 1.0e-7
DIM = 2
key_map = [[0 for nn in range(2)] for nnn in range(2)]
for nr in range(2):
    for np in range(2):
        key_map[nr][np] = int(nr + np * (1- nr*2) + 2 * nr)
print key_map

def Hawkes_predict(series, tr, tp):

    if len(series)<= tp:
        print "too short series"
        return [-1,-1],[-1,-1]
    orgvalue = [sum([series[i][n] for i in range(tp)]) for n in range(2)]
    orgref = [sum([series[i][n] for i in range(tr)]) for n in range(2)]
    if orgref[0] == 0 or orgref[1] == 0:
        print "zero series"
        return [-1,-1],[-1,-1]
    paras = Hawkes_fit_EM(series[:tr])


    # ########### according to TingTing He et al. paper #######
    # predict = [0 for n in range(2)]

    # for n in range(2):
    #     predict[n] = 1 - pow((tp-tr),-paras[n]) + paras[2*(n+1)]*(1-math.exp(-pow(tp-tr,2)/2/pow(paras[n*2+6],2)))\
    #     + paras[5-n*2]*(1-math.exp(-pow(tp-tr,2)/2/pow(paras[9-n*2],2)))
    #     predict[n] += orgref[n]
    ######## According to definition of Hawkes process #######
    predict = predict_cal(series, paras, tr, tp)

    print "origin value for series 0: ", orgvalue[0], ", predicted value", predict[0], ", origin ref value:", orgref[0]
    print "origin value for series 1: ", orgvalue[1], ", predicted value", predict[1], ", origin ref value:", orgref[1]
    abs_err = [abs(predict[i]-orgvalue[i]) for i in range(2)]
    rel_err = [0 for i in range(2)]
    for i in range(2):
        if orgvalue[i] > 0:
            rel_err[i] = abs_err[i]/orgvalue[i]
        else:
            rel_err[i] = abs_err[i]
    print "abs_err", abs_err
    print "rel_err", rel_err
    return abs_err, rel_err

def predict_cal(series, paras, tr, tp):
    dimension = len(series[0])
    orgref = [sum([series[i][n] for i in range(tr)]) for n in range(2)]
    predict = [orgref[n] for n in range(dimension)]
    for np in range(dimension):
        for nr in range(dimension):
            influence_factor = [influence_factor_cal(paras, t, tr, tp, nr, np) for t in range(tr)]
            # print "nr, np:", nr, np
            # print influence_factor
            predict[np] += sum([influence_factor[t] * series[t][nr] for t in range(tr)])
        predict[np] += influence_factor_cal(paras, 0, tr, tp, -1, np)
    return predict

def influence_factor_cal(paras, t, tr, tp, nr, np):

    if nr == -1:
        return  (pow(tr+T0, -paras[np]) - pow(tp+T0, -paras[np]))
    else:
        return cross_influence_cal(paras[key_map[nr][np]+2], paras[key_map[nr][np]+6], t, tr, tp)

def cross_influence_cal(yita, sigma, t, tr, tp):
    return yita * (math.exp(-pow(tr-t,2)/2/pow(sigma,2)) - math.exp(-pow(tp-t ,2)/2/pow(sigma,2)))



def Hawkes_fit_EM(series):

    # para_list
    # alpha1, alpha2, yita11, yita12, yita22, yita21,
    # sigma11, sigma12, sigma22, sigma21
    # derived para
    # lamda(t,i,t,j), prob(t,i,t,j)
    # series_para
    # length N,

    N = len(series)
    paras = [1.0 for i in range(10)]
    improve = [0.0 for i in range(10)]
    lamda = [[0 for i in range(N)] for j in range(6)] # 0 for 11, 1 for 12, 2 for 22, 3 for 21, 4 for 1, 5 for 2
    prob = [[0 for i in range(N)] for j in range(6)]

    paras_best = [1.0 for i in range(10)]
    likelihood_largest = 0
    for randomtry in range(10):
        paras = paras_init(paras)
        paras_old = paras
        lamda = lamda_update(paras, lamda)
        improve = [THRESHOLD[i]+0.1 for i in range(len(THRESHOLD)) ]
        likelihood_old = 0
        Q_old = 0
        iteration = 0
        while not converge(improve, THRESHOLD):

            ############## likelihood test #############
            if iteration == 0:
                pass
            else:
                likelihood = likelihood_update(lamda, series, paras)
                # print "likelihood value:", likelihood
                if likelihood_old != 0 and likelihood < likelihood_old:
                    print "error: decreasing likelihood"
                    print paras_old
                    # return paras_old
                likelihood_old = likelihood

            ############# E step #####################

            prob = prob_update(lamda, series, prob)

            ############# M step #####################

            paras_old[:] = paras[:]
            ########## Q test in M step #########
            # Q_old = Q_update(lamda, series, prob, paras)
            paras_new = paras_update(prob, series, paras)
            lamda = lamda_update(paras_new, lamda)
            improve = improve_cal(paras_old, paras_new)
            paras = paras_new
            ########### Q test in M step #########
            # Q = Q_update(lamda, series, prob, paras)
            # if Q< Q_old:
            #     print "error: decreasing Q"
            #     print "Q_old:", Q_old, " Q:", Q

            iteration += 1
            if iteration >= MAX_ITERATION:
                break
        if converge(improve, THRESHOLD):
            print "convergent"
            print "iteration:", iteration
            print "paras:", paras
        else:
            print "not convergent"
            print "paras:", paras
            print "paras_old", paras_old
        print "likelihood:", likelihood
        if likelihood_largest == 0:
            likelihood_largest = likelihood
            paras_best = paras[:]
        else:
            if likelihood > likelihood_largest:
                paras_best = paras[:]
                likelihood_largest = likelihood
        print "largest likelihood:", likelihood_largest
        print "best paras:", paras_best
    return paras_best

def paras_init(paras):

    paras[0] = random.random() # alpha1
    paras[1] = random.random() # alpha2
    paras[2] = random.random() # yita11
    paras[3] = random.random() # yita12
    paras[4] = random.random() # yita22
    paras[5] = random.random() # yita21
    paras[6] = random.random() # sigma11
    paras[7] = random.random() # sigma12
    paras[8] = random.random() # sigma22
    paras[9] = random.random() # sigma21
    return paras


def lamda_update(paras, lamda):

    for n in range(4):
        for i in range(len(lamda[n])-1):
            j = float(i+1)
            lamda[n][i+1] = paras[n+2]*j/paras[n+6]/paras[n+6]*math.exp(0-(j*j)/2/paras[n+6]/paras[n+6])
        lamda[n][0] = 0
    for n in [4,5]:
        for i in range(len(lamda[n])):
            #### here series begin from 0
            j = float(i+T0)
            try:
                lamda[n][i] = paras[n-4]/pow(j,(paras[n-4]+1))
            except OverflowError, e:
                print "overflow of immigrant term,", paras[n-4]
    return lamda


def prob_update(lamda, series, prob):
    """
    reference altitude at each time index is included in prob, so that later sum will be over time index
    """
    for i in range(len(series)):
        lamda_sum = [0 for n in range(2)]
        lamda_sum[0] = sum([series[j][0]*lamda[0][i-j] + series[j][1]*lamda[3][i-j] for j in range(i)]) + lamda[4][i]
        lamda_sum[1] = sum([series[j][1]*lamda[2][i-j] + series[j][0]*lamda[1][i-j] for j in range(i)]) + lamda[5][i]
        for n in range(4):
            if n == 0 or n == 1:
                tag = 0
            else:
                tag = 1
            if n == 0 or n == 3:
                tagtarget = 0
            else:
                tagtarget = 1
            prob[n][i] = [series[j][tag]*lamda[n][i-j]/lamda_sum[tagtarget] for j in range(i)]
        for n in [4,5]:
            prob[n][i] = lamda[n][i]/lamda_sum[n-4]
    return prob


def paras_update(prob, series, paras):

    paras_new = [paras[i] for i in range(len(paras))]
    for n in [0,1]:
        denominator = sum([math.log(float(i+T0))*prob[n+4][i]*series[i][n] for i in range(len(series))])
        if denominator == 0:
            paras_new[n] = 0
        else:
            paras_new[n] = sum([prob[n+4][j]*series[j][n] for j in range(len(series))])/denominator

    series_sum = [sum([series[j][i] for j in range(len(series))]) for i in [0, 1]]
    for n in range(2, 6):
        denominator = series_sum[(n-2)/2]
        if denominator ==0:
            paras_new[n] = 0
        else:
            if n == 2 or n == 5:
                tag = 0
            else:
                tag = 1
            paras_new[n] = sum([sum(prob[n-2][i])*series[i][tag] for i in range(len(prob[n-2]))])/denominator
    for n in range(6,10):
        if n == 6 or n == 9:
            tag = 0
        else:
            tag = 1
        denominator = sum([sum(prob[n - 6][i])*series[i][tag] for i in range(len(prob[n - 6]))])
        if denominator == 0:
            paras_new[n] = paras[n]
        else:
            paras_new[n] = math.sqrt(sum([sum([prob[n - 6][i][j] * math.pow(i - j, 2) / 2.0 for j in range(i)])*series[i][tag]\
                                for i in range(len(prob[n - 6]))]) / denominator)
    return paras_new


def improve_cal(paras_old, paras_new):

    improve = [0 for i in range(len(paras_old))]
    for i in range(len(paras_old)):
        if paras_old[i] > EPSILON:
            improve[i] = abs(paras_new[i]-paras_old[i])/paras_old[i]
        else:
            improve[i] = abs(paras_new[i] - paras_old[i])
    return improve


def converge(improve, THRESHOLD):
    for i in range(len(improve)):
        if(improve[i]>THRESHOLD[i]):
            return 0
    return 1

def series_prep():
    series = [[1,2],[3,4],[5,5],[6,9],[12,18],[20,36],[55,98],[150,140],[34,55],[13,34],[11,11],[10,10],[10,10],[9,9],[8,8],[8,8],[7,7],[7,7],[7,7],[6,6],[6,6],[6,6],[6,6],[4,4],[4,4]]
    # for i in range(len(series)):
    #     temp = series[i][0]
    #     series[i][0] = series[i][1]
    #     series[i][1] = temp
    return series

def likelihood_update(lamda, series, paras):
    likelihood = 0
    lamda_sum = [[0 for n in range(2)] for i in range(len(series))]
    series_sum = [sum([series[j][i] for j in range(len(series))]) for i in [0, 1]]
    for i in range(len(series)):
        lamda_sum[i][0] = sum([series[j][0] * lamda[0][i - j] + series[j][1] * lamda[3][i - j] for j in range(i)]) + lamda[4][i]
        lamda_sum[i][1] = sum([series[j][1] * lamda[2][i - j] + series[j][0] * lamda[1][i - j] for j in range(i)]) + lamda[5][i]
    for n in range(2):
        try:
            likelihood += (sum([series[i][n]*math.log(lamda_sum[i][n]) for i in range(len(series))])-
                  series_sum[n]*(1+paras[n*2+2])- series_sum[1-n]*paras[5-n*2])
        except ValueError,e:
            print series_sum[n], paras[n*2+2], series_sum[1-n], paras[5-n*2]
    return likelihood

def Q_update(lamda, series, prob, paras):
    Q = 0
    N = len(series)
    term = [[0 for j in range(N)] for i in range(6)]
    for nn in range(4):
        for i in range(N):
            term[nn][i] = sum([prob[nn][i][j]*math.log(lamda[nn][i-j]) for j in range(i)])
    for nn in [4,5]:
        for i in range(N):
            term[nn][i] = prob[nn][i]*math.log(lamda[nn][i])
    for nn in range(6):
        if nn == 0 or nn == 3 or nn == 4:
            tag = 0
        else:
            tag = 1
        Q += sum([term[nn][i] * series[i][tag] for i in range(N)])
    series_sum = [sum([series[j][i] for j in range(len(series))]) for i in [0, 1]]
    for n in range(2):
        Q = Q - series_sum[n]*(1+paras[n*2+2])- series_sum[1-n]*paras[5-n*2]
    return Q

if __name__ == "__main__":
    Hawkes_predict(series_prep(), 20,21)
