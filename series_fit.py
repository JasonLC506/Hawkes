import cPickle
import Hawkes_TT

def series_fetch():
    file1 = open("viewCount","r")
    list_series_1 = cPickle.load(file1)
    file1.close()
    file2 = open("likeCount","r")
    list_series_2 = cPickle.load(file2)
    file2.close()
    list_series = []
    for i in range(len(list_series_2)):
        series_1 = list_series_1[i]
        series_2 = list_series_2[i]
        L = len(series_1)
        if L > len(series_2):
            L = len(series_2)
        series = [[series_1[j], series_2[j]] for j in range(L)]
        list_series.append(series)
    return list_series

def batch_predict(list_series):
    abs_errs =[]
    rel_errs = []
    N = len(list_series)
    for series in list_series:
        abs_err, rel_err = Hawkes_TT.Hawkes_predict(series, 12,24)
        if abs_err[0]<0:
            N -= 1
            continue
        abs_errs.append(abs_err)
        rel_errs.append(rel_err)

    mean_abs_err = [sum(abs_errs[:][i]) for i in range(2)]
    mean_rel_err = [sum(rel_errs[:][i]) for i in range(2)]
    print mean_abs_err, mean_rel_err
    resultfile = open("Hawkes_TT_result","a")
    resultfile.write("tr: %d, tp: %d\n" %(12,24))
    resultfile.write(str(mean_abs_err)+"\n")
    resultfile.write(str(mean_rel_err)+"\n")
    resultfile.close()
    
if __name__ == "__main__":
    batch_predict(series_fetch())