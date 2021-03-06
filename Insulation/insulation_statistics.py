from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
start_time = time.process_time()

binsize = 25000

# Необходимо: 
# 1) размеры хромосом .chrom.sizes (указать путь в переменной chrfile); 
# 2) набор доменов (указать путь в переменной domains_from);
# 3) профиль инсуляции .BedGraph (указать путь в переменной feature_from).

# Для каждой хромосомы оценивается распределение коэф. инсуляции относительно границ доменов.
# Интервалы = домены + междоменные промежутки, нужны только для определения границ. 
# Правее от границы расстояние положительно, левее - отрицательно.
# При равенстве расстояний до границ выбирается положительное значение.

def domains_shuffler(intervals):

    start_position = int(intervals.iat[0,0])
    
    distances = pd.DataFrame([])    
    for index, row in intervals.iterrows():
        i = int(row.iat[0])
        j = int(row.iat[1])
        diff = abs(i-j)
        data = pd.Series([diff],dtype='int')
        distances = distances.append(data,ignore_index=True)
    distances = np.array(distances)
    
    np.random.shuffle(distances)
    
    shuffled_distances = [start_position]
    for u in distances:
        start_position = start_position+u
        shuffled_distances.insert(len(shuffled_distances),start_position[0])

    shuffled = pd.DataFrame(shuffled_distances)

    shuffled = pd.concat([shuffled, shuffled], axis=1, ignore_index=True)
    shuffled[1].index = shuffled[1].index-1
    shuffled = pd.concat([shuffled[0], shuffled[1]], axis=1,ignore_index=True)
    shuffled.drop(shuffled.head(1).index,inplace=True)
    shuffled.drop(shuffled.tail(1).index,inplace=True)
    shuffled = np.array(shuffled)
    return shuffled

def distant_lights(intervals,feature_track,N):
    
    intervals = np.array(intervals)
    
    feature_values = pd.DataFrame([])

    for index, row in feature_track.iterrows():
        i = row.iat[3]
        data = pd.Series([i])
        feature_values = feature_values.append(data,ignore_index=True)

    feature_values = np.array(feature_values[0])
    feature_values = feature_values.astype(np.float)

    distdict = defaultdict(int)
    result_dict = defaultdict(int)
    
    for istart in np.arange(0,(intervals[0][0]),1):
        distance = istart-(intervals[0][0])
        distdict[istart] += distance

    for interval in intervals:
        a = interval[0]
        b = interval[1]
        c = (a+b)/2
        rangeofi = np.arange(a, b+1, 1)
        
    # Правее от границы расстояние положительно, левее - отрицательно.
    # При равенстве расстояний до границ выбирается положительное значение.

        for i in rangeofi:
            if i <= c:
                distance = i-a
                distdict[i] = distance
            else:
                distance = i-b
                distdict[i] = distance

    for iend in np.arange((intervals[intervals.shape[0]-1][1]),N):
        distance = iend - (intervals[intervals.shape[0]-1][1])
        distdict[iend] += distance   
        
    c = Counter(distdict.values())

    for key in distdict:
        key = int(key)
        result_dict[distdict[key]] += feature_values[key]/c[distdict[key]]

    df = pd.DataFrame([result_dict])
    df = df.transpose()
    df[1] = df[0]
    df[0] = df.index
    df = df.reset_index(drop=True)
    # D - количество бинов справа и слева от границы.
    # Выбрано [-4,4], т.к. средний размер выделенных доменов около 5 бинов. 
    Dx = [-4,4]
    df = df.loc[(df[0] >= Dx[0]) & (df[0] <= Dx[1])]
    df = df.reset_index(drop=True)
    return df

######################################################################

chrfile = pd.read_csv('/home/al/Anop/ChrSizesRenew/Acol_VARYA.chrom.sizes',sep='\t',header = None)
for chrom in chrfile[0]:
    chrom = str(chrom)
    domains_from = '/home/al/Anop/Domains/Acol_' + chrom + '_domains.tsv'
    feature_from = '/home/al/Anop/Insulation/Acol_' + chrom + '_ins.bedGraph'
    input_domains = pd.read_csv(domains_from,sep='\t',header = None)
    feature_track = pd.read_csv(feature_from,sep='\t',header = None)

    N = chrfile.loc[(chrfile[0] == chrom)]
    N = N.reset_index(drop=True)
    N = int(round(N[1][0]/binsize)-1)
    
    binned_domains = pd.DataFrame([])
    for index, row in input_domains.iterrows():
        i = int(row.iat[1])
        j = int(row.iat[2])
        j += 1
        i /= 25000
        j /= 25000
        data = pd.Series([i,j],dtype='int')
        binned_domains = binned_domains.append(data,ignore_index=True)
    
    binned_domains = np.array(binned_domains)
    rolled_domains = np.roll(binned_domains, 1)
    for x in rolled_domains:
        if x[0] >= x[1]:
            rolled_domains = np.delete(rolled_domains, np.where(x == rolled_domains)[0][0], 0)

    binned_domains = pd.DataFrame(binned_domains)
    rolled_domains = pd.DataFrame(rolled_domains)
    intervals = pd.concat([binned_domains, rolled_domains], axis=0, ignore_index=True)
    intervals = intervals.sort_values(0)
    intervals = intervals.reset_index(drop=True)

    control = pd.DataFrame([])
    
    iterations = 10

    while iterations > 0:
        print(iterations)
        iterations -= 1
        shuffled_domains = domains_shuffler(intervals)
        control1 = distant_lights(shuffled_domains, feature_track,N)
        control = pd.concat([control, control1[1]], axis=1, ignore_index=True)
        control.columns = range(control.shape[1])

    mean_control = control.mean(axis=1)
    std_control = control.std(axis=1)
    plus3std = mean_control+3*std_control
    minus3std = mean_control-3*std_control

    true_result = distant_lights(intervals,feature_track,N)

    result = pd.concat([true_result, mean_control,plus3std,minus3std], axis=1, ignore_index=True)
    result.columns = (chrom, 'avgtrue','avgctrl', '+3s', '-3s')

    print(result)

    fig, ax = plt.subplots(1, 1, sharex=True)
    
    fig.patch.set_facecolor('#FFFFFF')
    ax.plot(result[chrom],result['avgctrl'],color = "k", linewidth = 1)
    ax.plot(result[chrom],result['+3s'],color = "#C0C0C0", linewidth = 0.8)
    ax.plot(result[chrom],result['-3s'],color = "#C0C0C0", linewidth = 0.8)
    ax.fill_between(result[chrom],result['+3s'],result['-3s'],color = "#D3D3D3")
    ax.plot(result[chrom],result['avgtrue'],color = "k", linestyle="-",linewidth = 2)
    ax.tick_params(direction='inout')
    plt.show()
    figure_to = '/home/al/Anop/Graphs/Acol_' + chrom + '_Insul.png'
    df_to = '/home/al/Anop/InsStat/Acol_' + chrom + '_Insul.tsv'
    fig.savefig(figure_to,dpi=180)
    result.to_csv(df_to, sep='\t', float_format='%s',index=False)

print ("{:g} s".format(time.process_time() - start_time))
