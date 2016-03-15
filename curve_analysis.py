#Created 2016-03-11 by QF
#This file analyzes a frequency curve in a way as objective and reproducible as possible.
#Takes for input a file query.csv as produced by the corpiscamis.py script. 

def curve_analysis(filename):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    occ_time = pd.read_csv(filename)

    a_k = 4 #number of additional decades on each side of the time range selected on which the algorithm tries to optimize the logit.

    min_nb = 5 #minimal number of points for the logit

    tol = 2 #number of successive decades allow to decrease in growth phase

    min_growth = 0.1 #threshold of the relative derivative above which the growth is assumed to be significative

    good_fit = 0.95 #minimal r^2 value above which the fit is considered to be good


    frequency_data = occ_time['Frequence moyenne']

    relative_der = [ (frequency_data[i + 1] - frequency_data[i - 1]) / (2.0 * frequency_data[i]) for i in range(1,len(frequency_data)-1) ]
    relative_der = [(frequency_data[1] - frequency_data[0]) / frequency_data[0]] + relative_der + [(frequency_data[len(frequency_data) - 1] - frequency_data[len(frequency_data) - 2]) / frequency_data[len(frequency_data) - 1]]

    growth = []
    for k in range(len(relative_der)):
        if relative_der[k] > min_growth:
            growth.append(k)

    growth_parts = []
    growth_part = []
    while growth != []:
        if growth_part != []:
            if growth[0] - growth_part[-1] <= tol:
                growth_part.append(growth[0])
                growth.pop(0)
            else:
                growth_parts.append(growth_part)
                growth_part= []
        else:
            growth_part.append(growth[0])
            growth.pop(0)
    growth_parts.append(growth_part)

    growth_parts = [p for p in growth_parts if len(p) >= min_nb]


    for part in range(len(growth_parts)):
        pairs = []
        r_values = []
        slopes = []
        born_1 = max(growth_parts[part][0] - a_k, 0)
        born_2 = min(growth_parts[part][-1] + a_k,len(frequency_data) - 1)
        for i in range(born_1,born_2-min_nb):
            for j in range(i+min_nb,born_2):
                i_min = i  
                i_max = j
                if all(frequency_data[k] < frequency_data[i_max] for k in range(i_min, i_max)) & all(frequency_data[k] > frequency_data[i_min] for k in range(i_min + 1, i_max)):
                    pairs.append([i_min,i_max])
                    logit = [np.log((frequency_data[i] - frequency_data[i_min]) / (frequency_data[i_max] - frequency_data[i])) for i in range(i_min + 1,i_max)]
                    x = np.arange(len(logit))
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x,logit)
                    r_values.append(r_value ** 2)
                    slopes.append(slope)

        good_pairs = [pairs[j] for j in range(len(pairs)) if r_values[j] > good_fit]
        if good_pairs == []:
            print('Period of investigation : %s - %s' % (occ_time['Decade'][born_1][:4], occ_time['Decade'][born_2][-4:]))
            print('The linear fit of the logit is rather poor.')
        else:
            good_rs = [r_values[j] for j in range(len(r_values)) if r_values[j] > good_fit]
            width = [good_pairs[j][1] - good_pairs[j][0] for j in range(len(good_pairs))]
            r_left = [good_rs[j] for j in range(len(width)) if width[j] == max(width)]
            indices = [j for j in range(len(width)) if width[j] == max(width)]
            selected_pair = good_pairs[indices[r_left.index(max(r_left))]]
            i_min, i_max = selected_pair[0], selected_pair[1]
            slope = slopes[pairs.index(selected_pair)]
            logit = [np.log((frequency_data[i] - frequency_data[i_min]) / (frequency_data[i_max] - frequency_data[i])) for i in range(i_min+1,i_max)]
            print('Period of rise : %s - %s' % (occ_time['Decade'][i_min][:4], occ_time['Decade'][i_max][-4:]))
            print('Rising time = %s' %(len(logit) + 2))
            print('r = %s' %max(r_left))
            print('slope = %.2f' %slope)
            occ_time[born_1:born_2].plot(kind='bar', x=['Decade'], y=['Frequence moyenne'])
            plt.title('Overall period of growth')
            occ_time[i_min:i_max].plot(kind='bar', x=['Decade'], y=['Frequence moyenne'], color='r')
            plt.title('Period of sigmoidal growth')
            plt.show()
            plt.plot(logit,'o')
            plt.title('Logit of the sigmoidal growth')
            plt.show()



