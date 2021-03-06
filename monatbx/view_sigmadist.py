from iotbx import reflection_file_reader
import sys
import matplotlib.pyplot as plt
import numpy as np
from cctbx.array_family import flex

hklin = sys.argv[1]
reflection_file = reflection_file_reader.any_reflection_file(hklin)
miller_arrays = reflection_file.as_miller_arrays()
num_bins = 20
for miller_array in miller_arrays:
  if miller_array.sigmas() is not None:
    miller_array.show_summary()
    mean_sigI = np.mean(miller_array.sigmas())
    mean_I = np.mean(miller_array.data())
    std_I = np.std(miller_array.data())
    ma_1 = miller_array.select((miller_array.sigmas() > 0) & (miller_array.sigmas() <= mean_sigI))
    ma_2 = miller_array.select((miller_array.sigmas() > mean_sigI))
    ma_3 = miller_array.select((miller_array.data() < (mean_I - (0*std_I))))
    #plot
    plt.subplot(321)
    x = miller_array.sigmas().as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
    plt.ylabel('Frequencies')
    plt.title('sigI distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
    plt.subplot(322)
    plt.scatter(miller_array.data(), miller_array.sigmas(), s=10, c='b', marker='o')
    plt.subplot(323)
    s = 1/(ma_1.d_spacings().data()**2)
    x = s.as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
    plt.title('Which resolution have low-medium error (sigI <= %5.2f)'%(mean_sigI))
    plt.subplot(324)
    s = 1/(ma_2.d_spacings().data()**2)
    x = s.as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
    plt.title('Which resolution have high error (sigI > %5.2f)'%(mean_sigI))
    plt.subplot(325)
    s = 1/(ma_3.d_spacings().data()**2)
    x = s.as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
    plt.title('Where are the weak reflections')
    plt.subplot(326)
    x = ma_3.sigmas().as_numpy_array()
    mu = np.mean(x)
    med = np.median(x)
    sigma = np.std(x)
    plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
    plt.ylabel('Frequencies')
    plt.title('Do weak reflections have high error?\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
    #plt.show()

    #bin reflection set and view histogram of sigI and I vs. sigI scater plot
    binner = miller_array.setup_binner(n_bins=10)
    binner_indices = binner.bin_indices()
    one_over_dsqr_set = flex.double()
    mean_i_over_sig_i_set = flex.double()
    med_i_over_sig_i_set = flex.double()
    mean_ios_strong_set = flex.double()
    med_ios_strong_set = flex.double()
    mean_ios_weak_set = flex.double()
    med_ios_weak_set = flex.double()
    for i in range(1, binner.n_bins_used()+1):
      i_binner = (binner_indices == i)
      ma_bin = miller_array.select(i_binner)
      plt.subplot(5,4,i*2-1)
      x = ma_bin.sigmas().as_numpy_array()
      mu, med, sigma = (np.mean(x), np.median(x), np.std(x))
      plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
      plt.ylabel('Frequencies')
      plt.title('Bin %2.0f mu %4.1f med %4.1f sig %4.1f' %(i, mu, med, sigma))
      plt.xlim([0,np.max(miller_array.sigmas())])
      plt.subplot(5,4,i*2)
      plt.scatter(ma_bin.data(), ma_bin.sigmas(), s=10, c='b', marker='o')
      cc_matrix = np.corrcoef(ma_bin.data(), ma_bin.sigmas())
      plt.title('Bin %2.0f: (%4.1f-%4.1f) CC=%5.2f'%(i, np.max(ma_bin.d_spacings().data()), np.min(ma_bin.d_spacings().data()),cc_matrix[0,1]))
      plt.xlim([0,np.max(miller_array.data())])
      plt.ylim([0,np.max(miller_array.sigmas())])
      one_over_dsqr_set.append(1/(ma_bin.d_min()**2))
      mean_i_over_sig_i_set.append(np.mean(ma_bin.data()/ma_bin.sigmas()))
      med_i_over_sig_i_set.append(np.median(ma_bin.data()/ma_bin.sigmas()))

      #view strong reflections
      ma_bin_strong = ma_bin.select(ma_bin.data()>(np.mean(ma_bin.data())+(np.std(ma_bin.data())*0.5)))
      mean_ios_strong_set.append(np.mean(ma_bin_strong.data()/ma_bin_strong.sigmas()))
      med_ios_strong_set.append(np.median(ma_bin_strong.data()/ma_bin_strong.sigmas()))
      #view weak reflections
      ma_bin_weak = ma_bin.select(ma_bin.data()<(np.mean(ma_bin.data())-(np.std(ma_bin.data())*0.5)))
      mean_ios_weak_set.append(np.mean(ma_bin_weak.data()/ma_bin_weak.sigmas()))
      med_ios_weak_set.append(np.median(ma_bin_weak.data()/ma_bin_weak.sigmas()))
      """
      #view strong reflections
      ma_bin_strong = ma_bin.select(ma_bin.sigmas()>(np.mean(ma_bin.sigmas())+(np.std(ma_bin.sigmas())*0.5)))
      mean_ios_strong_set.append(np.mean(ma_bin_strong.data()/ma_bin_strong.sigmas()))
      med_ios_strong_set.append(np.median(ma_bin_strong.data()/ma_bin_strong.sigmas()))
      #view weak reflections
      ma_bin_weak = ma_bin.select(ma_bin.sigmas()<(np.mean(ma_bin.sigmas())-(np.std(ma_bin.sigmas())*0.5)))
      mean_ios_weak_set.append(np.mean(ma_bin_weak.data()/ma_bin_weak.sigmas()))
      med_ios_weak_set.append(np.median(ma_bin_weak.data()/ma_bin_weak.sigmas()))
      """
    plt.show()
    #i/sigi by bins
    plt.plot(one_over_dsqr_set, mean_i_over_sig_i_set, color='navy')
    plt.plot(one_over_dsqr_set, med_i_over_sig_i_set, color='cornflowerblue')
    plt.plot(one_over_dsqr_set, mean_ios_strong_set, color='firebrick')
    plt.plot(one_over_dsqr_set, med_ios_strong_set, color='lightcoral')
    plt.plot(one_over_dsqr_set, mean_ios_weak_set, color='darksage')
    plt.plot(one_over_dsqr_set, med_ios_weak_set, color='yellowgreen')
    plt.grid(True)
    plt.show()
