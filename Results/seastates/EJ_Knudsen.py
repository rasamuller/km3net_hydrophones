import sys
import math
import string
from matplotlib import pyplot as plt

ss = [44.5, 50, 55., 61.5, 64.5, 66.5, 70]

def Knudsen_noise(seastate = 0, frequency=1000):
    return ss[seastate] -17.* math.log10(frequency/1000.)

def main(argv):
    
    for seastate in range(len(ss)):
        print ''
        print ss[seastate]
        seastate = seastate + 1
        print 56 + 19*math.log10(seastate)
    
    sys.exit()    

    #seastate = string.atoi(argv[1])
    for seastate in range(len(ss)):
        F = []
        S = []
        for freq in xrange(1000, 25000,10):
            #print freq, Knudsen_noise(seastate, freq)
            F.append(freq)
            S.append(Knudsen_noise(seastate, freq))
        plt.plot(F,S)
    plt.show()

if __name__ == "__main__":
    main(sys.argv)
