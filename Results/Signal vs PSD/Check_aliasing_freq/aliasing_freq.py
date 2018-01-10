

Fs = (25e6/128)

#for freq in [2e3, 1e4, 2e4, 3e4, 5e4]:
#    print '\nFREQ:', freq/1000, 'kHz'
#    for i in range(10):
#        for n in range(10):
#            if 0 < (int((i*Fs - n * freq)/1000.)) < 100:
#                print 'for', i, 'times Fs, n = ', n, ': ', int((i*Fs - n * freq)/1000.)
#            if 0 < (int((i*freq - n * Fs)/1000.)) < 100:
#                print 'for', i, 'times Fs, n = ', n, ': ', int((i*freq - n * Fs)/1000.)


#(mail Robert 05-05-2017): "There always seem to be some resonances at frequencies m/n*f_sample"

for freq in [2e3, 1e4, 2e4, 3e4, 5e4]:
    print '\nFREQ:', freq/1000, 'kHz'
    for m in range(10):
        for n in range(10):
            n +=1
            if 0 < (int((float(m)/float(n)*Fs)/1000.)) < 100:
                print 'for m/n = ', m,'/',n, ': ', int((float(m)/float(n)*Fs)/1000.)