import numpy as np

f = open('Commands.txt','w')
dl_b = np.linspace(0.1, 1.0, 10)

for i in range(10):
	o = dl_b[i]
	cmdline = 'nohup '+'./simDER' + ' option.txt '  + '--' + ' hangLength ' + '%.4f'% o +'\n'
	f.write(cmdline)

f.close()
