fp_write = open("plot_mod.txt","w")
with open('plot.txt') as fp:
	for ln in fp:
		word = ln.split(' ')
		fp_write.write(word[1]+' '+word[2]+' '+word[3])
fp_write.close()