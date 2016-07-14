import numpy
import seaborn
import matplotlib.pyplot as plt

lines=[i.strip().split(',') for i in open('es.csv').readlines()]

task=[]
region=[]
roisize=[]
cohend=[]
pctsig=[]
desc=[]

for i in range(len(lines)):
    task.append(lines[i][0])
    region.append(lines[i][1])
    desc.append(lines[i][0]+': '+lines[i][1])
    roisize.append(int(lines[i][2]))
    cohend.append([float(j) for j in lines[i][3:6]])
    pctsig.append([float(j) for j in lines[i][6:9]])

fig = plt.figure(figsize=(12,6))
ax=plt.gca()

cohend=numpy.array(cohend)
cohend_range=cohend[:,[0,2]]
cohend_range[:,0]=cohend[:,1]-cohend_range[:,0]
cohend_range[:,1]=cohend_range[:,1]-cohend[:,1]

pos=numpy.arange(cohend.shape[0])
plt.barh(pos,cohend[:,1],xerr=cohend_range.T,align='center',ecolor='red')
plt.yticks(pos,desc)
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontname('Arial')
    label.set_fontsize(16)
plt.xlabel("effect size (Cohen's d)",fontsize=18)
fig.tight_layout()
plt.savefig('cohensd_barplot.pdf')
