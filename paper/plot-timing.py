import pandas as pa
import sys
from matplotlib import pyplot as plt
import numpy as np

plt.close()
f, ax = plt.subplots(1, figsize=(5,3))

df = pa.read_table(sys.argv[1])

pf = df.ix[~df.parallel]
pt = df.ix[df.parallel]


ax.set_ylabel('Intervals / Second')

x = 0.1 + np.arange(3)
width = 0.35

rects0 = ax.bar(x, 3327. / pt.time, width, color='r', edgecolor='r', alpha=0.91)
rects1 = ax.bar(x[:len(pf.time)] + width + 0.03, 3327. / pf.time, width,
        color='b', edgecolor='b', alpha=0.91)

ax.set_xticks(x + width)
ax.set_ylim(0, (3327. / min(pf.time)) + 50)
ax.set_xticklabels(["%s %s" % tup for tup in zip(pt['loc'], pt.instance)])

ax.legend((rects0[0], rects1[0]), ('parallel', 'not-parallel'))

plt.savefig('manuscript-latex/figure1.eps')
plt.savefig('figure1.pdf')
