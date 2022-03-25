import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import sys
from sys import argv
from collections import OrderedDict

if len(argv) != 5:
	print("python structure_plot.py prefix sample.list 2 8")
	sys.exit()

color = ['#3CB371', '#CD5C5C', '#FF7F50', '#EE82EE', '#FFE4C4', '#4876FF', '#00B2EE', '#9F79EE', '#EE7AE9', '#EE9A00']

script, prefix, sample, start_num, end_num = argv
start_num = int(start_num)
end_num = int(end_num)
sample_list = OrderedDict()
x_major_tick = []
x_tick_pos = []
x_tick_label = []
with open(sample) as f:
	index = 0
	for line in f:
		index += 1
		info = line.strip('\n').split('\t')
		if info[1] not in sample_list:
			if not sample_list:
				sample_list[info[1]] = [index-0.5]
				group = info[1]
			else:
				sample_list[group].append(index - 0.5)
				group = info[1]
				sample_list[group] = [index-0.5]
	sample_list[group].append(index+0.5)
for key, value in sample_list.items():
	x_tick_label.append(key)
	x_tick_pos.append((value[0]+value[1])/2)
	for i in value:
		if i in x_major_tick: continue
		else: x_major_tick.append(i)
figure = plt.figure(figsize=(20,6))
n = 0
for num in range(start_num, end_num + 1):
	n += 1
	data = np.loadtxt('./%s.%s.Q' % (prefix, num), unpack=True, delimiter=' ', dtype=np.dtype('f8'))
	ax = figure.add_subplot(end_num - start_num + 1, 1 , n)
	for k_num in range(num):
		if k_num == 0:
			x_data = np.arange(1, len(data[k_num]) + 1)
			ax.bar(x_data, data[k_num], color=color[k_num], width=1)
			bottom_data = data[0]
		else:
			x_data = np.arange(1, len(data[k_num]) + 1)
			ax.bar(x_data, data[k_num], bottom=bottom_data, color=color[k_num], width=1)
			bottom_data += data[k_num]
	ax.set_yticks([0, 0.5, 1])
	ax.set_yticklabels([])
	ax.set_xticks([])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_linestyle(':')
	ax.spines['bottom'].set_position(('data', -0.1))
	ax.margins(x=0, y=0.2)
	ax.set_ylim(top=1)
	ax.set_ylabel('k=%s' % num, rotation='horizontal', labelpad=20, size=20, position=(-0.5, 0.35))
	ax.axhline(y=-0.1, linestyle='--', linewidth=4)
ax.spines['bottom'].set_visible(False)
ax.spines['bottom'].set_position(('data', 0.02))
ax.set_xticks(x_major_tick)
ax.set_xticklabels([])
ax.tick_params(axis='x', which='major', length=8, width=2.5)
ax.set_xticks(x_tick_pos, minor=True)
ax.set_xticklabels(x_tick_label, minor=True)
ax.tick_params(axis='x', which='minor', length=0, width=0, labelsize=20, pad=15)
#print(ax.get_xticks())

figure.subplots_adjust(bottom=0.1, top=0.95, left=0.05, right=0.98, hspace=0.2)
plt.savefig('structure.png', format='png')

