#!/usr/bin/env python

import os
import shutil

header = "OUTPUT AlTi.train\n\nTYPES\n2\nAl 3461.4286  ! eV\nTi 3396.2232  ! eV\n\nSETUPS\nAl   Al.fingerprint.stp\nTi Ti.fingerprint.stp\n\nFILES\n"

files = os.listdir('gen/xsf')

f = open('generate.in', 'w')

try:
	os.mkdir('train')
	os.mkdir('test')
except:
	print "Directories: train, test already made...skipping"

f.write(header)
f.write(str(len(files)) + "\n")
train_count = 0
test_count = 0
for file in files:
	file_count = int(file[:-4].split('_')[-1])
	if file_count >= 5000 and file_count < 5700:
		print file_count
		new_fname = 'train/trainingConfig_' + str(train_count) + '.xsf'
		train_count += 1
		shutil.copyfile('gen/xsf/' + file, new_fname)
		f.write(new_fname + "\n")
	elif file_count >= 5700:
		new_fname = 'test/testingConfig_' + str(test_count) + '.xsf'
		test_count += 1
		shutil.copyfile('gen/xsf/' + file, new_fname)
f.close()