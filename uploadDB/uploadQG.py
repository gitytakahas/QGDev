#! /usr/bin/env python

import os
import shutil
import subprocess
import urllib
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--file', type='string', action='store', dest='file', help='the database file to be uploaded')
parser.add_option('-o', '--offline', action='store_true', dest='offline', default=False, help='Upload to the offline DB instead of PREP')
(options, args) = parser.parse_args()

# Check if input file exists"
if not os.path.isfile(options.file): 
  print file + ' does not exist!'
  exit(1) 

# Choose the database
if options.offline: database = 'oracle://cms_orcon_prod/CMS_COND_31X_PHYSICSTOOLS'
else:               database = 'oracle://cms_orcoff_prep/CMS_COND_TEMP'
print "Uploading to database: " + database
print

# Getting tags in the .db file
tags = subprocess.check_output(['cmscond_list_iov', '-c', 'sqlite_file:' + options.file, '-a']).split('\n')
tags.pop()
print 'Tags in ' + options.file + ':'

# Preparing JSON files
for tag in tags:
  print tag
  templateFile = open('QGL_template.txt', 'r')
  outputFile   = open(tag + '.txt', 'w')
  contents     = templateFile.read().replace('TAGNAME', tag).replace('DATABASE',database)
  outputFile.write(contents)
print

# Uploading
uploadScript = os.path.abspath('upload.py')
for tag in tags :
  shutil.copy(options.file, tag + '.db')
  try:
    if options.offline: print subprocess.check_output( [uploadScript, tag + '.db', '--backend', 'offline'], stderr=subprocess.STDOUT)
    else:               print subprocess.check_output( [uploadScript, tag + '.db'], stderr=subprocess.STDOUT)
  except subprocess.CalledProcessError as e:
    print e.output
    exit(1)

# Clean up
for tag in tags :
  os.remove(tag + '.db')
  os.remove(tag + '.txt')
