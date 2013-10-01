#!/usr/bin/python
import fileinput
import re
for line in fileinput.input():
    line = re.sub('\n', '', line)
    if (re.search('<->', line)):
        line = re.sub(r'(\S\')', r'\1 ', line)
        line = re.sub(r'<->', r' <-> ', line)
        line = re.sub(r', +([A-Z])', r',\1', line)
        els = re.split(' +', line)
        print(els[0] + " " + els[1] + " -> " + els[3] + " @ " + els[5])
        print(els[0] + " " + els[3] + " -> " + els[1] + " @ " + els[7])
        ## print(els)
    else:
        if (re.search('\A\s*%init', line)):
            line = re.sub('\*', '', line)
            els = re.split(' +', line)
            m = re.match('%init:.*\s+\((.*)\)', line)
            agent = m.group(1)
            agent = re.sub(r'~u,', r'~u~p,', agent)
            agent = re.sub(r'~u\)', r'~u~p)', agent)
            print('%agent: ' + agent)
            line = re.sub(r'~u~p', r'~u', line)
            print(line)
        else:
            print(line)
