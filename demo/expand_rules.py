#!/usr/bin/python
import fileinput
import re
out = ''
agent_str = ''
init_str = ''
obs_str = ''
misc_str = ''
for line in fileinput.input():
    line = re.sub('\n', '', line)
    if (re.search('<->', line)):
        line = re.sub(r'(\S\')', r'\1 ', line)
        line = re.sub(r'<->', r' <-> ', line)
        line = re.sub(r', +([A-Z])', r',\1', line)
        els = re.split(' +', line)
        forward_name = els[0]
        backward_name = re.sub(r'\'\Z', r"_diss'", forward_name)
        out += forward_name  + " " + els[1] + " -> " + els[3] + " @ " + els[5] + '\n'
        out += backward_name + " " + els[3] + " -> " + els[1] + " @ " + els[7] + '\n'
        ## print(els)
    else:
        if (re.search('\A\s*%init', line)):
            line = re.sub('\*', '', line)
            els = re.split(' +', line)
            m = re.match('%init:.*\s+\((.*)\)', line)
            agent = m.group(1)
            agent = re.sub(r'~u,', r'~u~p,', agent)
            agent = re.sub(r'~u\)', r'~u~p)', agent)
            agent_str += '%agent: ' + agent + '\n'
            init_str += re.sub(r'~u~p', r'~u', line) + '\n'
        else:
            if (re.search('%obs', line)):
                obs_str += line + '\n'
            else: 
                if (re.search('->', line)):
                    out += line + '\n'
                else:
                    misc_str += line + '\n'

print(agent_str)
print(out)
print(init_str)
print(obs_str)
print(misc_str)
