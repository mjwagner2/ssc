import re

f1 = open("mw.lk", 'r').read()
f2 = open("mb.lk", 'r').read()

objs = re.findall("var\( '(.*?)', (.*?)\)", f1, re.M)

comps = {}

for name,val in objs:

    m = re.search("var\( '{:s}', (.*?)\)".format(name), f2, re.M)

    if m:
        if m.group(1) != val:
            comps[name] = [val, m.group(1)]

fout = open('compare.csv','w')
keys = sorted( comps.keys() )


for key in keys:
    fout.write("{:s},{:s},{:s}\n".format(key, comps[key][0], comps[key][1]))