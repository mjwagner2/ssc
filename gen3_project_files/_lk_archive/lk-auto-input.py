import re

ssc_src = open("../ssc/cmod_tcsmolten_salt.cpp", 'r').read()

all_ssc = re.findall(" \{ SSC_INPUT, .*?\"(.*?)\"", ssc_src) + re.findall(" \{ SSC_INOUT, .*?\"(.*?)\"", ssc_src)

lk_src = open("mw.lk",'r').read()

all_lk = re.findall("\n(var\( '(.*?)',.*?;)", lk_src, re.MULTILINE)
all_lk_vars = [v[1] for v in all_lk]

fout = open('mw_mod.lk', 'w')

for line,var in all_lk:
    if "adjust:" in var:
        continue
    if not var in all_ssc:
        line = line.replace("\n","\n//")
        lk_src = lk_src.replace(line, "//"+line)

fout.write(lk_src)
fout.close()