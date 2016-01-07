#!/usr/bin/env python
import re,sys,os

fn = open(sys.argv[1]) #blast result fmt 11
o = open(sys.argv[2], "w") #save result

#Query patter can not all " for more than one line
Query = re.compile("local str \"(Query_.*)\"")
s_name = re.compile("local str \"([0-9]+\|.*)")
sta = re.compile("([0-9-]{1,10})")

count_sname = 0
starts = {}
lens = {}
start_flag = False
len_flag = False
for f in fn:
	f = f.strip()
	#starts[query][count_sname][sname] = [starts]
	#lens[query][count_sname][sname] = [lens]
	if Query.match(f):
		query = Query.match(f).group(1)
		if query not in starts:
			starts[query] = {}
		if query not in lens:
			lens[query] = {}

        #extract sname: 1607466|Q7MWE4|Q7MWE4
        if s_name.match(f):
		#count here is the same as numbers of matching pairs (query vs subject)
                count_sname = count_sname + 1
		#serverl results of sname
		#"1234|LADJFL|ALKSDF" in one line
		#"123124|LADOUFOAISDUF|ALDFUOQWEI| haersfasd uiadsfu" in one line
		#"123124|ALDKFUAOISFOAISDUF00012309481|ER0000
		#01284" in two lines
                subject_name = s_name.match(f).group(1).split("\"")[0].split()[0]
                #print "%s %s %s" % (query, count_sname, subject_name)
                #print subject_name     
                #also correct with query
                #Here is the problem: For example for nid = 1009329
                #Query_6 VS 1009329 (the same query seq mathcing the same subject seq with differnt part of seq)
                #But you can't distinguish only from qname and sname.
                #So I add another count_sname in order to distinguish this difference
                starts[query][count_sname] = {}
                lens[query][count_sname] = {}
                starts[query][count_sname][subject_name] = ''
                lens[query][count_sname][subject_name] = ''
                #1471228|Q7QJZ2|Q7QJZ2
                #833 too...
                #print count_sname

	#The start points for detail mapping starting points of one sequence
        if f.startswith("starts {"):
                start_flag = True
                len_flag = False
                #correct! 833 paired sequence all together
                #print count
                continue

        #Extract starts sites (distinguish from lens using two flag)
        if sta.match(f) and start_flag==True and len_flag==False:
                #print "starts success %s" % sta.match(f).group(1)
                #if count not in starts:
                #       starts[count] = ''
                #print "In starts %s %s %s" % (query, count_sname, subject_name)
                #print seq
		try:
                	if len(starts[query][count_sname][subject_name]) == 0:
                        	starts[query][count_sname][subject_name] = str(sta.match(f).group())
                	else:
                        	starts[query][count_sname][subject_name] = starts[query][count_sname][subject_name] + ' ' + str(sta.match(f).group())
		except:
			print "starts %s %s %s" % (query, count_sname, subject_name)
                #print starts[count]

        #Extract lens: between lens and the next starts is dangerous
	#Becareful matching happens of pattern sta, so I add len(f) < 6
        if sta.match(f) and start_flag==False and len_flag==True and len(f) < 5:
                #print "successful %s" % sta.match(f).group(1)          
                #if count not in lens:
                #       lens[count] = ''
                #print "In lens %s %s %s" % (query, count_sname, subject_name)
		try:
                	if len(lens[query][count_sname][subject_name]) == 0:
                        	lens[query][count_sname][subject_name] = str(sta.match(f).group())
                	else:
                        	lens[query][count_sname][subject_name] = str(lens[query][count_sname][subject_name]) + ' ' + str(sta.match(f).group())
                except:
			print "lens %s %s %s" % (query, count_sname, subject_name)

        #Change flag
        if f.startswith("lens {"):
                start_flag = False
                len_flag = True
                continue
fn.close()

#write result
o.write("#Query\tsname\tstarts\tlens\n")

keys = starts.keys()
for q in keys:
	count = starts[q].keys()
	for s in count:
		ke = starts[q][s].keys()
		for k in ke:
			start1 = starts[q][s][k]
			lens1 = lens[q][s][k]
			o.write("%s\t%s\t%s\t%s\n" % (q, k, start1, lens1))






	
