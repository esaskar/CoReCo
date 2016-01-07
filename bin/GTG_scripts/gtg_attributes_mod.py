# grep gtg GS_GI.mapping | cut -f 3,6-7,13-16 | time python gtg_attributes.py | cut -f 1-3,7 | time perl parse_gtg_attributes.pl 

import sys,os,string,re
import array

###############################################################################

def getCAA(leaf_vertex):
        # get CAA_id and aaix for leaf_vertex
        #indexArray=array.array('L')
        indexArray=array.array('I')
        if param_arch==64: 
                address=leaf_vertex*8*2 # 64-bit
	else: 
                address=leaf_vertex*8
        #print "getCAA", address, fh_CAA_vertex_ptrin
        try: 
                fh_CAA_vertex_ptrin.seek(address)
        except: 
                return((0,0)) # vertex not in database
        indexArray.fromfile(fh_CAA_vertex_ptrin,2)
        if param_arch==64: data_address=indexArray[0]*4*2 # 64-bit
	else: data_address=indexArray[0]*4

        #print leaf_vertex, data_address
        #sys.exit()

        try: fh_CAA_vertex_datain.seek(data_address)
        except: return((0,0)) # address not in database
        #dataArray=array.array('L')
        dataArray=array.array('I')
        dataArray.fromfile(fh_CAA_vertex_datain,indexArray[1]) # list of CAAs containing query vertex
	#print "# leaf_vertex, dataArray",leaf_vertex,dataArray
        # return list of CAAs with aaix
        tuplelist=[]
        for CAA_id in dataArray:
                # lookup aaix of CAA_id
                aaix=get_aaix_for_CAAid(CAA_id)
                tuplelist.append( (CAA_id, aaix) )
        return(tuplelist) # CAA_id, aaix

def get_aaix_for_CAAid(CAA_id):
                # get all members of CAA
                #indexArray=array.array('L')
                indexArray=array.array('I')
                if param_arch==64: address=(CAA_id*12-12)*2 # 64-bit
		else: address=CAA_id*12-12
                #print "# CAA_id=%i address=%i" %(CAA_id,address)
                try: fh_CAA_member_ptrin.seek(address)
                except:
                        #print "CAA_id not in database:",CAA_id
                        return(0)
                indexArray.fromfile(fh_CAA_member_ptrin,3)
                #print "query CAA_id retrieved", CAA_id,indexArray
                return(indexArray[2])

###############################################################################

class MyGetSequence: # global fh_nid_ptrin/datain
	"""emulates gtglib's g.getSequence object

	x=MyGetSequence(nid) -> x.front, x.back, x.aasequence()
	
	"""
	def __init__(self,nid):
		if param_arch==64: address=(nid-1)*16*2 # 64-bit
		else: address=(nid-1)*16
		fh_nid_ptrin.seek(address)
		#indexArray=array.array('L') 
                indexArray=array.array('I') 
                #print "arraysize:",indexArray.itemsize

		indexArray.fromfile(fh_nid_ptrin,4)
		self.front=indexArray[0]
		self.back=indexArray[1]
		self.seqstart=indexArray[2]
		self.nres=indexArray[3]
		#print "MyGetSequence(%i): " %nid,indexArray
		
	def aasequence(self):
		fh_nid_datain.seek(self.seqstart)
		seqArray=array.array('c')
		seqArray.fromfile(fh_nid_datain,self.nres)
		aasequence=seqArray.tostring()
		return(aasequence)

class MyVertex: 
	"""
	x=MyVertex(leaf vertex) -> x.asClusterIndex, x.asNid, x.asResidue, x.arrayElement (points to cluster topology table) 
	"""
	def __init__(self,vertex):
		if param_arch==64: address=vertex*16*2 # 64-bit
		else: address=vertex*16
		fh_vertex_ptrin.seek(address)
		#indexArray=array.array('L')
                indexArray=array.array('I')
		try: indexArray.fromfile(fh_vertex_ptrin,4)
		except: 
#			print "#MyVertex error (MyVertex)",address,vertex	# 137727259, 2203636144
			indexArray.fromlist([0,0,0,0]) # HACK
		self.asClusterIndex=indexArray[0]
		self.arrayElement=indexArray[1]
		self.asNid=indexArray[2]
		self.asResidue=indexArray[3]

def input_one_line():
#	line=raw_input("Query,v4_nid,query_starts,sbjct_starts,blocklengths,v2_rep40,query_starts,sbjct_starts,blocklengths,query_sequence:")
#	print
	line=raw_input()
	return(line)

def get_alignment(query_starts,sbjct_starts,block_length):
	# ali1[query_residue]=sbjct_residue
	# ali2[sbjct_residue]=query_residue
	ali1={}
	ali2={}
	n=len(block_length)
	for i in range(0,n):
		q=query_starts[i]
		s=sbjct_starts[i]
		for k in range(0,block_length[i]):
			q+=1
			s+=1
			ali1[q]=s
			ali2[s]=q
	return(ali1,ali2)

def get_attributes_intracluster(query_name,gtg_nid,ali2,query_sequence):
	gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
        #print "# gtg_nid=%i front=%i back=%i" %(gtg_nid,front,back)
        for gtg_vertex in range (front,back):
		x=MyVertex(gtg_vertex)
		sbjct_residue=x.asResidue
		if not ali2.has_key(sbjct_residue): continue # not aligned
		cluster_id=x.asClusterIndex
		array_element=x.arrayElement
		query_residue=ali2[sbjct_residue]
		print "%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,gtg_vertex,cluster_id,array_element)

def get_attributes(query_name,gtg_nid,ali2,query_sequence):
	gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
	#print "get_attributes: gtg_nid=%i front=%i back=%i nres=%i" %(gtg_nid,front,back,len(query_sequence))
	nres=len(query_sequence)
        for gtg_vertex in range (front,back):
		x=MyVertex(gtg_vertex)
		sbjct_residue=x.asResidue
		if not ali2.has_key(sbjct_residue): continue # not aligned
		query_residue=ali2[sbjct_residue]
		if query_residue>nres: break # security hack
		try: query_aaix=aa2ix[query_sequence[query_residue-1]]
		except: break # hack: mysql error 
		tuplelist=getCAA(gtg_vertex) # match CAA-aaix
		if len(tuplelist)>param_exclude_CAAs_larger: continue # exclude oversized CAAs
		for CAA_vertex,CAA_aaix in tuplelist:
			#print "# ",CAA_vertex,CAA_aaix,gtg_vertex,query_aaix,sbjct_residue,query_residue
			if CAA_aaix<>query_aaix: continue
			print "%s\t%i\t%s\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,CAA_aaix,CAA_vertex)

def get_simple_attributes(query_name,gtg_nid,ali2,query_sequence):
	# attribute == clusid.aaix
        gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
        for gtg_vertex in range (front,back):
                x=MyVertex(gtg_vertex)
                sbjct_residue=x.asResidue
                if not ali2.has_key(sbjct_residue): continue # not aligned
                query_residue=ali2[sbjct_residue]
                query_aaix=aa2ix[query_sequence[query_residue-1]]
		clusid=x.asClusterIndex
		print "%s\t%i\t%s\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,query_aaix,clusid)

param_exclude_CAAs_larger = 10000

def initGTG(datadir, CAAlevel):
        global param_arch, DATADIR, param_exclude_CAAs_larger
        global aa2ix
        global fh_nid_ptrin, fh_nid_datain
        global fh_vertex_ptrin
        global fh_CAA_vertex_ptrin, fh_CAA_vertex_datain
        global fh_CAA_member_ptrin, fh_CAA_member_datain
        global param_CAA_level
        param_CAA_level = CAAlevel
        param_arch = 32
        DATADIR = datadir + "/"
	aa2ix={'A':1, 'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20,'X':0, 'Z': 0, 'B': 0, 'a':2,'b':2,'c':2,'d':2,'e':2,'f':2,'g':2,'h':2,'i':2,'j':2,'k':2,'l':2,'m':2,'n':2,'o':2,'p':2,'q':2,'r':2,'s':2,'t':2,'u':2,'v':2,'w':2,'x':2,'y':2,'z':2}
        fh_nid_ptrin=open(DATADIR+'nidindex.ptr','rb')
	fh_nid_datain=open(DATADIR+'nidindex.store','rb')
	fh_vertex_ptrin=open(DATADIR+'vertex.ptr','rb') 
        if param_CAA_level==0:
                fh_CAA_vertex_datain=open(DATADIR+"CAA0_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA0_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA0_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA0_members.ptr","rb")
        elif param_CAA_level==1:
                fh_CAA_vertex_datain=open(DATADIR+"CAA1_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA1_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA1_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA1_members.ptr","rb")
        elif param_CAA_level==2:
                fh_CAA_vertex_datain=open(DATADIR+"CAA2_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA2_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA2_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA2_members.ptr","rb")
        

if __name__ == "__main__":
	param_arch=32 #64

	try:
		DATADIR=sys.argv[2] # '/data/liisa/'
		DATADIR+='/' # just in case
		param_CAA_level=int(sys.argv[1]) # 0 or 1 or 2 or 3 (3=simple clusid.aaix attributes)
	except:
                print "USAGE: python gtg_attributes.py <CAA_level> <DATADIR>"
                sys.exit(1)

	param_exclude_CAAs_larger=10000

	aa2ix={'A':1, 'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20,'X':0, 'Z': 0, 'B': 0,
		'a':2,'b':2,'c':2,'d':2,'e':2,'f':2,'g':2,'h':2,'i':2,'j':2,'k':2,'l':2,'m':2,'n':2,'o':2,'p':2,'q':2,'r':2,'s':2,'t':2,'u':2,'v':2,'w':2,'x':2,'y':2,'z':2}
	
	# connect databases
	#print "DATADIR=%s" %DATADIR
	fh_nid_ptrin=open(DATADIR+'nidindex.ptr','rb')
	fh_nid_datain=open(DATADIR+'nidindex.store','rb')
	fh_vertex_ptrin=open(DATADIR+'vertex.ptr','rb') 

	# CAA databases
        if param_CAA_level==0:
                fh_CAA_vertex_datain=open(DATADIR+"CAA0_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA0_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA0_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA0_members.ptr","rb")
        elif param_CAA_level==1:
                fh_CAA_vertex_datain=open(DATADIR+"CAA1_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA1_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA1_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA1_members.ptr","rb")
        elif param_CAA_level==2:
                fh_CAA_vertex_datain=open(DATADIR+"CAA2_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA2_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA2_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA2_members.ptr","rb")

	# input: grep 'gtg' from pairsdb_mapper.pl-output
	while(1):
		try:
			line=input_one_line()
		except:
			#print " Exit"
			break

		data=line.split("\t")
		#print "data: ",data
		query_name=data[0]
		try: gtg_nid=int(data[2])
		except: continue
		#print "gtg_nid: ",gtg_nid
		if gtg_nid<1: continue # no mapping
		try:
			query_starts=map(int,data[3].split(None))
			sbjct_starts=map(int,data[4].split(None))
			block_lengths=map(int,data[5].split(None))
			query_sequence=data[6]
			#print "query_starts: ",query_starts
			#print "sbjct_starts: ",sbjct_starts
			#print "block_lengths: ",block_lengths
		except:
			continue
		ali1,ali2=get_alignment(query_starts,sbjct_starts,block_lengths)
	        #print "ali1: ",ali1
		#print "ali2: ",ali2
		if param_CAA_level<3: get_attributes(query_name,gtg_nid,ali2,query_sequence)
		else: get_simple_attributes(query_name,gtg_nid,ali2,query_sequence)
	
	if param_CAA_level<3:
		fh_nid_ptrin.close()
		fh_nid_datain.close()
		fh_vertex_ptrin.close()
	        fh_CAA_vertex_datain.close()
	        fh_CAA_vertex_ptrin.close()
	        fh_CAA_member_datain.close()
	        fh_CAA_member_ptrin.close()

