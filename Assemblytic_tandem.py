#! /usr/bin/python2.7
import argparse
import gzip
import os
import pandas as pd
import re
import collections as cl
import sys
import getopt
import multiprocessing as mul

header='#reference	ref_start	ref_stop	ID	size	strand	type	ref_gap_size	query_gap_size	query_coordinates	method'.split()

def takeinput(thefile):
        
    table=pd.read_csv(thefile,sep='\t',index_col=False)
    
    if '#reference' not in [x for x in table.columns.tolist()]:
        table=pd.read_csv(thefile,sep='\t',names=header,index_col=None)
    

    table['index']=range(len(table))
        
    #table=table[(table['ref_gap_size']<0) | (table['query_gap_size']<0)]
    #table=table[(table['ref_gap_size']>0) & (table['query_gap_size']>0)]
    

    allcordinates=[[qcoordi.split(':')[0]]+map(int,qcoordi.split(':')[1].split('-'))+[ref,int(rend),int(rstart),index] for ref,rstart,rend, qcoordi,index in tuple(table[['#reference','ref_start','ref_stop', 'query_coordinates','index']].values.tolist())]
    

    
    return allcordinates


def find_intercept(list1,list2):
    
    list1=[sorted(x) for x in list1]
    list2=[sorted(x) for x in list2]
    
    
    l1=sum(list1,[])
    l2=sum([[x[0]-10,x[1]+10] for x in list2],[])
    
    len1=len(l1)
    
    lists=l1+l2
    
    lists_index=sorted(range(len(lists)), key= lambda x:lists[x])
    
    lists_sort=[lists[x] for x in lists_index]
    
    span_record={}
    current_spans=[]
    current_seq=[]
    for x in lists_index:
        
        
        if x>=len1:
            
            pair_index1=(x-len1)/2

            if x%2==0:
                current_spans.append(pair_index1)

                span_record[pair_index1]=[x for x in current_seq]
            else:
                current_spans=[x for x in current_spans if x != pair_index1]
        else:
            pair_index0=x/2
            
            if x%2==0:

                current_seq.append(pair_index0)
                for current_spans0 in current_spans:

                    span_record[current_spans0].append(pair_index0)
            else:

                current_seq=[x for x in current_seq if x != pair_index0]
            
    return span_record
    


def findoverlap(titles,allcordinates_list):
    
    titles_dict=cl.defaultdict(list)
    for x in titles:
        titles_dict[x[0]+'_'+x[3]].append(x)
    

    allcordinates=cl.defaultdict(list)
    for x in allcordinates_list:
        allcordinates[x[0]+"_"+x[3]].append(x)
    
    overlap_delta=[]
    overlap_sv=[]
    for contig,coordis in titles_dict.items():
        
        if contig not in allcordinates.keys():
            continue
        
        span_record=find_intercept( [x[4:6] for x in allcordinates[contig]], [x[4:6] for x in coordis])
                
        span_record={key:[x for x in index if max(allcordinates[contig][x][1:3]+coordis[key][1:3])- min(allcordinates[contig][x][1:3]+coordis[key][1:3]) - (max(coordis[key][1:3])-min(coordis[key][1:3]))-(max(allcordinates[contig][x][1:3])-min(allcordinates[contig][x][1:3]))<=10] for key, index in span_record.items()}
        
        overlap_delta.extend([x for i,x in enumerate(coordis) if len(span_record[i])>0])
        
        overlap_sv.extend([[allcordinates[contig][x0] for x0 in x] for x in span_record.values() if len(x)>0])


    return zip(overlap_delta,overlap_sv)


def update_coordi(rbreaks,qbreaks,rfind,qfind):
    
    
    print '\nout',rbreaks,qbreaks,rfind,qfind,'out\n'
    allrbreaks=rbreaks+rfind
    
    l=len(rbreaks)
    
    allrbreaks_index=sorted(range(len(allrbreaks)), key=lambda x:allrbreaks[x])
    if qbreaks[-1]>=qbreaks[0]:
        
        std=1
    else:
        std=-1
    
    
    allrbreaks_sort=[allrbreaks[x] for x in allrbreaks_index]
    
    find_results=cl.defaultdict(list)
    anchor_index=-1
    for index in allrbreaks_index:
        if index >=l:
            true_index=index-l

            find_results[true_index]=[anchor_index,std*(allrbreaks[index]-allrbreaks[anchor_index])]
        else:
            anchor_index=index
    
    
    find_results={key:qbreaks[x[0]]+x[1] for key,x in find_results.items() if x[0]>=0 and x[0]<len(rbreaks)}
    
    
    if len(find_results.keys())==1:
        #print rbreaks,qbreaks,rfind,qfind
        #exit()
        pass
    
    return find_results

    

def findcoordi(title, delta, coordi):
    
        
    print title, coordi
        
    
    delta=[int(x) for x in delta]
    
    print delta
    qstart,qend,rchr,rstart,rend=tuple(title[1:6])
    
    if rend<rstart:
        std='-'
        std0=-1
    else:
        std='+'
        std0=1
        
    print0=0
    if 1==2 and 16928 in [x for y in coordi for x in y] and 17850 in [x for y in coordi for x in y]:
        print0=1
    
    query_coordi=qstart
    ref_coordi=rstart
    
    query_coordis=[query_coordi]
    ref_coordis=[ref_coordi]
    
    for tick in delta:
        
        if tick>0:
            ref_change0=tick
            query_change0=tick-1
        else:
            ref_change0=-tick-1
            query_change0=-tick
            
            
        query_coordi+=std0*query_change0
        ref_coordi+=ref_change0
        
        query_coordis.append(query_coordi)
        ref_coordis.append(ref_coordi)
        
    rfind=[x for y in coordi  for x in y[4:6]]
    qfind=[x for y in coordi  for x in y[1:3]]

    if print0:
        print ref_coordis,query_coordis,rfind,qfind

    findresults=update_coordi(ref_coordis+[rend],query_coordis+[qend],rfind,qfind)
    
    if len([x[-1] for x in coordi if x[-1] in [65,4152,6703,6823,7828,8244,10308,12597,12966,12981,13032,13046,13053,13061]]):
        print [x[-1] for x in coordi],(ref_coordis,query_coordis,rfind,qfind), title,coordi
    
    print findresults

    organized_results={}
    for index, values in findresults.items():
        
        sv_index=coordi[index/2][-1]
        
        
        if sv_index not in organized_results.keys():
            organized_results[sv_index]=[[],[]]
        
        organized_results[sv_index][index%2]=values
    
    print organized_results
    
    return organized_results
        
            
        
                


def run(args):
    deltafile = args.delta
    inputfile = args.inputfile
    outputfile = args.outputfile
    if outputfile=='':
        outputfile=inputfile+'_tandem_fixed.tsv'
    #minimum_variant_size=1
    
    allcordinates=takeinput(inputfile)
    
    with open(deltafile, mode='r') as f:
        delta_read=f.read().split(">")[1:]
    f.close()
    
    
    
    
    titles=[x.splitlines()[0].split()+[i] for i,x in enumerate(delta_read)]
    
    
    
        
    delta_read=[('\n'.join(x.splitlines()[1:])).split('\n0\n') for x in delta_read]
    
    
    
    delta_segments=[[segment.splitlines()[0].split()+[j] for j,segment in enumerate(contig)] for contig in delta_read]
  

    if 1==1:
        titles=[titles[14]]
        delta_segments=[delta_segments[14]]


  
    titles=[[title[1],int(segment[2]), int(segment[3]), title[0],int(segment[0]),int(segment[1]),title[-1], segment[-1]] for title,contig in zip(titles,delta_segments) for segment in contig]
    


    
    
    delta_read=[[segment.splitlines()[1:] for j,segment in enumerate(contig)] for contig in delta_read]
    

    
    overlap_groups=findoverlap(titles,allcordinates)
    
    
    
    
    
    print 'groups',overlap_groups
    
    
    full_results=cl.defaultdict(list)
    for delta,svs in  overlap_groups:
        
        i0,j0=delta[-2:]
        
        print 'delta',delta
        
        results=findcoordi(delta, delta_read[i0][j0], svs)
        
        for k,v in results.items():
            
            full_results[k].append('-'.join([str(x) for x in v]))
    
    

    
    table=pd.read_csv(inputfile,sep='\t')
    
    if '#reference' not in [x for x in table.columns.tolist()]:
        table=pd.read_csv(inputfile,sep='\t',names=header)
    
    output=[]
    for i in xrange(len(table)):
        
        if i not in  full_results.keys():
            output.append('.')
        
        else:
            output.append(';'.join([str(x) for x in full_results[i]]))
            
    table['tandem_fix']=output
    print table[table['type'].str.contains('Tandem')]
    table.to_csv(outputfile,index=True, mode='w' , sep='\t')
    
    print table
    exit()
    

def main():

    parser=argparse.ArgumentParser(description="Calculate tandem coordinates")
    parser.add_argument("-d","--delta",help="delta file" ,dest="delta", type=str, required=True)
    parser.add_argument("-i","--input",help="input sv file" ,dest="inputfile",type=str,required=True)
    parser.add_argument("-o","--output",help="output sv file" ,dest="outputfile",type=str,default='')
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()




