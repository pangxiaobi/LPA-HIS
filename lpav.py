import collections
import random
from graph_helper import load_graph
from graph_helper import clone_graph
import string

'''
paper : <<Near linear time algorithm to detect community structures in large-scale networks>>
'''
neighbor={}
label={}
def readdata(path):
    fin=open(path,'r')
    i=0
    for line in fin.readlines():
        if i==0:
            i+=1
            continue
        line=line.rstrip().split(' ')
        s=string.atoi(line[0])
        t=string.atoi(line[1])
        if not neighbor.has_key(s):
            tmp=[]
            neighbor[s]=tmp[:]
        if not neighbor.has_key(t):
            tmp=[]
            neighbor[t]=tmp[:]
        neighbor[s].append(t)
        neighbor[t].append(s)
    return neighbor
def can_stop():
    # all node has the label same with its most neighbor
    for (i,n) in neighbor.items():
        l = label[i]
        max_labels = get_max_neighbor_label(n)
        if(l not in max_labels):
            return False
    return True
    
def get_max_neighbor_label(nei):
    m ={}
    for i in nei:
        if m.has_key(label[i]):
            m[label[i]] += 1
        else:
            m[label[i]]=1
    max_v = max(m.values())
    return [item[0] for item in m.items() if item[1] == max_v]


'''asynchronous update'''
def populate_label():
    #random visit
    visitSequence = random.sample(neighbor.keys(),len(neighbor))
    for i in visitSequence:
        l=label[i]
        max_labels = get_max_neighbor_label(neighbor[i])
        if(l not in max_labels):
            newLabel = random.choice(max_labels)
            label[i] = newLabel
    
def get_communities():
    communities ={}
    for node in neighbor.keys():
        l=label[node]
        if not communities.has_key(l):
            tmp=[]
            communities[l]=tmp[:]
        communities[l].append(node)
    return communities.values()

def execute(path):
    #initial label
    readdata(path)
    for i in neighbor.keys():
        label[i] = i
    iter_time = 0
    #populate label
    while(not can_stop() and iter_time<10):
        populate_label()
        iter_time += 1
    return get_communities()
    
    
def main():
    fout=open(r'comunity_big.txt','w')
    print "haha"
    communities = execute(r'edge_big.txt')
    print len(communities)
    sum=0
    for community in communities:
        for term in community:
            fout.writelines (str(term)+' ')
        fout.writelines('\n')
main()