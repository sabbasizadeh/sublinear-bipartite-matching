import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import random
import time

eps = 0.05
c2=3
rep=30

# arXiv cond-mat
r=22015
n=16726
reverse=False
fname="out.opsahl-collaboration"

#youtube
#r=94238
#n=30087
#reverse=True
#fname="out.youtube-groupmemberships"

#actor RANKING_MOvie
#r=383640
#n=127823
#reverse=False
#fname="out.actor-RANKING_MOvie"

#stackoverflow
#r=545196
#n=96678
#reverse=True
#fname="out.stackexchange-stackoverflow"

#flickr
#r=395979
#n=103631
#reverse=True
#fname="out.flickr-groupmemberships"

S = []
greedy_S = []
pi = []
sigma = []
computed_nodes = [-1]*(n+r)

greedy_visited_edges=0
ranking_visited_edges=0

def generate_graph():
    global pi, sigma, S, greedy_S
    global greedy_visited_edges, ranking_visited_edges
    global computed_nodes 
    computed_nodes = [-1]*(n+r)
    greedy_visited_edges=0
    ranking_visited_edges=0
    S = random.sample(range(n), c2*int((1/eps)*(1/eps)))
    greedy_S = random.sample(range(n+r), c2*int((1/eps)*(1/eps)))
    pi = [i for i in range(n)]
    sigma = [i for i in range(n,n+r)]
    random.shuffle(pi)
    random.shuffle(sigma)
    sigma = [0]*n+sigma

    with open(fname) as f:
        lines = [line.strip().split()[0:2] for line in f]
    edges=[]
    if reverse is True:
        for edge in lines:
            edges.append((int(edge[1])-1,int(edge[0])-1+n))
    else:
        for edge in lines:
            edges.append((int(edge[0])-1,int(edge[1])-1+n))
    G = nx.Graph()
    G.add_nodes_from([i for i in range(n)], bipartite=0)
    G.add_nodes_from([i for i in range(n,n+r)], bipartite=1)
    G.add_edges_from(edges)
    v_L = [i for i in range(n)]
    v_R = [i for i in range(n,n+r)]
    edge_pi = [i for i in range(len(G.edges))]
    random.shuffle(edge_pi)
    global greedy_computed_edges
    greedy_computed_edges = [-1]*len(G.edges)
    i=0
    for e in G.edges(data=True):
        G[e[0]][e[1]]['rank'] = edge_pi[i]
        G[e[1]][e[0]]['rank'] = edge_pi[i]
        i = i+1
    return G, v_L, v_R

def hopcroft_karp_matching(G, v_L):
    print("started to run matching")
    matching = nx.bipartite.hopcroft_karp_matching(G, v_L)
    print(len(matching)/2)


def RANKING_get_sorted_neighbors(G, s, perm, bound):
    neigh = []
    for t in G.neighbors(s):
        try:
            if perm[t]<bound:
                neigh.append([t,perm[t]])
        except:
            print(t)
            print(bound)
            print(perm)
            print(s)
    neigh = sorted(neigh, key=lambda tup: tup[1])
    return [x[0] for x in neigh]

def RANKING_matching_estimator(G, v_L, v_R):
    print("started to run RANKING matching")
    x = 0
    for s in S:
        res = RANKING_MO(s, float('inf'), G, v_L, v_R)
        if res==0:
            computed_nodes[s]=-2
        if res==1:
            x = x + 1
    estimated_size = n*x/len(S)
    print("visited", ranking_visited_edges)
    return estimated_size

def RANKING_MO(s, rank_sigma, G, v_L, v_R):
    global ranking_visited_edges
    if computed_nodes[s]==-2:
        return 0
    if computed_nodes[s] > -1 and computed_nodes[s] < rank_sigma:
        #ranking_visited_edges=ranking_visited_edges+1
        return 1
    if computed_nodes[s] > -1 and computed_nodes[s] >= rank_sigma:
        return 0
    neighbors = RANKING_get_sorted_neighbors(G,s,sigma,rank_sigma)
    if len(neighbors)==0:
        return 0
    if len(neighbors) >= len(v_L):
        computed_nodes[s]=1
        return 1
    for j in neighbors:
        ranking_visited_edges = ranking_visited_edges+1
        if computed_nodes[j] > -1 and computed_nodes[j] == pi[s]:
            return 1
        if computed_nodes[j] > -1 and computed_nodes[j] < pi[s]:
            continue
        K = RANKING_get_sorted_neighbors(G, j, pi, pi[s])
        if len(K)==0:
            computed_nodes[s]=sigma[j]
            computed_nodes[j]=pi[s]
            return 1
        evaluated = 0
        for k in K:
            ranking_visited_edges = ranking_visited_edges+1
            if RANKING_MO(k, sigma[j], G, v_L, v_R)==0:
                break
            else:
                evaluated = evaluated+1
        if evaluated == len(K):
            computed_nodes[j]=pi[s]
            computed_nodes[s]=sigma[j]
            return 1
    return 0

def Greedy_node_get_sorted_neighbors(G, s, bound):
    neigh = []
    for t in G.neighbors(s):
        if G[s][t]['rank']<bound:
            neigh.append([t, G[s][t]['rank']])
    neigh = sorted(neigh, key=lambda tup: tup[1])
    return [x[0] for x in neigh]  

def Greedy_edge_get_sorted_neighbors(G, e, bound):
    neigh = []
    for t in G.neighbors(e[0]):
        if G[e[0]][t]['rank']<bound:
            neigh.append([(e[0],t), G[e[0]][t]['rank']])
    for t in G.neighbors(e[1]):
        if G[e[1]][t]['rank']<bound:
            neigh.append([(e[1],t), G[e[1]][t]['rank']])
    neigh = sorted(neigh, key=lambda tup: tup[1])
    return [x[0] for x in neigh]  
    
def Greedy_VO(v, G):
    global greedy_visited_edges
    neighbors = Greedy_node_get_sorted_neighbors(G, v, float('inf'))
    for j in neighbors:
        greedy_visited_edges=greedy_visited_edges+1
        if Greedy_MO(G, (v,j)) == 1:
            return 1
    return 0
    
def Greedy_MO(G, e):
    global greedy_visited_edges
    if greedy_computed_edges[G[e[0]][e[1]]['rank']] > -1:
        return greedy_computed_edges[G[e[0]][e[1]]['rank']]
    neighbors = Greedy_edge_get_sorted_neighbors(G,e, G[e[0]][e[1]]['rank'])
    for e1 in neighbors:
        greedy_visited_edges=greedy_visited_edges+1
        if Greedy_MO(G, e1)==1:
            greedy_computed_edges[G[e[0]][e[1]]['rank']]=0
            return 0
    greedy_computed_edges[G[e[0]][e[1]]['rank']]=1
    return 1

def Greedy_matching_estimator(G):
    print("started to run Greedy matching")
    x = 0
    for s in greedy_S:
        if Greedy_VO(s, G)==1:
            x = x + 1
    estimated_size = n*x/len(S)
    return estimated_size

def main():
    G, v_L, v_R = generate_graph()
    start = time.time()
    hopcroft_karp_matching(G, v_L)
    end = time.time()
    print("hopcroft_karp time: ", end-start)
    
    rank_alg_visited=[]
    rank_alg_estimate=[]
    
    greedy_alg_visited=[]
    greedy_alg_estimate=[]
    
    for i in range(rep):
        print(i)
        G, v_L, v_R = generate_graph()
        
        start = time.time()
        est_size = RANKING_matching_estimator(G, v_L, v_R)
        end = time.time()
        print(est_size)
        print(ranking_visited_edges)
        rank_alg_visited.append(ranking_visited_edges)
        rank_alg_estimate.append(est_size)

        start = time.time()
        est_size = Greedy_matching_estimator(G)
        end = time.time()
        print(est_size)
        print(greedy_visited_edges)
        greedy_alg_visited.append(greedy_visited_edges)
        greedy_alg_estimate.append(est_size)

    print(greedy_alg_visited)
    print(rank_alg_visited)
    
    print(greedy_alg_estimate)
    print(rank_alg_estimate)
    
    print('${0} \pm {1}$'.format(round(np.mean(greedy_alg_estimate),1), round(np.std(greedy_alg_estimate),1)))
    print('${0} \pm {1}$'.format(round(np.mean(rank_alg_estimate),1), round(np.std(rank_alg_estimate),1)))

    print('${0} \pm {1}$'.format(round(np.mean(greedy_alg_visited),1), round(np.std(greedy_alg_visited),1)))
    print('${0} \pm {1}$'.format(round(np.mean(rank_alg_visited),1), round(np.std(rank_alg_visited),1)))
    
    print(max(greedy_alg_estimate))
    print(max(rank_alg_estimate))
    
main()