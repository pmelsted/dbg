import collections, sys


def parseGFA(f):
    rev = dict(('+-','-+'))
    nodes = []
    G = collections.defaultdict(dict)
    for line in f:
        if len(line) == 0:
            pass
        c = line[0]
        if c == 'H':
            pass
        elif c == 'S':
            l = line.split()
            nodes.append(l[1])
        elif c == 'L':
            l = line.split()
            i,j = l[1],l[3]
            o = (l[2],l[4])
            if j not in G[i]:
                G[i][j] = o
            if i not in G[j]:
                G[j][i] = rev[o[1]],rev[o[0]]
    
    for i in nodes:
        if i not in G:
            d = G[i]
    
    return G

def dotfromGFA(G):
    lines = ["digraph G {", "graph [rankdir=LR, fontname=\"Courier\"];", "node [shape=record];"]
    
    for i in G:
        lines.append("%s[label=\"<F> %s+ | <R> %s-\"];"%(i,i,i))

    done = set()

    t = dict(('+F','-R'))
    for i in G:
        for j in G[i]:
            o = G[i][j]
            lines.append("%s:%s -> %s:%s;"%(i,t[o[0]],j,t[o[1]]))


    lines.append("}")
    return '\n'.join(lines);



if __name__ == "__main__":
    f = sys.stdin
    if len(sys.argv) > 1:
        f = open(sys.argv[1])

    G = parseGFA(f)
    print dotfromGFA(G)
