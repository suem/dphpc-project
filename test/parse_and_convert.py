import argparse
import sys

parser = argparse.ArgumentParser(description='Parse graph')
parser.add_argument("--directed", 
        help="input is a directed graph", action="store_true")
args = parser.parse_args()

counter = 0
unique_ids = {}
def getNodeId(node_id):
    global counter
    if node_id in unique_ids:
        return unique_ids[node_id]
    unique_ids[node_id] = counter
    counter += 2
    return unique_ids[node_id]



for line in sys.stdin:
    if line[0] == "#":
        continue

    [u,v] = line.split()
    u = int(u)
    v = int(v)

    u_left = getNodeId(u)
    u_right = u_left + 1

    v_left = getNodeId(v)
    v_right = v_left + 1

    if args.directed:
        print("%d %d" % (u_left, v_right)) 
    else:
        print("%d %d" % (u_left, v_right)) 
        print("%d %d" % (u_right, v_left)) 

        


