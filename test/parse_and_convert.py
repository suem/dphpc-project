import argparse
import sys

parser = argparse.ArgumentParser(description='Parse graph')
parser.add_argument("--directed", 
        help="input is a directed graph", action="store_true")
args = parser.parse_args()

# counter = 0
# unique_ids = {}
# def getNodeId(node_id):
#     global counter
#     if node_id in unique_ids:
#         return unique_ids[node_id]
#     unique_ids[node_id] = counter
#     counter += 2
#     return unique_ids[node_id]


# left_ids = {}
# counter_left = 0
# def getLeftId(node_id):
#     global counter_left
#     if node_id in left_ids:
#         return left_ids[node_id]
#     left_ids[node_id] = counter
#     counter += 1
#     return left_ids[node_id]


# right_ids = {}
# counter_right = 0


class IdGenerator(object):

    def __init__(self):
        self.counter = 0
        self.unique_ids = {}

    def getNodeId(self, node_id):
        if node_id in self.unique_ids:
            return self.unique_ids[node_id]
        self.unique_ids[node_id] = self.counter
        self.counter += 1
        return self.unique_ids[node_id]

    def getCounter(self):
        return self.counter


left_generator = IdGenerator()
right_generator = IdGenerator()

edges = []

for line in sys.stdin:
    if line[0] == "#":
        continue

    [u,v] = line.split()
    u = int(u)
    v = int(v)

    u_left = left_generator.getNodeId(u)
    u_right = right_generator.getNodeId(u)

    v_left = left_generator.getNodeId(v)
    v_right = right_generator.getNodeId(v)

    if args.directed:
        edges += [(u_left, v_right)] 
    else:
        edges += [(u_left, v_right)] 
        edges += [(v_left, u_right)] 


n_left = left_generator.getCounter() 
n_right = right_generator.getCounter()
n_total = n_left + n_right

print(n_total) # total number of edges
print(n_left)  # starting id of right edges

for (l, r) in edges:
    # print("%d %d" % (l, r)) 
    print("%d %d" % (l, r + n_left)) 

