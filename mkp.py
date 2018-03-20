from gurobipy import *

f = open("gk11.dat", "r")
data = f.read()

words = data.split()

num_vars = int(words[0])
num_consts = int(words[1])
# words[2] : ignore

profits = [ float(i) for i in words[3:(num_vars+1)] ]
print profits
end_consts = num_vars+1+num_consts*num_vars+1

constraints = [ float(i) for i in words[(num_vars+1):end_consts] ]
capacities = [ float(i) for i in words[end_consts:(end_consts+num_consts+1)] ]


try:
    m = Model("mkp")

    x = m.addVars(num_vars, lb=0.0, ub=1.0, obj=0.0, vtype=GRB.CONTINUOUS)
    
    

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError:
    print('Encountered an attribute error')

