from gurobipy import *

num_voxels = 25 # range = 0-24
range_voxels = range(num_voxels)

tumor_voxels = [7, 8, 12, 13, 17, 21, 22] # Vt
organs_risk_voxels = [6, 10, 11, 15, 16, 20] # Vor
normal_voxels = [0, 1, 2, 3, 4, 5, 9, 14, 18, 19, 23, 24] # Vs

angle_index = { 0 : 0, 90 : 1, 180 : 2, 270 : 3 }

selected_angles = [0, 180]
selected_angles_indices = [angle_index[i] for i in selected_angles]

num_beamlets = 6
range_beamlets = range(num_beamlets)

influence_dose_matrix = [] # matrix A, indexed by beams x beamlets x voxels
                           # In this case: 4 x 6 x 25

A_0   = [ [], [], [], [], [], [] ]
A_90  = [ [], [], [], [], [], [] ]
A_180 = [ [], [], [], [], [], [] ]
A_270 = [ [], [], [], [], [], [] ]

A_0[0] = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.8, 0.8, 0.7, 0.8, 0.8, 0.85, 0.85, 0.75, 0.8, 0.75, 0.8, 0.7, 0.7, 0.65, 0.7, 0.7, 0.5, 0.5]
A_0[1] = A_0[0]
A_0[2] = A_0[0]
A_0[3] = A_0[0]
A_0[4] = A_0[0]
A_0[5] = A_0[0]

A_90[0] = [0.6, 0.7, 0.75, 0.7, 0.6, 0.6, 0.7, 0.8, 0.75, 0.5, 0.6, 0.7, 0.8, 0.75, 0.6, 0.6, 0.7, 0.8, 0.7, 0.6, 0.6, 0.75, 0.8, 0.7, 0.6]
A_90[1] = A_90[0]
A_90[2] = A_90[0]
A_90[3] = A_90[0]
A_90[4] = A_90[0]
A_90[5] = A_90[0]

A_180[0] = [0.6, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.7, 0.8, 0.8, 0.85, 0.85, 0.7, 0.7, 0.7, 0.8, 0.7, 0.7, 0.65, 0.7, 0.7, 0.6, 0.6]
A_180[1] = A_180[0]
A_180[2] = A_180[0]
A_180[3] = A_180[0]
A_180[4] = A_180[0]
A_180[5] = A_180[0]

A_270[0] = [0.6, 0.65, 0.75, 0.7, 0.6, 0.5, 0.7, 0.8, 0.8, 0.6, 0.6, 0.7, 0.8, 0.8, 0.6, 0.6, 0.7, 0.8, 0.7, 0.6, 0.65, 0.75, 0.8, 0.7, 0.6]
A_270[1] = A_270[0]
A_270[2] = A_270[0]
A_270[3] = A_270[0]
A_270[4] = A_270[0]
A_270[5] = A_270[0]

influence_dose_matrix.append(A_0  )
influence_dose_matrix.append(A_90 )
influence_dose_matrix.append(A_180)
influence_dose_matrix.append(A_270)

# Next, the arrays that indicate the indexes of non-zero elements of A

non_zero_tumor_voxels = tumor_voxels
non_zero_or_voxels = organs_risk_voxels
non_zero_normal_voxels = normal_voxels

Sor = 30 # upper bound on organs at risk
Ss = 10 # upper bound on healthy tissue
D = 50 # dose prescription for tumor

importance_factors = [0.2, 0.1, 0.3, 0.4]

try:

    m = Model("Obal_model")

    # Beamlets dose variables
    x = m.addVars(len(angle_index), num_beamlets, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="beamlet_dose")

    """
      Deviation of prescribed dose variables.
      theta:   deviation on organs at risk
      delta:   deviation on normal tissue
      epsilon: deviation on tumor
    """    
    theta_plus    = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="theta_plus"   )
    theta_minus   = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="theta_minus"  )
    delta_plus    = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="delta_plus"   )
    delta_minus   = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="delta_minus"  )
    epsilon_plus  = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="epsilon_plus" )
    epsilon_minus = m.addVars(num_voxels, lb=0.0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="epsilon_minus")

    m.addConstrs( (quicksum(quicksum(influence_dose_matrix[b][q][v] * x[b,q] for b in selected_angles_indices) for q in range_beamlets) == Sor + theta_plus[v]   - theta_minus[v]   for v in non_zero_or_voxels    ), "deviation_or")
    m.addConstrs( (quicksum(quicksum(influence_dose_matrix[b][q][v] * x[b,q] for b in selected_angles_indices) for q in range_beamlets) == Ss  + delta_plus[v]   - delta_minus[v]   for v in non_zero_normal_voxels), "deviation_s" )
    m.addConstrs( (quicksum(quicksum(influence_dose_matrix[b][q][v] * x[b,q] for b in selected_angles_indices) for q in range_beamlets) == D   + epsilon_plus[v] - epsilon_minus[v] for v in non_zero_tumor_voxels ), "deviation_t" )

    """
    f1 = m.addVars(lb = 0, ub = GRB.INFINITY, obj = 0, vtype = GRB.CONTINUOUS, name="f1")
    f2 = m.addVars(lb = 0, ub = GRB.INFINITY, obj = 0, vtype = GRB.CONTINUOUS, name="f2")
    f3 = m.addVars(lb = 0, ub = GRB.INFINITY, obj = 0, vtype = GRB.CONTINUOUS, name="f3")
    f4 = m.addVars(lb = 0, ub = GRB.INFINITY, obj = 0, vtype = GRB.CONTINUOUS, name="f4")
    """
    f = m.addVars(range(4), lb = 0, ub = GRB.INFINITY, obj = 0, vtype = GRB.CONTINUOUS, name="f")

    m.addConstr(f[0] == quicksum(theta_plus[i]    for i in range_voxels) )
    m.addConstr(f[1] == quicksum(delta_plus[i]    for i in range_voxels) )
    m.addConstr(f[2] == quicksum(epsilon_plus[i]  for i in range_voxels) )
    m.addConstr(f[3] == quicksum(epsilon_minus[i] for i in range_voxels) )

    g = LinExpr()
    for i in range(4):
        g.addTerms(importance_factors[i], f[i])

    m.setObjective(g, GRB.MINIMIZE)

    m.optimize()

    print(m.getObjective().getValue())
    for i in range(4):
        print(f[i])

    # Compute the deposited dose in each voxel

    dose_on_voxels = [LinExpr() for i in range(num_voxels)]
        
    for v in tumor_voxels:
        for b in selected_angles_indices:
            for q in range_beamlets:
                dose_on_voxels[v] += influence_dose_matrix[b][q][v] * x[b,q]


    for v in organs_risk_voxels:
        for b in selected_angles_indices:
            for q in range_beamlets:
                dose_on_voxels[v] += influence_dose_matrix[b][q][v] * x[b,q]


    for v in normal_voxels:
        for b in selected_angles_indices:
            for q in range_beamlets:
                dose_on_voxels[v] += influence_dose_matrix[b][q][v] * x[b,q]

    print("Dose on voxels:")
    i = 0
    for v in dose_on_voxels:
        print(i, v.getValue())
        i+=1

    print("\n\n")
    for b in selected_angles_indices:
            for q in range_beamlets:
                print(x[b,q].getValue())


except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError:
    print('Encountered an attribute error')