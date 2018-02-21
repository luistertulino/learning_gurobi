from gurobipy import *

num_voxels = 25 # range = 0-24
tumor_voxels = [7, 8, 12, 13, 17, 21, 22] # Vt
organs_risk_voxels = [6, 10, 11, 15, 16, 20] # Vor
normal_voxels = [0, 1, 2, 3, 4, 5, 9, 14, 18, 19, 23, 24] # Vs

angle_index = { 0 : 0, 90 : 1, 180 : 2, 270 : 3 }

selected_angles = [0, 180]

num_beamlets = 6

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

range_0_24 = range(25)

Q_0   = [ range_0_24, range_0_24, range_0_24, range_0_24, range_0_24 ]
Q_90  = Q_0
Q_180 = Q_0
Q_270 = Q_0


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



except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError:
    print('Encountered an attribute error')