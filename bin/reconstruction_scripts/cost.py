import math

#COST_EPSILON = 1e-10
#COST_EPSILON = 10
COST_EPSILON = 1
#COST_EPSILON = 0.1

#MIN_PR = 1e-17
MIN_PR = 1e-6
#MIN_PR = 1e-10
MAX_PR = 1.0

def rcost(x):
    if x < MIN_PR:
        x = MIN_PR
    return -math.log(x) / math.log(2) + COST_EPSILON

def cost_to_p(c):
    return 2**(-c)

def min_cost():
    return rcost(MAX_PR)

def max_cost():
    return rcost(MIN_PR)
