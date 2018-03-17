import numpy

# functions for sampling n_out for no land cases to 
# correct for bias in model failure in these cases
# sampling skewed toward higher n_out values ensures
# final distribution for n_out is appproximately uniform

def no_land_no_methane():
    n_out=numpy.random.uniform(0.0,0.73)
    if numpy.random.uniform(0.0,1.0) < 0.45: # for no land case
        if n_out<0.66:
            n_out=numpy.random.uniform(0.0,0.73)
        if n_out<0.6:
            n_out=numpy.random.uniform(0.0,0.73)
        if n_out <0.45:
            n_out=numpy.random.uniform(0.0,0.73)
    return n_out

def no_land_with_methane():
    n_out=numpy.random.uniform(0.0,0.73)
    if numpy.random.uniform(0.0,1.0) < 0.6: # for no land and CH4 case
        if n_out<0.65:
            n_out=numpy.random.uniform(0.0,0.73)
        if n_out <0.45:
            n_out=numpy.random.uniform(0.0,0.73)
        if n_out <0.30:
            n_out=numpy.random.uniform(0.0,0.73)
    return n_out

