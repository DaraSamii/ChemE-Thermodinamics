from math import exp
from scipy.optimize import root


def Antonie(T, param):
    '''
    T : Temprature in Centigrade
    param: list of [A, B, C]

    This Function Computes saturated Pressure at given Temprature
    '''
    if len(param) != 3:
        raise Exception("param must be list of [A, B, C]")

    A = param[0]
    B = param[1]
    C = param[2]

    ln_P_sat = A - B/(T + C)
    P_sat = exp(ln_P_sat)

    return P_sat
#
#
#
#
#
#
def BUBL_P(X, T, all_params):
    '''
    X = list of mollar fractions of liquid like [0.2 ,0.8] or [0.1 0.2 0.7]
        Sumation of X list must be 1.0

    T = Temprature in Centigrade

    all_params = list of parameters for Antonie equations

    example for all params:
        all_params = [[A1, B1, C1],
                      [A2, B2, C2],
                      [A3, B3, C3]] 
    '''
    
    # checking function inputs to be correct
    if len(X) != len(all_params):
        raise Exception(
            "count of ellements of X and all_params must be equal!")

    if sum(X) != 1.0:
        raise Exception("Summation of all X list must be 1.0")

    # Computing Pressure
    P = 0
    for i in range(len(X)):
        xi = X[i]
        param = all_params[i] #[A, B, C]

        # Computing Saturated Pressureof i_th Component at Tempreture given with Antonie Equation
        P_sat_i = Antonie(T, param)

        # Computing the General Pressure with Raoult's law
        P += P_sat_i * xi

    # Computing Y mollar fractions of vapor
    Y = []
    for i in range(len(X)):
        xi = X[i]
        param = all_params[i]

        # Computing Saturated Pressureof i_th Component at Tempreture given with Antonie Equation
        P_sat_i = Antonie(T, param)

        # Computing i_th vapor mollar fraction with raoult's law
        yi = (xi * P_sat_i)/P
        Y.append(yi)

    # returning results
    return P, Y
#
#
#
#
#
#
def BUBL_T(X, P, all_params):
    """
    X = list of mollar fractions of liquid like [0.2 ,0.8] or [0.1 0.2 0.7]
        Sumation of X list must be 1.0

    P = Pressure in kPa

    all_params = list of parameters for Antonie equations

    example for all params:
        all_params = [[A1, B1, C1],
                      [A2, B2, C2],
                      [A3, B3, C3]]
    """
    # creating root finding function
    def func(T):
        return (P - BUBL_P(X, T, all_params)[0])

    # solving and finding Temprature
    solve = root(func, 20, method='lm')
    T = solve['x'][0]

    # computing Y mollar fra
    Y = BUBL_P(X, T, all_params)[1]

    # Computing Y mollar fractions of vapor
    return T, Y
#
#
#
#
#
#
def DEW_P(Y, T, all_params):
    '''
    Y = list of mollar fractions of vapor like [0.2 ,0.8] or [0.1 0.2 0.7]
        Sumation of X list must be 1.0

    T = Temprature in Centigrade

    all_params = list of parameters for Antonie equations

    example for all params:
        all_params = [[A1, B1, C1],
                      [A2, B2, C2],
                      [A3, B3, C3]] 
    '''
    # Checking function inputs are correct
    if len(Y) != len(all_params):
        raise Exception(
            "count of ellements of Y and all_params must be equal!")

    if sum(Y) != 1.0:
        raise Exception("Summation of all Y list must be 1.0")

    # Computing Pressure
    sum_yi_Psat = 0
    for i in range(len(Y)):
        yi = Y[i]
        param = all_params[i]

        # Computing Saturated Pressureof i_th Component at Tempreture given with Antonie Equation
        P_sat_i = Antonie(T, param)

        sum_yi_Psat += yi/P_sat_i

    # Computing the General Pressure with Raoult's law
    P = 1/sum_yi_Psat

    # Computing X mollar fractions of liquid
    X = []
    for i in range(len(Y)):
        yi = Y[i]
        param = all_params[i]

        # Computing Saturated Pressureof i_th Component at Tempreture given with Antonie Equation
        P_sat_i = Antonie(T, param)

        # Computing i_th liquid mollar fraction with raoult's law
        xi = (yi * P)/P_sat_i
        X.append(xi)

    # returning results
    return P, X
#
#
#
#
#
#
def DEW_T(Y, P, all_params):
    """
    Y = list of mollar fractions of vapor like [0.2 ,0.8] or [0.1 0.2 0.7]
        Sumation of X list must be 1.0

    P = Pressure in kPa

    all_params = list of parameters for Antonie equations

    example for all params:
        all_params = [[A1, B1, C1],
                      [A2, B2, C2],
                      [A3, B3, C3]]
    """
    # creating root finding function
    def func(T):
        return (P - DEW_P(Y, T, all_params)[0])

    # solving and finding Temprature
    solve = root(func, 20, method='lm')
    T = solve['x'][0]

    # Computing X mollar fractions of liqui
    X = DEW_P(Y, T, all_params)[1]

    return T, X