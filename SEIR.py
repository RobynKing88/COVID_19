from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
from openpyxl import Workbook, load_workbook
from datetime import datetime
from dateutil.relativedelta import relativedelta

def deriv(y, t, N, beta, gamma, alpha, mu, i):

    # N = total population
    # gamma = 1/ (sickness period)
    # alpha = 1/ (incubation period)
    # beta = (R0 * alpha) ie ( Replication Rate / Incubation period )
    # beta is probability of a successful transmission
    # mu = Case Fatality Rate
    # i = Rate of NEW infections entering the population, as a rate of total population

    S, E, I, R, D = y

    dSdt = -beta * S * I / (S+I+R)
    dEdt = beta * S * I / (S+I+R) - alpha*E
    dIdt = i*N + alpha*E - gamma * I
    dRdt = gamma * I * (1 - mu)
    dDdt = gamma * I * mu

    return dSdt, dEdt, dIdt, dRdt,dDdt

def get_dailynewcases(CaseHistory):

    # This function takes an array = [[date_array],[total_infections_array]]
    # and returns an array = [[date_array],[daily_new_infections_array]]

    # created for us to model the progression for existing past cases separately to new infections
    # because we can then assume existing cases are isolating, and only new cases have been spreading the virus
    # without knowing

    CaseDiff_arr = []
    CaseDiff = np.diff(CaseHistory[1])
    CaseDiff_arr.append(CaseHistory[1][0])
    for n in CaseDiff:
        CaseDiff_arr.append(n)
    CaseHistory.append(CaseDiff_arr)
    CaseHistory.pop(1)

    return CaseHistory

def SEIR(N,I,R,D,R0,Tinc,Tinf,CFR,Forecast_days,i):

    # This function integrates over the set of differential equations to give running total of
    # S = Susceptible population
    # E = Exposed population (incubating, not infectious or infected yet)
    # I = Infected population (infected and infectious)
    # R = Recovered and immune population
    # D = Deceased population

    # Inputs are as follows:
    # N = total population
    # I = initial Infected population
    # R = initial Recovered population
    # D = initial Deaths
    # R0 = Replication Rate
    # Tinc = Incubation period
    # Tinf = Infectious period
    # CFR = Case Fatality Rate
    # Forecast days = period of simulation
    # i = Rate of NEW infections entering the population, as a rate of total population

    # ** We haven't included the natural birth and mortality rates in the population
    # ** We can easily add but may just be adding noise

    alpha = 1/Tinc
    gamma = 1/Tinf
    beta = R0/Tinc
    mu = CFR
    E = beta*I
    S = N - I - R - D - E

    t = np.linspace(0, Forecast_days, Forecast_days)

    y0 = S, E, I, R, D
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, alpha, mu, i))
    S, E, I, R, D = ret.T

    return S, E, I, R, D

def SEIR_array(N,I,R,D,R0,Tinc,Tinf,CFR,Forecast_days,i):

    #Takes output and formats into array to be output easily

    S, E, I, R, D = SEIR(N,I,R,D,R0,Tinc,Tinf,CFR,Forecast_days,i)

    SEIRD = [[],[],[],[],[]]

    for sus in S:
        SEIRD[0].append(round(sus))
    for exp in E:
        SEIRD[1].append(round(exp))
    for inf in I:
        SEIRD[2].append(round(inf))
    for rec in R:
        SEIRD[3].append(round(rec))
    for dea in D:
        SEIRD[4].append(round(dea))

    return SEIRD

test = SEIR_array(59000000,709,0,0,2.3,4,14,0.02,14,0)
print(test)
