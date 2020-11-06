#!/usr/bin/env python

import numpy as np
import random

from parameters import *

def spiketrain(rate, tmax):
  t = 0
  train = []
  while t<=tmax:
    t+= random.expovariate(rate)
    if t<=tmax:
      train.append(t)
  return train

def spiketrainBurst(rate1, rate2, tBurst, duration, tmax):
  if tBurst > tmax:
     print 'Warning: does not burst during simulation period'
  if rate1 > rate2:
     print 'Error: burst rate must be greater than tonic rate'
  t = tmax
  train = []
  while t > tBurst+duration:
    t-= random.expovariate(rate1)
    train.append(t)
  while t > tBurst:
    t-= random.expovariate(rate2)
    train.append(t)
  while t > 0:
    t-= random.expovariate(rate1)
    train.append(t)
  return train


