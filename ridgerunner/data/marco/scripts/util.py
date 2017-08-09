import numpy as np

def rotate(v, angle, anchor=(0,0)):
  x,y = v
  x -= anchor[0]
  y -= anchor[1]
  cos_t = np.cos(angle)
  sin_t = np.sin(angle)

  nx = x * cos_t - y * sin_t
  ny = x * sin_t + y * cos_t
  return (nx + anchor[0], ny + anchor[1])

def angle(v1, v2):
  return norm(np.arctan2(-det(v1, v2), -dot(v1, v2)) + 3 * np.pi)

def dot(v1, v2):
  return v1[0] * v2[0] + v1[1]*v2[1]

def det(v1, v2):
  return v1[0] * v2[1] - v1[1]*v2[0]

def norm(a):
  return np.mod(a, np.pi * 2)

def interior(a):
  a = norm(a)
  if a > np.pi:
    return np.pi * 2 - a
  return a