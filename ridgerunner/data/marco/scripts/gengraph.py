from util import *
from pprint import pprint
import numpy as np
import sys
import copy

OPS_YSHIFT = 2
UP = (0,1)

def get_basename(filename):
  doti = filename.rfind('.')
  return filename[:doti]

def read_infile(filename):
  infile = open(filename, 'r')
  _, num_starting = infile.readline().split(' ')
  num_starting = int(num_starting)
  starting = {}
  for i in range(num_starting):
    x,y,ax,ay,pside = [int(v) for v in infile.readline().split(' ')]
    starting[(x,y)] = ((ax, ay), pside)
  
  _, num_ending = infile.readline().split(' ')
  num_ending = int(num_ending)
  ending = []
  for i in range(num_ending):
    x,y = [int(v) for v in infile.readline().split(' ')]
    ending.append((x, y))

  typemap = {}
  maxtype = 0
  adj = {}
  _, num_adj = infile.readline().split(' ')
  num_adj = int(num_adj)
  for i in range(num_adj):
    x,y,ax,ay,bx,by, t = [int(v) for v in infile.readline().split(' ')]
    typemap[(x,y)] = t
    if t > maxtype:
      maxtype = t
    adj[(x,y)] = [(ax, ay), (bx, by)]

  rot = {}
  _, num_rot = infile.readline().split(' ')
  num_rot = int(num_rot)
  for i in range(num_rot):
    x,y,r = [float(v) for v in infile.readline().split(' ')]
    rot[(int(x), int(y))] = r

  incoming = {}
  _, num_incoming = infile.readline().split(' ')
  for i in range(int(num_incoming)):
    x,y,lx,ly,rx,ry = [int(v) for v in infile.readline().split(' ')]
    incoming[(x,y)] = [(lx,ly), (rx, ry)]

  infile.close()
  return starting, ending, adj, rot, incoming, maxtype + 1, typemap

def pos_shift(pos):
  return float(-6 + 4 *(pos))

def pad_path(path, pos, ry, rz):
  path.append((pos_shift(pos), ry, rz))

def apply_cross(ops_path, pos, ry):
  ry += OPS_YSHIFT
  if pos == 1:

    pad_path(ops_path, pos + 0.5, ry, 2)
    pos = 2
  elif pos == 2:
    pad_path(ops_path, pos - 0.5, ry, -2)
    pos = 1
  else:
    pass
    pad_path(ops_path, pos, ry, 0)
  
  ry += OPS_YSHIFT
  pad_path(ops_path, pos, ry, 0)
  return pos, ry

def apply_twist(ops_path, pos, ry):
  ry += OPS_YSHIFT
  if pos % 2 == 0:
    pad_path(ops_path, pos + 0.5, ry, -2)
    pos += 1
  else:
    pad_path(ops_path, pos - 0.5, ry, 2)
    pos -= 1
  ry += OPS_YSHIFT
  pad_path(ops_path, pos, ry, 0)
  return pos, ry

def add_pin(ops_path, pos, location, ry):
  ry += OPS_YSHIFT
  pad_path(ops_path, pos, ry, 0)
  pinpoints = [(-8 + 4 * (location - 1), ry, i) for i in (-20,0,20)]
  ry += OPS_YSHIFT
  pad_path(ops_path, pos, ry, 0)
  return ry, pinpoints

def run_ops(pos, ops, blocksize):
  # treat like centered at 0,0 and going upwards
  pins = []
  rx, ry = pos_shift(pos),-blocksize / 2
  ops_path = []
  
  if pos == 1 or pos == 2:
    pad_path(ops_path, pos, ry - 5, 0)
  pad_path(ops_path, pos, ry, 0)

  for op in ops:
    if op == 'T':
      pos,ry = apply_twist(ops_path, pos, ry)
    elif op == 'C':
      pos,ry = apply_cross(ops_path, pos, ry)
    elif int(op) in range(1,5+1):
      ry,pinpoints = add_pin(ops_path, pos, int(op), ry)
      pins.append(pinpoints) 
    else:
      assert 0, 'unknown OP ' + op

  if pos == 1 or pos == 2:
    pad_path(ops_path, pos, ry + 5, 0)
  return ops_path,pos, pins

linepaths = []
pins = {}

#ops = 'CT3CT'
startingtrans, ending, adj, rot, incoming, num_types, typemap = read_infile(sys.argv[1])

maxopsl = 0
assert len(sys.argv) > 2, 'need at least 1 pattern'

ops_paths = []
print num_types
for i in range(num_types):
  if len(sys.argv) == 3:
    pattern = sys.argv[2]
  else:
    pattern = sys.argv[2+i]

  newpath = [run_ops(x, pattern, len(pattern) * OPS_YSHIFT * 2) for x in range(4)]
  ops_paths.append(newpath)
  if len(pattern) > maxopsl:
    maxopsl = len(pattern)

gridsize = maxopsl * OPS_YSHIFT * 2 * 2.1

#startingtrans = {(1, 0): ((1, 1), 1)}
#startingtrans = {(3, 9): ((3, 10), 0)}
#print startingtrans
#assert 0
hitcount = {}
for sx,sy in startingtrans:
  pside = startingtrans[(sx, sy)][1]
  for pos in (0 + 2 * pside, 1 + 2* pside):
    x,y = sx, sy
    
    line_path = []
    dpath = []
    while True:
      gx, gy = x * gridsize, y * gridsize

      if (x,y) not in hitcount:
        hitcount[(x,y)] = -2
      hitcount[(x,y)] += 2
      #print (x,y), 'wow!'

      if (x,y) in ending and (x,y) != (sx, sy):
        rx, ry = rotate((0, 0), rot[(x,y)], (0,0))
        rz = hitcount[(x,y)]
        line_path.append((rx + gx, ry + gy, rz))
        dpath.append((x,y))
        break
      elif (x,y) in startingtrans:
        rx, ry = rotate((0, 0), rot[(x,y)], (0,0))
        rz = hitcount[(x,y)]
        line_path.append((rx + gx, ry + gy, rz))

        dpath.append((x,y))
        nx,ny = startingtrans[(x, y)][0]
        #print x, y, nx, ny, pos,
        pos = ((pos) % 2) + 2 * (not incoming[(nx,ny)].index((x,y)))
        #print pos
        x,y = nx, ny

        continue
      
      #pad_path(line_path, x * gridsize, y * gridsize, 0)
      #ops_path = [(lambda tz, tx, ty: (tx + gx, ty + gy, tz))(t[2], *rotate(t, rot[(x,y)])) for t in ops_paths[pos][0]]
      #pprint([(t[0] + gx, t[1] + gy, t[2]) for t in ops_paths[pos][0]])
      ops_path = [(lambda tx, ty, tz: (rotate((tx + gx, ty + gy), rot[(x,y)], (gx,gy))) + (tz,))(*t) for t in ops_paths[typemap[(x,y)]][pos][0]]

      line_path += ops_path
      pos = ops_paths[typemap[(x,y)]][pos][1]
      print 'sdf'
      pins[(x,y)] = []
      for pinpath in ops_paths[typemap[(x,y)]][pos][2]:
        pins[(x,y)].append([(lambda tx, ty, tz: (rotate((tx + gx, ty + gy), rot[(x,y)], (gx,gy))) + (tz,))(*t) for t in pinpath])
      dpath.append((x,y))

      lmove = adj[(x,y)][0]
      rmove = adj[(x,y)][1]
      if pos <= 1:
        nx,ny= rmove
      else:
        nx,ny = lmove

      #print x, y, nx, ny, pos,
      if (nx, ny) not in ending: # not perfect
        pos = (pos) % 2 + 2 * (not incoming[(nx,ny)].index((x,y)))

      #print pos
      x,y = nx, ny

    #print '____'
    #pprint(line_path)
    #print dpath
    linepaths.append(line_path)

print 'done'

flattened_pins = []
for t in pins:
  for e in pins[t]:
    flattened_pins.append(e)

#pins = [pins[k] for k in pins]

outfile = open(get_basename(sys.argv[1]) + '.vect', "w")
outfile.write("VECT\n")

# header info
outfile.write('{} {} {}\n'.format(
  len(linepaths) + len(flattened_pins), 
  sum([len(x) for x in linepaths]) + sum([len(x) for x in flattened_pins]),
  len(linepaths) + (len(flattened_pins) if len(flattened_pins) > 0 else 0)))

for lpath in linepaths:
  outfile.write("{}\n".format(len(lpath)))

for pinpath in flattened_pins:
  outfile.write("{}\n".format(len(pinpath)))

#pprint(pins)

for l in range(len(linepaths)):
  outfile.write("1\n")

if len(flattened_pins) > 0:
  for p in range(len(flattened_pins)):
    outfile.write("1\n")

# data entries
xyz_format = '{} {} {}\n'
for linepath in linepaths:
  for x,y,z in linepath:
    outfile.write(xyz_format.format(round(x,6),round(y,6),round(z,6)))

for pinpath in flattened_pins:
  for x,y,z in pinpath:
    outfile.write(xyz_format.format(round(x,6),round(y,6),round(z,6)))

for l in range(len(linepaths)):
  cintensity = float(l) /len(linepaths)
  outfile.write('{} {} {} {}\n'.format(cintensity,cintensity,1,1))

if len(flattened_pins) > 0:
  for p in range(len(flattened_pins)):
    outfile.write('0 0 0 0\n')

outfile.close()
