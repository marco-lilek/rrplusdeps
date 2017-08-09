from pprint import pprint
from scipy.spatial import distance
'''
ops = 'CTC'
num_lines = 4
outfile.write('{} {} {}\n'.format(num_lines, num_lines * nppl, 1))
for l in range(num_lines):
  outfile.write("{}\n".format(nppl))

outfile.write("1\n")
for l in range(num_lines - 1):
  outfile.write("0\n")


'''

point_moves = {
  (2,1): [(1,2), (3,2)],
  (1,2): [(0,3), (2,3)],
  (3,2): [(2,3), (4,3)],
  (2,3): [(1,4), (3,4)],
  (1,0): [None, (2,1)],
  (2,-1): [None, (3, 0)],
  (4,-1): [(3, 0), None],
  (3,0): [(2,1), (4,1)],
  (0,1): [None, (1,2)],
  (4,1): [(3,2), (5,2)], # sdf
  (5,0): [(4,1), None],
  (5,2): [(4,3), (6,3)],
  (4,3): [(3,4), (5,4)],
  (6,1): [(5,2), None],
  (3,4): [(2,5), (4,5)]
}

starting_points = [
  (1,0,2),  (0,1,2), (5,0,0), (6,1,0), (2,-1, 2), (4,-1, 0)
]

ending_points = [
  (0,3), (1,4), (2,5), (4,5), (5,4), (6,3)
]

def adjusted(pos):
  return -6 + 4 * pos

def pad(path_points, line_points, x, y, pos, ry, rz):
  path_points.append((x,y,pos))
  line_points.append((x * g + adjusted(pos) * 1,ry,rz))

def apply_cross(path_points, line_points, x, y, pos, ry):
  ry += s
  if pos == 1:
    pad(path_points, line_points, x, y, pos + 0.5, ry, 2)
    pos = 2
  elif pos == 2:
    pad(path_points, line_points, x, y, pos - 0.5, ry, -2)
    pos = 1
  else:
    pass
    pad(path_points, line_points, x, y, pos, ry, 0)

  ry += s
  pad(path_points, line_points, x, y, pos, ry, 0)
  return pos, ry

def apply_twist(path_points, line_points, x, y, pos, ry):
  ry += s
  if pos % 2 == 0:
    pad(path_points, line_points, x, y, pos + 0.5, ry, -2)
    pos += 1
  else:
    pad(path_points, line_points, x, y, pos - 0.5, ry, 2)
    pos -= 1
  ry += s
  pad(path_points, line_points, x, y, pos, ry, 0)
  return pos, ry

ops = 'CTCT'
block_size = (len(ops)) * 4

g = block_size * 1.5
s = 2
start_scale = 4
pathpoints = []
linepaths = []
for sx,sy,sp in starting_points:
  for pos in (sp, sp+1):
    x,y = sx, sy

    rx, ry = float(x * g), float(y * g)
    z = 0
    path_points = []
    line_points = []
    while True:
      ryp = ry
      if (x,y) in ending_points:
        pad(path_points, line_points, x, y, pos, ry, 0)
        break

      if (x,y, sp) not in starting_points:
        ry -= block_size / 2
        pad(path_points, line_points, x, y, pos, ry, 0)
        for op in ops:
          if op == 'T':
            pos,ry = apply_twist(path_points, line_points, x, y, pos, ry)
          elif op == 'C':
            pos,ry = apply_cross(path_points, line_points, x, y, pos, ry)
          else:
            pass
        ry = ry - block_size / 2
      else:
        pad(path_points, line_points, x, y, pos, ry, 0)
      
      #print (x,y,pos)
      ry = ry + g

      leftmove = point_moves[(x,y)][0]
      rightmove = point_moves[(x,y)][1]
      if pos <= 1:
        x,y = leftmove
        pos += 2
      else:
        x,y = rightmove
        pos -= 2

    pathpoints.append(path_points)
    linepaths.append(line_points)

print 'done'

#pprint(pathpoints)

def midpoint(p1, p2):
  x1,y1,z1 = p1
  x2,y2,z2 = p2
  return (float(x1 + x2) / 2, float(y1 + y2) / 2,float(z1 + z2) / 2)

# adjustment
def smooth(linepaths, safe=False):
  for linepath in linepaths:
    i = 0
    while i < len(linepath) - 1:
      if not safe or distance.euclidean(linepath[i],linepath[i+1]) > 5:
        linepath.insert(i+1, midpoint(linepath[i], linepath[i+1]))
      i += 2

smooth(linepaths)
smooth(linepaths, True)
pprint(linepaths)

#print sum([len(x) for x in linepaths])

# file ops
outfile = open("test.vect", "w")
outfile.write("VECT\n")

# header info
outfile.write('{} {} {}\n'.format(
  len(linepaths), 
  sum([len(x) for x in linepaths]),
  1))

for lpath in linepaths:
  outfile.write("{}\n".format(len(lpath)))

outfile.write("1\n")
for l in range(len(linepaths) - 1):
  outfile.write("0\n")

# data entries
xyz_format = '{} {} {}\n'
for linepath in linepaths:
  for x,y,z in linepath:
    outfile.write(xyz_format.format(round(x,6),round(y,6),round(z,6)))

outfile.write('{} {} {} {}'.format(1,1,1,1))
outfile.close()

'''
ls = 3 # space between lines
ll = 7 # line length
zmag = 2
for l in range(num_lines):
  x,y,z = 0, ls * l, 0
  outfile.write(xyz.format(x, y, z))

  for op in ops:
    if op == 'T':
      d = (-1) ** (y / ls)
      y = y + ls * 0.5 * d
      z -= d * zmag
    else:
      movez = 0
      move = 0
      if y == ls:
        move = 0.5 * ls
        movez = zmag
      elif y == ls * 2:
        move = -(0.5 * ls)
        movez = -zmag
      y += move
      z += movez
    x += ll
    outfile.write(xyz.format(x, y, z))

    if op == 'T':
      y = y + ls * 0.5 * d
      z += d
    else:
      y += move
      z -= movez
    x += ll
    outfile.write(xyz.format(x, y, z))


outfile.write('{} {} {} {}'.format(1,1,1,1))
outfile.close()
'''