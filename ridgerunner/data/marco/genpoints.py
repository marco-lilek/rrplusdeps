from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
import sys

def dot(v1, v2):
  return v1[0] * v2[0] + v1[1]*v2[1]

def det(v1, v2):
  return v1[0] * v2[1] - v1[1]*v2[0]

def rotate(v, angle, anchor):
  x,y = v
  x -= anchor[0]
  y -= anchor[1]
  cos_t = np.cos(angle)
  sin_t = np.sin(angle)

  nx = x * cos_t - y * sin_t
  ny = x * sin_t + y * cos_t
  #print (nx + anchor[0], ny + anchor[1])
  return (nx + anchor[0], ny + anchor[1])

def norm(a):
  return np.mod(a, np.pi * 2)

def interior(a):
  a = norm(a)
  if a > np.pi:
    return np.pi * 2 - a
  return a

def angle(v1, v2):
  return norm(np.arctan2(-det(v1, v2), -dot(v1, v2)) + 3 * np.pi)

UP = (0,1)
adj = {} # (x,y) -> [ad, adj] TO not from
rot = {} # (x,y) -> rotation from UP
incoming = {} #incoming into each vertex

def read_infile(filename):
  infile = open(filename, 'r')
  _, pattern_h, pattern_w = infile.readline().rstrip().split('\t')

  pattern = []
  for line in infile:
    triplets = line.rstrip().split('\t')
    for t in triplets:
      pattern_t = []
      print t
      coords = [int(x) for x in t[1:-1].split(',')]    
      for x,y in zip(coords[0::2], coords[1::2]):
        pattern_t.append([x,y])
      pattern.append(pattern_t)

  infile.close()
  return pattern, int(pattern_w), int(pattern_h)

pattern, pattern_w, pattern_h = read_infile(sys.argv[1]) 

#pattern = [[[0,0],[1,0],[0,1]],[[0,1],[-1,1],[0,2]]]
print pattern, pattern_h, pattern_h

def shift(points, sx, sy):
  global pattern_w, pattern_h
  xd = sx * pattern_w
  yd = sy * pattern_h
  x,y = points[0]
  xa,ya = points[1]
  xb,yb = points[2]
  cur = (x + xd, y + yd)
  a = (xa + xd,ya + yd)
  b = (xb + xd,yb + yd)
  return cur, a, b

tilew = 4
for tilex in range(tilew):
  for tiley in range(tilew):
    for points in pattern:
      cur,a,b = shift(points, tilex, tiley)

      # need to bootstrap correct lr adj for starters
      adj[cur] = [a, b]
      for t in (a, b):
        if t not in incoming:
          incoming[t] = []
        incoming[t].append(cur)



#pprint(adj)
corrected = set()
for tilex in range(tilew):
  for tiley in range(tilew):
    for points in pattern:
      t,_,_ = shift(points, tilex, tiley)
      if t in incoming and len(incoming[t]) == 2 and t in adj:
        assert len(incoming[t]) <= 2, '' + str(t) + str(incoming[t])
        corrected.add(t)
        #print 'ok',t
        #print incoming[t]
        #print adj[t]
        a = adj[t][0]
        b = adj[t][1]
        c = incoming[t][0]
        d = incoming[t][1]
        x,y = t
        out1 = (a[0] - x, a[1] - y)
        out2 = (b[0] - x, b[1] - y)
        in1 = (c[0] - x, c[1] - y)
        in2 = (d[0] - x, d[1] - y)
        in1out1a = angle(in1, out1)
        in1out2a = angle(in1, out2)
        #print in1out1a, in1out2a
        if in1out1a < in1out2a:
          adj[t] = [a,b]
          outl = out1
        else:
          adj[t] = [b,a]
          outl = out2

        in1outla = angle(in1, outl)
        in2outla = angle(in2, outl)
        if in1outla > in2outla:
          incoming[t] = [d,c]
        else:
          incoming[t] = [c,d]


        li = incoming[t][0][0]-x,incoming[t][0][1]-y
        ri = incoming[t][1][0]-x,incoming[t][1][1]-y
        lo = adj[t][0][0]-x,adj[t][0][1]-y
        ro = adj[t][1][0]-x,adj[t][1][1]-y
        midin = interior(angle(li, ri)) / 2 + angle(UP, ri) #(angle(UP, li) + angle(UP, ri))
        #midin = interior(midin) / 2
        midout = interior(angle(lo, ro)) / 2 + angle(UP, lo) #(angle(UP, lo) + angle(UP, ro))
        #midout = interior(midout) / 2
        if midout > midin:
          btw = (midout - midin)/ 2 + midin + np.pi / 2
        else:
          btw = (midout - midin)/ 2 + midin - np.pi / 2
        #print  midin, midout, btw, angle(UP, lo)
        rot[t] = btw #norm(btw - np.pi / 2) # + angle(UP, lo) #btw + angle(UP,lo) + 0.1

for x in range(tilew * pattern_w):
  for y in range(tilew * pattern_h):
    cur = x,y
    if cur not in corrected:
      #print 'del', cur
      if cur in incoming:
        del incoming[cur]
      if cur in adj:
        del adj[cur]

#pprint(adj)

bbox = [[0,tilew * pattern_w],[0,tilew * pattern_h]]
bbox = [[0,17],[0,17]]
#bbox = [[1,1],[1,1]]
pts = adj.keys()
activepts = []
xs = []
ys = []

for a in pts:
  x,y = a
  if x < bbox[0][0] or x > bbox[0][1] or y < bbox[1][0] or y > bbox[1][1]:
    continue  

  activepts.append(a)
  xs.append(x)
  ys.append(y)
  xl,yl = incoming[a][0]
  xr,yr = incoming[a][1]
  plt.arrow(xl, yl, x- xl, y - yl, shape='full', lw=0.3, color='g', length_includes_head=True,head_width=.05)
  plt.arrow(xr, yr, x- xr, y - yr, shape='full', lw=0.3, color='turquoise', length_includes_head=True, head_width=.05)

  #print x,y
  rotx, roty = rotate((x,y+1), rot[a], a)
  #print rotx, roty
  plt.arrow(x, y, rotx - x, roty - y, shape='full', lw=0.6, color='black', length_includes_head=True, head_width=.05)

  xl,yl = adj[a][0]
  xr,yr = adj[a][1]
  plt.arrow(x, y, xl - x, yl - y, shape='full', lw=0.5, color='r', length_includes_head=True,head_width=.05, zorder=10)
  plt.arrow(x, y, xr - x, yr - y, shape='full', lw=0.5, color='b', length_includes_head=True, head_width=.05, zorder=10)

plt.scatter(xs, ys)
plt.xticks(range(tilew * pattern_w))
plt.yticks(range(tilew * pattern_h))

###########################################################

# point_moves = adj

starting_pts = []
ending_pts = {} # map of node to inputs []
for pt in activepts:
  #print pt
  if incoming[pt][0] not in activepts or incoming[pt][1] not in activepts:
    #print incoming[pt]
    assert pt in incoming and len(incoming[pt]) == 2
    x,y = pt
    a = incoming[pt][0]
    b = incoming[pt][1]
    for v in (a,b):
      if v not in activepts:
        starting_pts.append((v, pt, int(v==a)))
        vx,vy = v
        rot[v] = angle(UP, (x-vx,y-vy))
        plt.scatter([vx], [vy], color='orange')

for pt in activepts:
  x,y = pt
  lout = adj[pt][0]
  rout = adj[pt][1]
  for v in (lout, rout):
    if v not in activepts:
      if v not in ending_pts:
        ending_pts[v] = []
      ending_pts[v].append(pt)

for v in ending_pts:
  inps = ending_pts[v]
  vx,vy = v
  assert len(inps) > 0 and len(inps) <= 2
  if len(inps) == 1:
    rot[v] = angle(UP, (vx-inps[0][0], vy-inps[0][1]))
  else:
    x1,y1 = inps[0]
    x2,y2 = inps[1]
    dir1 = (vx - x1, vy - y1)
    dir2 = (vx - x2, vy - y2) 
    angle1 = angle(UP, dir1)
    angle2 = angle(UP, dir1)
    if angle1 <= angle2:
      rot[v] = angle1 + angle(dir1, dir2) / 2
    else:
      rot[v] = angle2 + angle(dir1, dir2) / 2

  rotx, roty = rotate((vx,vy+1), rot[v], v)
  plt.arrow(vx, vy, rotx - vx, roty - vy, shape='full', lw=0.6, color='black', length_includes_head=True, head_width=.05)
  plt.scatter([vx], [vy], color='purple')

plt.gca().invert_yaxis()

def get_basename(filename):
  doti = filename.rfind('.')
  return filename[:doti]

base_filename = get_basename(sys.argv[1])

plt.savefig(base_filename + '.pattern.png')
#print starting_pts
#print ending_pts

outfile = open(base_filename + '.pattern', 'w')
outfile.write('STARTING {}\n'.format(len(starting_pts)))
for e in starting_pts:
  fout = list(sum([e[0], e[1]], ()))
  outfile.write('{} {} {} {} {r}\n'.format(*fout, r=e[2]))

outfile.write('ENDING {}\n'.format(len(ending_pts)))
for pt in ending_pts:
  outfile.write('{} {}\n'.format(*pt))

outfile.write('ADJ {}\n'.format(len(activepts)))
for pt in activepts:
  fout = list(sum([pt, adj[pt][0], adj[pt][1]], ()))
  outfile.write('{} {} {} {} {} {}\n'.format(*fout))

outfile.write('ROT {}\n'.format(len(activepts) + len(starting_pts) + len(ending_pts)))
for pt in activepts:
  outfile.write('{} {} {}\n'.format(pt[0], pt[1], rot[pt]))
for pt in starting_pts:
  outfile.write('{} {} {}\n'.format(pt[0][0], pt[0][1], rot[pt[0]]))
for pt in ending_pts:
  outfile.write('{} {} {}\n'.format(pt[0], pt[1], rot[pt]))

outfile.write('INCOMING {}\n'.format(len(activepts)))
for pt in activepts:
  outfile.write('{} {} {} {} {} {}\n'.format(pt[0], pt[1],
    incoming[pt][0][0], incoming[pt][0][1], incoming[pt][1][0], incoming[pt][1][1]))

outfile.close()
