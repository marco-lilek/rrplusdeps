from util import *
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
import sys

TILETIMES = 4
UP = (0,1)

def read_infile(filename):
  infile = open(filename, 'r')
  _, pattern_h, pattern_w = infile.readline().rstrip().split('\t')

  pattern = []
  for line in infile:
    triplets = line.rstrip().split('\t')
    for t in triplets:
      pattern_t = []
      coords = [int(x) for x in t[1:-1].split(',')]    
      for x,y in zip(coords[0::2], coords[1::2]):
        pattern_t.append([x,y])
      pattern.append(pattern_t)

  infile.close()
  return pattern, int(pattern_w), int(pattern_h)

def shift(points, sx, sy, pattern_w, pattern_h):
  xd = sx * pattern_w
  yd = sy * pattern_h
  x,y = points[0]
  xa,ya = points[1]
  xb,yb = points[2]
  cur = (x + xd, y + yd)
  a = (xa + xd,ya + yd)
  b = (xb + xd,yb + yd)
  return cur, a, b

def add_line(plt, x1,y1,x2,y2, color, middle_arrow=True, zorder=0):
  plt.arrow(x1, y1, x2- x1, y2 - y1, shape='full', lw=0.5, color=color, length_includes_head=True,head_width=0 if middle_arrow else 0.1, zorder=zorder)
  if middle_arrow:
    plt.arrow(x1, y1, float(x2- x1)/2, float(y2 - y1)/2, shape='full', lw=0.3, color=color, length_includes_head=True,head_width=0.1, zorder=zorder)

def get_basename(filename):
  doti = filename.rfind('.')
  return filename[:doti]

pattern, pattern_w, pattern_h = read_infile(sys.argv[1]) 

adj = {}        # (x,y) -> [out left, out right] TO not from
rot = {}        # (x,y) -> rotation from UP
incoming = {}   # incoming into each vertex
typemap = {}    # (x,y) -> type identifier

for tilex in range(TILETIMES):
  for tiley in range(TILETIMES):
    for i, points in enumerate(pattern):
      cur,a,b = shift(points, tilex, tiley, pattern_w, pattern_h)
      typemap[cur] = i
      # for now just guess which incoming line is left and right, we can correct next loop
      adj[cur] = [a, b] 
      for t in (a, b):
        if t not in incoming:
          incoming[t] = []
        incoming[t].append(cur)

corrected = set()   # some points will not have corrected incoming
for t in adj:
  if t in incoming and len(incoming[t]) == 2 and t in adj:
    assert len(incoming[t]) <= 2, '' + str(t) + str(incoming[t])
    x,y = t
    
    # outgoing 2
    a = adj[t][0]
    b = adj[t][1]

    # incoming 2
    c = incoming[t][0]
    d = incoming[t][1]

    # vector to outging 2
    out1 = (a[0] - x, a[1] - y)
    out2 = (b[0] - x, b[1] - y)

    # vector to incoming 2
    in1 = (c[0] - x, c[1] - y)
    in2 = (d[0] - x, d[1] - y)
    in1out1a = angle(in1, out1)
    in1out2a = angle(in1, out2)

    # clockwise angle between input line and output lines
    if in1out1a < in1out2a:
      adj[t] = [a,b]
      outl = out1
    else:
      adj[t] = [b,a]
      outl = out2

    # now that we know the left output, we want to figure out which is the left input
    # whichever one has the smaller angle to the left input
    in1outla = angle(in1, outl)
    in2outla = angle(in2, outl)
    if in1outla > in2outla:
      incoming[t] = [d,c]
    else:
      incoming[t] = [c,d]

    # now we need to figure out the angle the intersection would be going
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

    corrected.add(t) # this point has now a corrected left and right

# remove any points that haven't been corrected
for x in range(TILETIMES * pattern_w):
  for y in range(TILETIMES * pattern_h):
    cur = x,y
    if cur not in corrected:
      #print 'del', cur
      if cur in incoming:
        del incoming[cur]
      if cur in adj:
        del adj[cur]

bbox = [[0,TILETIMES * pattern_w],[0,TILETIMES * pattern_h]]

if len(sys.argv) == 6:
  bbox = [[int(sys.argv[2]), int(sys.argv[3])], [int(sys.argv[4]), int(sys.argv[5])]]
  print 'setting bbox to ', bbox

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
  add_line(plt, xl,yl,x,y,'g')
  add_line(plt, xr,yr,x,y,'turquoise')

  #print x,y
  rotx, roty = rotate((x,y+1), rot[a], a)
  #print rotx, roty
  add_line(plt, x,y,rotx,roty,'black', False)

  xl,yl = adj[a][0]
  xr,yr = adj[a][1]

  add_line(plt, x,y,xl,yl,'r', zorder=10)
  add_line(plt, x,y,xr,yr,'b', zorder=10)

plt.scatter(xs, ys)

# cleanup typemap for active points (number range only samples from these points)
typemapmap = {}
numtypes = 0
for pt in typemap:
  if pt in activepts:
    if typemap[pt] not in typemapmap:
      typemapmap[typemap[pt]] = numtypes
      numtypes += 1
    typemap[pt] = typemapmap[typemap[pt]]

for pt in typemap:
  if pt in activepts:    
    plt.text(pt[0], pt[1], typemap[pt], color='#ff00ff', weight='heavy') 
plt.xticks(range(TILETIMES * pattern_w))
plt.yticks(range(TILETIMES * pattern_h))

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
    # getting the angle on the output points
    x1,y1 = inps[0]
    x2,y2 = inps[1]
    dir1 = (x1 - vx, y1 -vy)
    dir2 = (x2- vx, y2 -vy) 
    angle1 = angle(UP, dir1)
    angle2 = angle(UP, dir2)
    midangl = interior(norm(angle(dir1, dir2))) / 2
    if angle1 <= angle2:
      rot[v] = angle1 + midangl + np.pi
    else:
      rot[v] = angle2 + midangl + np.pi

  rotx, roty = rotate((vx,vy+1), rot[v], v)
  plt.arrow(vx, vy, rotx - vx, roty - vy, shape='full', lw=0.6, color='black', length_includes_head=True, head_width=.05)
  plt.scatter([vx], [vy], color='purple')

plt.gca().invert_yaxis()

base_filename = get_basename(sys.argv[1])

plt.savefig(base_filename + '.pattern.png')
#print starting_pts
#print ending_pts

###########################################################
# writing all the data

outfile = open(base_filename + '.pattern', 'w')
outfile.write('STARTING {}\n'.format(len(starting_pts)))
for e in starting_pts:
  fout = list(sum([e[0], e[1]], ()))
  outfile.write('{} {} {} {} {r}\n'.format(*fout, r=e[2]))

outfile.write('ENDING {}\n'.format(len(ending_pts)))
for pt in ending_pts:
  outfile.write('{} {}\n'.format(*pt))

outfile.write('ACTIVE {}\n'.format(len(activepts)))
nexttype = 0
typemapmap = {}
for pt in activepts:

  # just numbering it appropriately
  fout = list(sum([pt, adj[pt][0], adj[pt][1], (typemap[pt],)], ()))
  outfile.write('{} {} {} {} {} {} {}\n'.format(*fout))

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
