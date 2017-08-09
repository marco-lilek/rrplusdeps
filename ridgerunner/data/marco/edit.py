# python
import fileinput

xp = 0
yp = -5
zp = 0

count = 0
for line in fileinput.input():
  x,y,z = line.split(' ')
  x = float(x)
  y = float(y)
  z = float(z)
  print x + xp, y + yp, z + zp
  count += 1

print count
