import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

from particle import ParticleBuffer, Particle
from box import Box3d, Grid


"""The Example Environment and Particles"""

box = Box3d(-2, -2, -2, 2, 2, 2, 0.8)
buffer = ParticleBuffer()      #position, acc, vel, density, pressure
for x in range(-3, 3):
  for y in range(4, 5):
    for z in range(-3, 3):
      buffer.addParticle([x/5,y/5,z/5], [0,0,0], [0,0,0], 0, 0)
# for x in range(-3, 1):
#   for y in range(2, 4):
#     for z in range(-3, 1):
#       buffer.addParticle([x/10,y/10,z/10], [0,0,0], [0,0,0], 0, 0)
#box, timestep, threshold, minIter, ext, miu, mass_rev, rou0
mass_rev = 50
print(buffer.particles[-1].position)
#data = open('position10.txt', 'w')
#data.write(str(len(buffer.particles))+'\n')
plt.ion()

ax = plt.subplot(1, 1, 1, projection='3d') # 创建一个三维的绘图工程
# 将数据点分成三部分画，在颜色上有区分度
ax.set_zlabel('Z') # 坐标轴
ax.set_ylabel('Y')
ax.set_xlabel('X')
plt.pause(1)

for i in range(200):
  plt.cla()
  buffer.update(box, 0.005, 0.1, 5, [0, -9.8/mass_rev, 0], 0.6, mass_rev, 1)
  print(buffer.particles[-1].position)
  x = []
  y = []
  z = []
  plt.xlim(-1, 1)
  plt.ylim(-1, 1)
  for particle in buffer.particles:
    #data.write(str(format(particle.position[0],'.4f')) + ' ' + str(format(particle.position[1],'.4f')) + ' ' + str(format(particle.position[2],'.4f')) + '\n')
    x.append(particle.position[0])
    y.append(particle.position[1])
    z.append(particle.position[2])
    ax.scatter([0,1], [0,1], [0,1], c='r')
    ax.scatter(x, z, y, c='b') # 绘制数据点

  plt.show()
  plt.pause(0.001)
# from google.colab import drive
# drive.mount('/content/drive')
