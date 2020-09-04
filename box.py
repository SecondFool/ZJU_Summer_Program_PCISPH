import numpy as np
import math
from particle import ParticleBuffer, Particle

class Grid:
  def __init__(self, gIdx):
    self.gIdx = gIdx
    self.particles = []
  def removeParticle(self, idx):
    if idx not in self.particles:
      print('Cannot find', idx,'th particle in', self.gIdx, 'th Grid!')
      return
    self.particles.remove(idx)
  def addParticle(self, idx):
    if idx not in self.particles:
      self.particles.append(idx)
print('Build Class Grid Successfully!')

class Box3d:
  def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax, unitlength):
    self.xmin, self.xmax = xmin, xmax
    self.ymin, self.ymax = ymin, ymax
    self.zmin, self.zmax = zmin, zmax
    self.unitlength = unitlength
    self.xnum = math.ceil((xmax - xmin)/unitlength)
    self.ynum = math.ceil((ymax - ymin)/unitlength)
    self.znum = math.ceil((zmax - zmin)/unitlength)       #the number of small grids on each axis
    self.grids = [-1]*(self.xnum*self.ynum*self.znum)

  def getgIdx(self, position):
    return (math.floor((position[0]-self.xmin)/self.unitlength) + \
                      math.floor((position[1]-self.ymin)/self.unitlength)*self.xnum + \
                      math.floor((position[2]-self.zmin)/self.unitlength)*self.xnum*self.ynum)
    # if ((position[0]<self.xmin) or (position[0]>=self.xmax) or
    #     (position[1]<self.ymin) or (position[1]>=self.ymax) or
    #     (position[2]<self.zmin) or (position[2]>=self.zmax)):
    #     print('Warning: position is out of range!(That would be fine if it works, just ignore it)')
    
  def findNeighbourGrids(self, position):               #it would be better to set limitations not too close to the bound
    if ((position[0]<self.xmin) or (position[0]>=self.xmax) or
        (position[1]<self.ymin) or (position[1]>=self.ymax) or
        (position[2]<self.zmin) or (position[2]>=self.zmax)):
        #print('position is out of range!')
        return []
    LDFgIdx = self.getgIdx([position[0]-0.5*self.unitlength,position[1]-0.5*self.unitlength,position[2]-0.5*self.unitlength])     #LDF stands for the left-down-front grid in the 8 neighbouring grids
    neighbours = [LDFgIdx, LDFgIdx+1, 
                  LDFgIdx+self.xnum, LDFgIdx+self.xnum+1, 
                  LDFgIdx+self.xnum*self.ynum, LDFgIdx+self.xnum*self.ynum+1, 
                  LDFgIdx+self.xnum*(self.ynum+1), LDFgIdx+self.xnum*(self.ynum+1)+1]
    if (position[0] - 0.5*self.unitlength < self.xmin) :      #these are for boundary conditions
      neighbours.remove(LDFgIdx)                  #if we use some proper size of container, then we can comment all these
      neighbours.remove(LDFgIdx+self.xnum)
      neighbours.remove(LDFgIdx+self.xnum*self.ynum)
      neighbours.remove(LDFgIdx+self.xnum*(self.ynum+1))
    if (position[1] - 0.5*self.unitlength < self.ymin):
      if LDFgIdx in neighbours:
        neighbours.remove(LDFgIdx)
      neighbours.remove(LDFgIdx+1)
      if LDFgIdx+self.xnum*self.ynum in neighbours:
        neighbours.remove(LDFgIdx+self.xnum*self.ynum)
      neighbours.remove(LDFgIdx+self.xnum*self.ynum+1)
    if (position[2] - 0.5*self.unitlength < self.zmin):
      if LDFgIdx in neighbours:
        neighbours.remove(LDFgIdx)
      if LDFgIdx+1 in neighbours:
        neighbours.remove(LDFgIdx+1)
      if LDFgIdx+self.xnum in neighbours:
        neighbours.remove(LDFgIdx+self.xnum)
      neighbours.remove(LDFgIdx+self.xnum+1)
    return neighbours

  def loadParticles(self, particleBuffer):
    self.grids = [-1]*(self.xnum*self.ynum*self.znum)
    for particle in particleBuffer.particles:
      gIdx = self.getgIdx(particle.position)
      if gIdx >= len(self.grids) or gIdx < 0:
        print(particle.position)
        continue
      if self.grids[gIdx] == -1:
        self.grids[gIdx] = Grid(gIdx)
      self.grids[gIdx].addParticle(particle.idx)
print('Build Class Box3d Successfully!')

#test 2-6
box = Box3d(0, 0, 0, 4, 4, 4, 1)
buffer = ParticleBuffer()
buffer.addParticle([0,0,0], [1,1,1], [1,1,1], 9, 5)
buffer.addParticle([2,2,2], [3,4,5], [6,7,8], 11, 8)
buffer.addParticle([2.1,2.1,2.1], [3,4,5], [6,7,8], 11, 8)
buffer.addParticle([1.2,2.1,2.1], [3,4,5], [6,7,8], 11, 8)
box.loadParticles(buffer)
if (box.grids[0] != -1) and (box.grids[42] != -1) and (box.grids[41] != -1):
  print('Pass Test2')
if (box.getgIdx([0, 0, 0]) == 0) and (box.getgIdx([2, 2, 2]) == 42) and (box.getgIdx([1.2, 2.1, 2.1]) == 41):
  print('Pass Test3')
if (box.findNeighbourGrids([0, 0, 0]) == [0]):
  print('Pass Test4')
if (box.findNeighbourGrids([2, 2, 2]).sort() == [21, 22, 37, 38, 25, 26, 41, 42].sort()):
  print('Pass Test5')
if (box.grids[box.getgIdx([2, 2, 2])].particles.sort() == [1, 2].sort()):
  print('Pass Test6')
