import numpy as np
import math

class Particle:
  def __init__(self, position, acc, vel, density, pressure, idx):
    self.position = position                        #all vectors are three dimensions ordered as [x,y,z] 
    self.acc = acc
    self.vel = vel
    self.density = density                        
    self.pressure = pressure
    self.idx = idx
    self.gIdx = -1
    self.Ppos, self.Pvel, self.Pden, self.Pdenv, self.Ppre = self.position, self.vel, 0, 0, 0
    self.forcep = [0,0,0]
    self.forcee = [0,0,0]
    self.forcev = [0,0,0]
    self.neighp = []
  
  def updateDensity(self, mass, rlist, h):
    self.Pden = 0.0
    for r in rlist:
      self.Pden += mass*(self.SKFunction(r, h))
    return

  def updatePressure(self, mass, timestep, rou0, rlist, h):
    self.Ppre = 0.0
    beta = 2*(timestep*mass/rou0)**2
    A = 0
    B = 0
    p = 0
    for r in rlist:
      A += self.derSKFunction(r, h)
      B += (self.derSKFunction(r, h))**2
    if (-A**2-B != 0):
      self.Ppre += (-self.Pdenv)/(beta*(-A**2-B))
    return

  def updateForcep(self, mass, plist, rlist, dirlist):          #pressure list, radius list, direction list
    self.forcep = [0.0, 0.0, 0.0]
    # if (len(plist) != len(rlist)) or (len(plist) != len(dirlist)):
    #   print('list lenth unmatching error!')
    #   return
    for i in range(len(plist)):
      # self.forcep[0] += dirlist[i][0]*plist[i]/(rlist[i]**3)
      # self.forcep[1] += dirlist[i][1]*plist[i]/(rlist[i]**3)
      # self.forcep[2] += dirlist[i][2]*plist[i]/(rlist[i]**3)
      self.forcep[0] += dirlist[i][0]*plist[i]/(rlist[i]**2)
      self.forcep[1] += dirlist[i][1]*plist[i]/(rlist[i]**2)
      self.forcep[2] += dirlist[i][2]*plist[i]/(rlist[i]**2)
    return
  
  def updateForcef(self, box, mass):
    shift = [1, box.xnum, box.xnum*box.ynum, -1, -box.xnum, -box.xnum*box.ynum]
    gidxcurr = box.getgIdx(self.position)
    shiftidx = 0
    if abs(self.vel[1]) > abs(self.vel[0]):
      shiftidx += 1
      if abs(self.vel[2]) > abs(self.vel[1]):
        shiftidx += 1
    else:
      if abs(self.vel[2]) > abs(self.vel[0]):
        shiftidx += 2
    if max([abs(self.vel[0]), abs(self.vel[1]), abs(self.vel[2])]) < 0:
      shiftidx += 3
    if (gidxcurr + shift[shiftidx] >= 0) and (gidxcurr + shift[shiftidx] <= len(box.grids)):
      if box.grids[gidxcurr + shift[shiftidx]] == -1:
        self.forcee[0] -= self.vel[0]*0.0005
        self.forcee[1] -= self.vel[1]*0.0005
        self.forcee[2] -= self.vel[2]*0.0005
    return 0
      
  def SKFunction(self, r, h):
    if (0 < r) and (r < h):
      return (315/(63*3.1416*(h**9)))*r*((h**2-r**2)**3)
    else:
      return 0

  def derSKFunction(self, r, h):
    if (0 < r) and (r < h):
      return (-945/(32*3.1416*(h**9)))*r*((h**2-r**2)**2)
    else:
      return 0

  def der2SKFunction(self, r, h):
    if (0 < r) and (r < h):
      return (45/(3.1416*(h**6)))*(h - r)
    else:
      return 0

print('Build Class Particle Successfully!')

class ParticleBuffer:
  def __init__(self):
    self.particles = []
    self.numParticles = 0

  def addParticle(self, position, acc, vel, density, pressure):
    p = Particle(position, acc, vel, density, pressure, self.numParticles)
    self.particles.append(p)
    self.numParticles += 1
  
  def getParticle(self, idx):
    if (idx < self.numParticles):
      return self.particles[idx]
    else:
      print('ParticleIdx out of range!')
      return
  
  def getdistance(self, position, pidx):
    return math.sqrt((position[0]-self.particles[pidx].position[0])**2 + \
                      (position[1]-self.particles[pidx].position[1])**2 + \
                      (position[2]-self.particles[pidx].position[2])**2)
  
  def getPdistance(self, Ppos, pidx):
    return math.sqrt((Ppos[0]-self.particles[pidx].Ppos[0])**2 + \
                      (Ppos[1]-self.particles[pidx].Ppos[1])**2 + \
                      (Ppos[2]-self.particles[pidx].Ppos[2])**2)
  
  def getPDirction(self, Ppos, pidx):
    return [Ppos[0]-self.particles[pidx].Ppos[0],
            Ppos[1]-self.particles[pidx].Ppos[1],
            Ppos[2]-self.particles[pidx].Ppos[2]]
  
  def update(self, box, timestep, threshold, minIter, ext, miu, mass_rev, rou0):
    box.loadParticles(self)
    for particle in self.particles:
      neighg = box.findNeighbourGrids(particle.position)
      particle.neighp = []
      for grid in neighg:
        if box.grids[grid] != -1:
          for pidx in box.grids[grid].particles:
            if self.getdistance(particle.position, pidx) <= box.unitlength/2:
              particle.neighp.append(pidx)        #update neighbours
      if particle.idx in particle.neighp:        
        particle.neighp.remove(particle.idx)
      particle.forcee = ext                 #update external forces
      particle.updateForcef(box, 1/mass_rev)
      particle.forcep = [0.0, 0.0, 0.0]          #initialize pressure force to 0
      particle.Ppre = 0.0                  #initialize pressure to 0
      particle.forcev = [0.0, 0.0, 0.0]
      for pidx in particle.neighp:           #Fvis = 6*pi*miu*v*r miu=0.6
        r = self.getdistance(particle.position, pidx)
        v = [self.particles[pidx].vel[0]-particle.vel[0], 
             self.particles[pidx].vel[1]-particle.vel[1], 
             self.particles[pidx].vel[2]-particle.vel[2]]
        particle.forcev[0] += 6*3.14159*miu*v[0]*r
        particle.forcev[1] += 6*3.14159*miu*v[1]*r
        particle.forcev[2] += 6*3.14159*miu*v[2]*r     #update vis force
      
    iter = 0
    while (iter < minIter):
      iter += 1
      for particle in self.particles:
        # particle.acc = [mass_rev*(particle.forcev[0]+particle.forcep[0]+particle.forcee[0]), 
        #                 mass_rev*(particle.forcev[1]+particle.forcep[1]+particle.forcee[1]),
        #                 mass_rev*(particle.forcev[2]+particle.forcep[2]+particle.forcee[2])] #update acceleration
        particle.acc = [mass_rev*(particle.forcev[0]/10+particle.forcep[0]/50+particle.forcee[0]*10), 
                        mass_rev*(particle.forcev[1]/10+particle.forcep[1]/50+particle.forcee[1]*10),
                        mass_rev*(particle.forcev[2]/10+particle.forcep[2]/50+particle.forcee[2]*10)] #update acceleration
        particle.Pvel = [particle.vel[0]+timestep*particle.acc[0],
                         particle.vel[1]+timestep*particle.acc[1],
                         particle.vel[2]+timestep*particle.acc[2]]              #update velocity
        particle.Ppos = [particle.position[0]+timestep*particle.Pvel[0],
                         particle.position[1]+timestep*particle.Pvel[1],
                         particle.position[2]+timestep*particle.Pvel[2]]           #update position

      for particle in self.particles:
        rlist = []
        for pidx in particle.neighp:
          rlist.append(self.getPdistance(particle.Ppos, pidx))
        particle.updateDensity(1/mass_rev, rlist, box.unitlength/2)                             #update density
        particle.Pdenv = particle.Pden - rou0
        particle.updatePressure(1/mass_rev, timestep, rou0, rlist, box.unitlength/2)                   #update pressure

      for particle in self.particles:
        plist = []
        rlist = []
        dirlist = []
        for pidx in particle.neighp:
          rlist.append(self.getPdistance(particle.Ppos, pidx))
          plist.append(self.particles[pidx].Ppre)
          dirlist.append(self.getPDirction(particle.Ppos, pidx))
        particle.updateForcep(1/mass_rev, plist, rlist, dirlist)                           #update force from pressure

    for particle in self.particles:
      particle.position = particle.Ppos
      particle.vel = particle.Pvel
      particle.density = particle.Pden                        
      particle.pressure = particle.Ppre                                 #update the predictive values
      if (particle.position[1] < 0):
        # particle.position[1] = -particle.position[1]
        # particle.vel[1] = -particle.vel[1]
        particle.position[1] = 0
        particle.vel[1] = -particle.vel[1]*0.9
      if (particle.position[0] > 1):
        # particle.position[1] = -particle.position[1]
        # particle.vel[1] = -particle.vel[1]
        particle.position[0] = 1
        particle.vel[0] = -particle.vel[0]*0.9
      if (particle.position[2] > 1):
        # particle.position[1] = -particle.position[1]
        # particle.vel[1] = -particle.vel[1]
        particle.position[2] = 1
        particle.vel[2] = -particle.vel[2]*0.9
      if (particle.position[0] < -1):
        # particle.position[1] = -particle.position[1]
        # particle.vel[1] = -particle.vel[1]
        particle.position[0] = -1
        particle.vel[0] = -particle.vel[0]*0.9
      if (particle.position[2] < -1):
        # particle.position[1] = -particle.position[1]
        # particle.vel[1] = -particle.vel[1]
        particle.position[2] = -1
        particle.vel[2] = -particle.vel[2]*0.9
    return
print('Build Class ParticleBuffer Successfully!')

#test 1
buffer = ParticleBuffer()
buffer.addParticle([0,0,0], [1,1,1], [1,1,1], 9, 5)
buffer.addParticle([2,2,2], [3,4,5], [6,7,8], 11, 8)
if buffer.getParticle(1).acc == [3,4,5]:
  print('Pass Test1')
