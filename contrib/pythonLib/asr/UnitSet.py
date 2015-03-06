import logging

class UnitSet(dict):
	'''A set of Units, with each unit in terms of sub-unit states.  The order in which 
	the units are read from a file or added is tracked.  The format is
	(index, name, subStateCount, components)'''
	def __init__(self):
		self._counter=0

	@classmethod
	def readUnitStates(cls,inSubUnitStateFname):
		su=cls()
		inf=file(inSubUnitStateFname)
		for l in inf:
			toks=l.strip().split()
			if len(toks)==2:
				toks.append(toks[0])
			(name, subStateCount)=toks[0:2]
			su.add(name,int(subStateCount), tuple(toks[2:]))
		inf.close()
		logging.info('Read %d unit definitions from file %s',len(su),inSubUnitStateFname)

		return su
		
	def writeUnitStates(self, outSubUnitStateFname):
		f=file(outSubUnitStateFname,'w')
		v=sorted(self.values())
		for u in v:
			f.write('%s %d %s\n'%(u[1],u[2], ' '.join(u[3])))
			
		f.close()
		logging.info('wrote %d units to %s',len(self),outSubUnitStateFname)

	def add(self,name,subStateCount, components):
		if name in self:
			raise KeyError('%s is already defined'%name)
		
		self[name]=(self._counter, name, subStateCount, components)
		self._counter += 1

	def indexToUnit(self):
		'''@return: a dict mapping the unit index to unit'''
		d= dict([ (i,name) for (i, name, subStateCount, components) in self])
		if len(d) != len(self):
			raise Exception("units do not have a unique index: len(d) = %d and len(self) = %d"%(len(d), len(self)))
		return d 
