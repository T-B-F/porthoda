#!/usr/bin/python3.2 -tt

import struct

#Domain Match Probability Matrix
class DMPM:
	
	
	def __init__(self, in_f):
		
		self.n_domains=0
		self.n_vals=0
		with open(in_f, "rb") as in_F:
			self.name=struct.unpack('<10s', in_F.read(10))[0].decode()
			self.n_domains=struct.unpack('<i', in_F.read(4))[0]
			self.n_vals=struct.unpack('<i', in_F.read(4))[0]
			ids = struct.unpack('<i'*self.n_domains, in_F.read(4*self.n_domains))
			self.ids={}
			for i in range(0,len(ids)):
				self.ids[ids[i]]=i
			self.row_ids = struct.unpack('<i'*self.n_domains, in_F.read(4*self.n_domains))
			self.col_ids = struct.unpack('<i'*self.n_vals, in_F.read(4*self.n_vals))
			self.values  = struct.unpack('<h'*self.n_vals, in_F.read(2*self.n_vals))



	def get_val(self, i, j):
		id1=-1
		id2=-1
		if i<j:
			id1=self.ids[i]
			id2=self.ids[j]
		else:
			id1=self.ids[j]
			id2=self.ids[i]

		id_col=self.row_ids[id1];
		if id_col != (self.n_domains-1):
			end=self.row_ids[id1+1]
		else:
			end=self.n_vals
		for i in range(id_col, end):
			if self.col_ids[i]==id2:
				return self.values[i]
		return -1;

def main():
	x=DMPM("/data/ckeme_01/projects/mda/dmpm_pfam27_crs.dat")
	print(x.get_val(483,1704))


if __name__ == "__main__":
    main()