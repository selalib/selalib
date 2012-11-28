class xdmf_2D:
   def __init__(self, prefix, index=0):
      self.index = index
      self.prefix = prefix
      self._f = open("%s%05d.xmf"%(self.prefix,self.index), "w")
      (self._f).write("<?xml version=\"1.0\" ?>\n")
      (self._f).write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
      (self._f).write("<Xdmf Version=\"2.0\">\n")
      (self._f).write(" <Domain>\n")
      (self._f).write("   <Grid Name=\"Mesh\" GridType=\"Uniform\">\n");
   def grid(self,nx, ny, meshprefix="Phi_3D_d"):
      self._nx = nx
      self._ny = ny
      self._mp = meshprefix
      (self._f).write("     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n"%(self._ny,self._nx));
      (self._f).write("     <Geometry GeometryType=\"VXVY\">\n");
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"% self._nx)
      (self._f).write("        %s%05d.h5:/rg\n" % (self._mp,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"% self._ny)
      (self._f).write("        %s%05d.h5:/thetag\n" % (self._mp,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Geometry>\n")
   def field(self,name):
      assert self._nx>0 
      assert self._ny>0
      (self._f).write("     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n" % name)
      (self._f).write("       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" "%(self._ny,self._nx))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n" )
      (self._f).write("        %s%05d.h5:/%s\n" % (self.prefix,self.index,name))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Attribute>\n")
   def close(self):
      (self._f).write("   </Grid>\n")
      (self._f).write(" </Domain>\n")
      (self._f).write("</Xdmf>\n")
      (self._f).close()


class xdmf_3D(xdmf_2D):
   def grid(self,nx, ny, nz, meshprefix="Phi_3D_d"):
      self._nx = nx
      self._ny = ny
      self._nz = nz
      self._mp = meshprefix
      (self._f).write("     <Topology TopologyType=\"3DRectMesh\" ")
      (self._f).write("        NumberOfElements=\"%d %d %d\"/>\n" % (self._nz, self._ny, self._nx))
      (self._f).write("     <Geometry GeometryType=\"VXVYVZ\">\n")
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" "%(self._nx))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n" )
      (self._f).write("        %s%05d.h5:/rg\n" % (self.prefix,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" "%(self._ny))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n")
      (self._f).write("        %s%05d.h5:/thetag\n" % (self.prefix,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("       <DataItem Dimensions=\"%d\" NumberType=\"Float\" "%(self._nz))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n")
      (self._f).write("        %s%05d.h5:/phig\n" % (self.prefix,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Geometry>\n")

   def field(self,name):
      assert self._nx>0 
      assert self._ny>0
      assert self._nz>0
      (self._f).write("     <Attribute Name=\"Phi3D\" AttributeType=\"Scalar\" Center=\"Node\">\n")
      (self._f).write("       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" "%(self._nx,self._ny,self._nx))
      (self._f).write("          Precision=\"8\" Format=\"HDF\">\n" )
      (self._f).write("        %s%05d.h5:/Phi_3D\n" % (self.prefix,self.index))
      (self._f).write("       </DataItem>\n")
      (self._f).write("     </Attribute>\n")

import glob

import h5py
import numpy as np
import pylab as pl

prefix = "Phi_3D_d"

nfile = len(glob.glob1('.',prefix+"*.h5"))

f = h5py.File(prefix+"00000.h5","r")
nr = f['Nr'].value+1
ntheta = f['Ntheta'].value+1
nphi = f['Nphi'].value+1
f.close()

#phi2d = xdmf_2D("Phi2D_d")
#phi2d.grid(nr,ntheta)
#phi2d.field("Phirth")
#phi2d.close()
#
#phi3d = xdmf_3D("Phi_3D_d")
#phi3d.grid(nr,ntheta,nphi)
#phi3d.field("Phi_3D")
#phi3d.close()

for ifile in range(nfile):
   f2drthv0 = xdmf_2D("f2D_d",ifile)
   f2drthv0.grid(nr,ntheta)
   f2drthv0.field("frthv0")
   f2drthv0.close()

