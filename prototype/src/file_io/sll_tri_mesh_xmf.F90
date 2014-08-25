#ifndef SORTIE_RESXDMF_HPP
#define SORTIE_RESXDMF_HPP
#include "RNM.hpp"
#include "typedef.hpp"
#include "ishNS.hpp"
#include "state2D.hpp"
#include <string>

/**
   sortie fichier format XDMF pour VisIt
*/


template<class D>
void sortie_resXDMF(const char* path, const SVU<D>& X, const R t,const int n, const Rn_name* valNode=0)
{

  char name[256];
  sprintf(name,"%s.%s",path,"xmf");
  std::ofstream file(name);
  file.precision(15);
  file << "<?xml version=\"1.0\"?>" << endl;
  file <<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> "<< endl;
  file << "<Xdmf Version=\"2.0\">" << endl;
  file << "<Domain>" << endl;
  file << "<Grid Name=\"Mesh\" GridType=\"Uniform\" >" << endl;
  file << "<Topology Type=\"Quadrilateral\" NumberOfElements=\""<<X.mesh()->nm()<<"\">" << endl;
  file << "<DataItem Name=\"Connections\" Format=\"XML\" DataType=\"Int\" Dimensions=\"";
  file << X.mesh()->nm()<<" 4\">" << endl;
  for (int i=0; i<X.mesh()->nm(); ++i)
  {
     for (int j=0; j<X.mesh()->ncote(i); ++j)
        file << X.mesh()->elementConnect(i).inode()[j] << ' ';
     file << endl;
  }
  file << "</DataItem>"<< endl;
  file << "</Topology>"<< endl;
  file << "<Geometry GeometryType=\"XY\">" << endl;
  file << "<DataItem Format=\"XML\" Dimensions=\"";
  file << X.mesh()->ns() << " 2\">" << endl;
  for (int i=0; i<X.mesh()->ns(); ++i)
     file << X.mesh()->sommet()[i].x() << ' ' 
          << X.mesh()->sommet()[i].y() << endl;
  file << "</DataItem>"<< endl;
  file << "</Geometry>"<< endl;

  // Variable contenue dans le vecteur
  if (X.Values().N()>1)
  {
     int ii=0;
     while(ii<X.Values().N())
     {
        file << "<Attribute Name=\"Cell centered "
             << X.name(ii) << "\" Center=\"Cell\">" << endl;
        file << "<DataItem Format=\"XML\" Dimensions=\""
             << X.mesh()->nm() << "\">" << endl;
        for (int i=0; i<X.mesh()->nm(); ++i)
           file << X.Values()(ii,i) <<  " ";
        ii++;
        file << endl;
        file << "</DataItem>" << endl;
        file << "</Attribute>" << endl;
     }
  }


  // Variable aux sommets
  if (valNode!=0)
  {
      //cerr << " valNode " << valNode->N() << " ns " << X.mesh()->ns() << endl;
      if (valNode->N()==X.mesh()->ns())
      {
          file << "<Attribute Name\"Node centered "<< valNode->name()
               << " Center=\"Nodes\">" << endl;
          file << "<DataItem Format=\"XML\" Dimensions=\"" << X.mesh()->ns() 
               << "\">" << endl;
          for (int i=0; i<X.mesh()->ns(); ++i) file << (*valNode)(i) << " ";
          file << "</DataItem>" << endl;
          file << "</Attribute>" << endl;
      }
  }

  file << "</Grid>" << endl;
  file << "</Domain>" << endl;
  file << "</Xdmf>" << endl;

}
#endif // SORTIE_RESXDMF_HPP
