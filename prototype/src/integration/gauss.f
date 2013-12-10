!!$ This comes from http://dl.acm.org, Algorithme 726 : ORTHPOL, appendices and supplements

!!$ To use those functions, READ the documentation beside and find more information 
!!$ about coefficients in paper *Algorithm xxx - ORTHPOL: A package of routines for 
!!$ generating orthogonal polynomials and Gauss-type quadrature rules* by _Walter 
!!$ Gautschi_ (here xxx is 726 in other references) formulas (1.1) to (1.3) page 2, 
!!$ and book **Numerical Mathematics** by _Alfio Quarteroni_, _Riccardo Sacco_ and 
!!$ _Fausto Saleri_ section 10.
!!$ If you just want to use Gauss-Lobatto inside Selalib, just use what is in 
!!$ sll_gausslobatto.F90 ans see selalib documentation

      subroutine gauss(n,alpha,beta,eps,zero,weight,ierr,e)
!!$ 
!!$ Given  n  and a measure  dlambda, this routine generates the n-point
!!$ Gaussian quadrature formula
!!$  
!!$     integral over supp(dlambda) of f(x)dlambda(x)
!!$ 
!!$        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!!$ 
!!$ The nodes are returned as  zero(k)=x(k) and the weights as 
!!$ weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion 
!!$ coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure 
!!$ dlambda. The routine computes the nodes as eigenvalues, and the 
!!$ weights in term of the first component of the respective normalized 
!!$ eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!!$ It uses a translation and adaptation of the algol procedure  imtql2,
!!$ Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified 
!!$ by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for 
!!$ Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!!$ routine  imtql2.
!!$ 
!!$        Input:  n - - the number of points in the Gaussian quadrature	
!!$                      formula; type integer
!!$                alpha,beta - - arrays of dimension  n  to be filled 
!!$                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!!$                      ...,n
!!$                eps - the relative accuracy desired in the nodes
!!$                      and weights
!!$ 
!!$        Output: zero- array of dimension  n  containing the Gaussian 
!!$                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
!!$                      ...,n
!!$                weight - array of dimension  n  containing the 
!!$                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!!$                ierr- an error flag equal to  0  on normal return,
!!$                      equal to  i  if the QR algorithm does not
!!$                      converge within 30 iterations on evaluating the 
!!$                      i-th eigenvalue, equal to  -1  if  n  is not in
!!$                      range, and equal to  -2  if one of the beta's is 
!!$                      negative.
!!$ 
!!$ The array  e  is needed for working space.
!!$ 
      dimension alpha(n),beta(n),zero(n),weight(n),e(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.) then
        ierr=-2
        return
      end if
      weight(1)=beta(1)
      if (n.eq.1) return
      weight(1)=1.
      e(n)=0.
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.) then
          ierr=-2
          return
        end if
        e(k-1)=sqrt(beta(k))
        weight(k)=0.
  100 continue
      do 240 l=1,n
        j=0
!!$ 
!!$ Look for a small subdiagonal element.
!!$ 
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(abs(e(m)).le.eps*(abs(zero(m))+abs(zero(m+1)))) goto 120
  110   continue
  120   p=zero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
!!$ 
!!$ Form shift.
!!$ 
        g=(zero(l+1)-p)/(2.*e(l))
        r=sqrt(g*g+1.)
        g=zero(m)-p+e(l)/(g+sign(r,g))
        s=1.
        c=1.
        p=0.
        mml=m-l
!!$ 
!!$ For i=m-1 step -1 until l do ...
!!$ 
        do 200 ii=1,mml
          i=m-ii
          f=s*e(i)
          b=c*e(i)
          if(abs(f).lt.abs(g)) goto 150
          c=g/f
          r=sqrt(c*c+1.)
          e(i+1)=f*r
          s=1./r
          c=c*s
          goto 160
  150     s=f/g
          r=sqrt(s*s+1.)
          e(i+1)=g*r
          c=1./r
          s=s*c
  160     g=zero(i+1)-p
          r=(zero(i)-g)*s +2.*c*b
          p=s*r
          zero(i+1)=g+p
          g=c*r-b
!!$ 
!!$ Form first component of vector.
!!$ 
          f=weight(i+1)
          weight(i+1)=s*weight(i)+c*f
          weight(i)=c*weight(i)-s*f
  200   continue
        zero(l)=zero(l)-p
        e(l)=g
        e(m)=0.
        goto 105
  240 continue
!!$ 
!!$ Order eigenvalues and eigenvectors.
!!$ 
      do 300 ii=2,n
        i=ii-1
        k=i
        p=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.p) goto 260
          k=j
          p=zero(j)
  260   continue
        if(k.eq.i) goto 300
        zero(k)=zero(i)
        zero(i)=p
        p=weight(i)
        weight(i)=weight(k)
        weight(k)=p
  300 continue
      print*," w3",weight
      do 310 k=1,n
        weight(k)=beta(1)*weight(k)*weight(k)
  310 continue
      return
!!$ 
!!$ Set error - no convergence to an eigenvalue after 30 iterations.
!!$ 
  400 ierr=l
      return
      end




      subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)
!!$ 
!!$ This is a double-precision version of the routine  gauss.
!!$ 
      double precision dalpha,dbeta,deps,dzero,dweigh,de,dp,dg,dr,
     *ds,dc,df,db
      dimension dalpha(n),dbeta(n),dzero(n),dweigh(n),de(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      dzero(1)=dalpha(1)
      if(dbeta(1).lt.0.d0) then
        ierr=-2
        return
      end if
      dweigh(1)=dbeta(1)
      if (n.eq.1) return
      dweigh(1)=1.d0
      de(n)=0.d0
      do 100 k=2,n
        dzero(k)=dalpha(k)
        if(dbeta(k).lt.0.d0) then
          ierr=-2
          return
        end if
        de(k-1)=dsqrt(dbeta(k))
        dweigh(k)=0.d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(dabs(de(m)).le.deps*(dabs(dzero(m))+dabs(dzero(m+1)))) 
     *      goto 120
  110   continue
  120   dp=dzero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(dzero(l+1)-dp)/(2.d0*de(l))
        dr=dsqrt(dg*dg+1.d0)
        dg=dzero(m)-dp+de(l)/(dg+dsign(dr,dg))
        ds=1.d0
        dc=1.d0
        dp=0.d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*de(i)
          db=dc*de(i)
          if(dabs(df).lt.dabs(dg)) goto 150
          dc=dg/df
          dr=dsqrt(dc*dc+1.d0)
          de(i+1)=df*dr
          ds=1.d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=dsqrt(ds*ds+1.d0)
          de(i+1)=dg*dr
          dc=1.d0/dr
          ds=ds*dc
  160     dg=dzero(i+1)-dp
          dr=(dzero(i)-dg)*ds+2.d0*dc*db
          dp=ds*dr
          dzero(i+1)=dg+dp
          dg=dc*dr-db
          df=dweigh(i+1)
          dweigh(i+1)=ds*dweigh(i)+dc*df
          dweigh(i)=dc*dweigh(i)-ds*df
  200   continue
        dzero(l)=dzero(l)-dp
        de(l)=dg
        de(m)=0.d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=dzero(i)
        do 260 j=ii,n
          if(dzero(j).ge.dp) goto 260
          k=j
          dp=dzero(j)
  260   continue
        if(k.eq.i) goto 300
        dzero(k)=dzero(i)
        dzero(i)=dp
        dp=dweigh(i)
        dweigh(i)=dweigh(k)
        dweigh(k)=dp
  300 continue
      do 310 k=1,n
        dweigh(k)=dbeta(1)*dweigh(k)*dweigh(k)
  310 continue
      return
  400 ierr=l
      return
      end
