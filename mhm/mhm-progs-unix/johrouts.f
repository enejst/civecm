c Copyright (c) James G. MacKinnon, 1998
c
c johrouts.f:  This is a set of subroutines to estimate asymptotic
c critical values and P values for Johansen's cointegration tests.
c It is written in Fortran 77.  Simply call johval, specifying the first
c six arguments. The result comes back in the last argument.
c
c A standalone program called johdist.f is also available.
c
c These routines and the associated data files may be used freely for
c non-commercial purposes, provided that proper attribution is made.
c Please cite the paper
c
c   James G. MacKinnon, Alfred Haug, and Leo Michelis, "Numerical
c   distribution functions of likelihood ratio tests for cointegration,"
c   Journal of Applied Econometrics, 14, 1999, 563-577.      
c
c The routines and data files may not be incorporated into any book
c or computer program without the express, written consent of the author.
c
c The routines must have access to the files probs.tab and joh-#.tab
c for # = 1, 2, 3 ..., 12. As currently written, these files must be
c in the current directory or in the directory /usr/local/urcdist.
c
      subroutine johval(npr,itt,itv,nc,isave,arg,val)
c
c npr = number of restrictions (p - r)
c itt = 1 or 2 for lambda-max test or trace test
c itv = 0, 1, 2, 3, or 4 for cases 0, 1, 2, 1*, or 2*.
c nc = 1 or 2 for critical value or P value
c isave = 0 or 1. If isave=1, this call to johval has the same values of
c                 npr, itt, and itv as the last one, so data files do not
c                 need to be read again. Be careful using this!
c arg = level of test (between .0001 and .9999) if nc = 1
c     = test statistic if nc = 2
c val = critical value if nc = 1 (returned by routine)
c     = P value if nc = 2 (returned by routine)
c
      implicit real*8 (a-h,o-z)
      character*12 fnames(12)
      character*40 usrfile
      real*8 probs(221), cnorm(221), crits(221), wght(221)
      common/saveall/ probs, cnorm, crits, wght
c
      fnames(1) = 'joh-1.tab'
      fnames(2) = 'joh-2.tab'
      fnames(3) = 'joh-3.tab'
      fnames(4) = 'joh-4.tab'
      fnames(5) = 'joh-5.tab'
      fnames(6) = 'joh-6.tab'
      fnames(7) = 'joh-7.tab'
      fnames(8) = 'joh-8.tab'
      fnames(9) = 'joh-9.tab'
      fnames(10) = 'joh-10.tab'
      fnames(11) = 'joh-11.tab'
      fnames(12) = 'joh-12.tab'
c
c Check that parameters are valid.
c
      if (npr.lt.1.or.npr.gt.12) then
       write(6,*) 'Error! Number of restrictions must be between 1',
     &   ' and 12.'
       stop
      end if
      if (itt.lt.1.or.itt.gt.2) then
       write(6,*) 'Error! The only valid options for itt are 1 and 2.'
       stop
      end if
      if (itv.lt.0.or.itv.gt.4) then
       write(6,*) 'Error! The valid options for itv are 0, 1, 2, 3,',
     &   ' and 4.'
       stop
      end if
      if (nc.lt.1.or.nc.gt.2) then
       write(6,*) 'Error! The only valid options for nc are 1 and 2.'
       stop
      end if
      if (isave.ne.0) go to 999
c
c These open statements assume that the files *.tab are in the current
c directory or in /usr/local/urcdist. They could easily be changed to
c relax that assumption.
c
      open (unit=12,file='probs.tab',status='OLD',err=600)
      go to 601
  600 continue
c
c Change the next line if files are not in /usr/local/urcdist
c
      open (unit=12,file='/usr/local/urcdist/probs.tab',
     &  status='OLD',err=602)
  601 continue
      open (unit=11,file=fnames(npr),status='OLD',err=603)
      go to 604
  603 continue
c
c Change the next line if files are not in /usr/local/urcdist
c Be sure that usrfile is declared long enough.
c
      usrfile = '/usr/local/urcdist/' // fnames(npr)
      open (unit=11,file=usrfile,status='OLD',err=605)
  604 continue
c
c read data from unit 11.
c
c skip copyright line
c
      read(11,301) junk
      iskip = 0
c
c skip groups of 222 lines as necessary.
c
      if (itt.ne.1) iskip = 1110
      iskip = iskip + itv*222
      if (iskip.gt.0) then
        do i=1,iskip
        read(11,301) junk
  301   format(1a1)
        end do
      end if
c
      read(11,301) junk
c
      do i=1,221
      ii = 222 - i
      read(11,303) crits(ii), wght(ii)
  303 format(2d15.8)
      end do
      do i=1,221
      read(12,*) probs(i), cnorm(i)
      end do
c
c Skip to here if data files do not have to be read again
c
  999 continue
c
      if (nc.eq.1) then
        size = arg
        if (size.lt.0.9999d-4.or.size.gt.0.99991d0) then
          write(6,*) 'Level must be between .0001 and .9999. '
          stop
        end if
      else
        stat = arg
      end if
c
      if (nc.eq.1) then
c
        np = 11
        precrt = 2.d0
        call fcrit(probs,cnorm,crits,wght,cval,size,precrt,np,nx)
        val = cval
      else
c
        np = 11
        precrt = 2.d0
        call fpval(crits,cnorm,wght,probs,pval,stat,precrt,np,nx)
        val = pval
      end if
      close (unit=11)
      close (unit=12)
      return
  602 continue
      write(6,*) 'Sorry. File probs.tab not found.'
      write(6,*) 'Please check current directory and ',
     &  '/usr/local/urcdist.'
      stop
  605 continue
      write(6,*) 'Sorry. File ', fnames(npr), ' not found.'
      write(6,*) 'Please check current directory and ',
     &  '/usr/local/urcdist.'
      stop
      end
      subroutine fcrit(probs,cnorm,crits,wght,cval,size,precrt,np,nx)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon, 1995
c Routine to find a critical value for any specified test size.
c Uses GLS to estimate approximating regression.
c
      real*8 probs(221), cnorm(221), crits(221), wght(221)
      real*8 yvect(20),xmat(20,4),xomx(4,4),resid(20),gamma(4)
      real*8 omega(20,20), fits(20)
      diffm = 1000.d0
      imin = 0
      do i=1,221
      diff = abs(size - probs(i))
      if (diff.lt.diffm) then
         diffm = diff
         imin = i
         if (diffm.lt.1.d-6) go to 100
      end if
      end do
  100 continue
c
      nph = np/2
      nptop = 221 - nph
      if (imin.gt.nph.and.imin.lt.nptop) then
c
c imin is not too close to the end. Use np points around stat.
c
        do i=1,np
          ic = imin - nph - 1 + i
          yvect(i) = crits(ic)
          xmat(i,1) = 1.d0
          xmat(i,2) = cnorm(ic)
          xmat(i,3) = xmat(i,2)*cnorm(ic)
          xmat(i,4) = xmat(i,3)*cnorm(ic)
        end do
c
c form omega matrix
c
        do i=1,np
          do j=i,np
            ic = imin - nph - 1 + i
            jc = imin - nph - 1 + j
            top = probs(ic)*(1.d0 - probs(jc))
            bot = probs(jc)*(1.d0 - probs(ic))
            omega(i,j) = wght(ic)*wght(jc)*sqrt(top/bot)            
          end do
        end do
        do i=1,np
          do j=i,np
            omega(j,i) = omega(i,j)
          end do
        end do
c
        nx = 4
        call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &    np,nx,20,4,0)
c
c check to see if gamma(4) is needed
c
        sd4 = sqrt((ssrt/(np-nx))*xomx(4,4))
        ttest = abs(gamma(4))/sd4
        if (ttest.gt.precrt) then
          call innorz(size,anorm)
          cval = gamma(1) + gamma(2)*anorm + gamma(3)*anorm**2
     &      + gamma(4)*anorm**3
          return
        else
          nx = 3
          call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &      np,nx,20,4,1)
          call innorz(size,anorm)
          cval = gamma(1) + gamma(2)*anorm + gamma(3)*anorm**2
          return
        end if
c
c imin is close to one of the ends. Use points from imin +/- nph to end.
c
      else
        if (imin.lt.np) then
          np1 = imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            yvect(i) = crits(i)
            xmat(i,1) = 1.d0
            xmat(i,2) = cnorm(i)
            xmat(i,3) = xmat(i,2)*cnorm(i)
            xmat(i,4) = xmat(i,3)*cnorm(i)
          end do
        else
          np1 = 222 - imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            ic = 222 - i
            yvect(i) = crits(ic)
            xmat(i,1) = 1.d0
            xmat(i,2) = cnorm(ic)
            xmat(i,3) = xmat(i,2)*cnorm(ic)
            xmat(i,4) = xmat(i,3)*cnorm(ic)
          end do
        end if
c
c form omega matrix
c
        do i=1,np1
          do j=i,np1
            if (imin.lt.np) then
              top = probs(i)*(1.d0 - probs(j))
              bot = probs(j)*(1.d0 - probs(i))
              omega(i,j) = wght(i)*wght(j)*sqrt(top/bot)
            else
c
c This is to avoid numerical singularities at the upper end
c
              omega(i,j) = 0.d0
              if (i.eq.j) omega(i,i) = 1.d0
            end if
          end do
        end do
        do i=1,np1
          do j=i,np1
            omega(j,i) = omega(i,j)
          end do
        end do
c
        nx = 4
        call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &    np1,nx,20,4,0)
c
c check to see if gamma(4) is needed
c
        sd4 = sqrt((ssrt/(np1-nx))*xomx(4,4))
        ttest = abs(gamma(4)/sd4)
        if (ttest.gt.precrt) then
          call innorz(size,anorm)
          cval = gamma(1) + gamma(2)*anorm + gamma(3)*anorm**2
     &      + gamma(4)*anorm**3
          return
        else
          nx = 3
          call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &      np1,nx,20,4,1)
          call innorz(size,anorm)
          cval = gamma(1) + gamma(2)*anorm + gamma(3)*anorm**2
          return
        end if
c
      end if
      end
      subroutine fpval(crits,cnorm,wght,probs,pval,stat,precrt,np,nx)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon, 1995
c Routine to find P value for any specified test statistic.
c
      real*8 crits(221), cnorm(221), wght(221), probs(221)
      real*8 yvect(20),xmat(20,4),resid(20),gamma(4)
      real*8 omega(20,20), fits(20), xomx(4,4)
c
c find critical value closest to test statistic
c
      diffm = 1000.d0
      imin = 0
      do i=1,221
      diff = abs(stat - crits(i))
      if (diff.lt.diffm) then
         diffm = diff
         imin = i
      end if
      end do
c
      nph = np/2
      nptop = 221 - nph
      if (imin.gt.nph.and.imin.lt.nptop) then
c
c imin is not too close to the end. Use np points around stat.
c
        do i=1,np
          ic = imin - nph - 1 + i
          yvect(i) = cnorm(ic)
          xmat(i,1) = 1.d0
          xmat(i,2) = crits(ic)
          xmat(i,3) = xmat(i,2)*crits(ic)
          xmat(i,4) = xmat(i,3)*crits(ic)
        end do
c
c form omega matrix
c
        do i=1,np
          do j=i,np
            ic = imin - nph - 1 + i
            jc = imin - nph - 1 + j
            top = probs(ic)*(1.d0 - probs(jc))
            bot = probs(jc)*(1.d0 - probs(ic))
            omega(i,j) = wght(ic)*wght(jc)*sqrt(top/bot)
          end do
        end do
        do i=1,np
          do j=i,np
            omega(j,i) = omega(i,j)
          end do
        end do
c
        nx = 4
        call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &    np,nx,20,4,0)
c
c check to see if gamma(4) is needed
c
        sd4 = sqrt((ssrt/(np-nx))*xomx(4,4))
        ttest = abs(gamma(4))/sd4
        if (ttest.gt.precrt) then
          crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
     &      + gamma(4)*stat**3
          call ddnor(crfit,pval)
          return
        else
          nx = 3
          call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &      np,nx,20,4,1)
          crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
          call ddnor(crfit,pval)
          return
        end if
      else
c
c imin is close to one of the ends. Use points from imin +/- nph to end.
c
        if (imin.lt.np) then
          np1 = imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            yvect(i) = cnorm(i)
            xmat(i,1) = 1.d0
            xmat(i,2) = crits(i)
            xmat(i,3) = xmat(i,2)*crits(i)
            xmat(i,4) = xmat(i,3)*crits(i)
          end do
        else
          np1 = 222 - imin + nph
          if (np1.lt.5) np1 = 5
          do i=1,np1
            ic = 222 - i
            yvect(i) = cnorm(ic)
            xmat(i,1) = 1.d0
            xmat(i,2) = crits(ic)
            xmat(i,3) = xmat(i,2)*crits(ic)
            xmat(i,4) = xmat(i,3)*crits(ic)
          end do
        end if
c
c form omega matrix
c
        do i=1,np1
          do j=i,np1
            if (imin.lt.np) then
              top = probs(i)*(1.d0 - probs(j))
              bot = probs(j)*(1.d0 - probs(i))
              omega(i,j) = wght(i)*wght(j)*sqrt(top/bot)
            else
c
c This is to avoid numerical singularities at the upper end
c
              omega(i,j) = 0.d0
              if (i.eq.j) omega(i,i) = 1.d0
            end if
          end do
        end do
        do i=1,np1
          do j=i,np1
            omega(j,i) = omega(i,j)
          end do
        end do
c
        nx = 4
        call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &    np1,nx,20,4,0)
c
c check to see if gamma(4) is needed
c
        sd4 = sqrt((ssrt/(np1-nx))*xomx(4,4))
        ttest = abs(gamma(4))/sd4
        if (ttest.gt.precrt) then
          crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
     &      + gamma(4)*stat**3
          call ddnor(crfit,pval)
        else
          nx = 3
          call gls(xmat,yvect,omega,gamma,xomx,fits,resid,ssr,ssrt,
     &      np1,nx,20,4,1)
          crfit = gamma(1) + gamma(2)*stat + gamma(3)*stat**2
          call ddnor(crfit,pval)
        end if
c
c check that nothing crazy has happened at the ends
c
        if (imin.eq.1.and.pval.gt.probs(1)) pval = probs(1)
        if (imin.eq.221.and.pval.lt.probs(221)) pval = probs(221)
        return
      end if
      end
      subroutine gls(xmat,yvect,omega,beta,xomx,fits,resid,ssr,ssrt,
     &  nobs,nvar,nomax,nvmax,ivrt)
c
c Copyright (c) James G. MacKinnon, 1995
c Subroutine to do GLS estimation the obvious way
c Use only when sample size is small (nobs <= 50)
c 1995-1-3
c
      implicit real*8 (a-h,o-z)
      real*8 xmat(nomax,nvmax), yvect(nomax), omega(nomax,nomax)
      real*8 beta(nvmax), xomx(nvmax,nvmax), fits(nomax), resid(nomax)
      real*8 xomy(50)
c
c xomx is covariance matrix of parameter estimates if omega is truly known
c First, invert omega matrix if ivrt=0. Original one gets replaced.
c
      if (ivrt.eq.0) call cholx(omega,nomax,nobs,kxx)
c
c form xomx matrix and xomy vector
c
      do j=1,nvar
        xomy(j) = 0.d0
        do l=j,nvar
          xomx(j,l) = 0.d0
        end do
      end do
c
      do 21 i=1,nobs
      do 21 k=1,nobs
      do 24 j=1,nvar
      xomy(j) = xomy(j) + xmat(i,j)*omega(k,i)*yvect(k)
      do 24 l=j,nvar
      xomx(j,l) = xomx(j,l) + xmat(i,j)*omega(k,i)*xmat(k,l)
   24 continue
   21 continue
c
      do j=1,nvar
        do l=j,nvar
          xomx(l,j) = xomx(j,l)
        end do
      end do
c
c invert xomx matrix
c
      call cholx(xomx,nvmax,nvar,kxx)
c
c  now form estimates of beta.
c
      do 5 i=1,nvar
      beta(i) = 0.d0
      do 5 j=1,nvar
      beta(i) = beta(i) + xomx(i,j)*xomy(j)
    5 continue
c
c find ssr, fitted values, and residuals
c
      ssr = 0.d0
      do i=1,nobs
        fits(i) = 0.d0
        do j=1,nvar
          fits(i) = fits(i) + xmat(i,j)*beta(j)
        end do
        resid(i) = yvect(i) - fits(i)
        ssr = ssr + resid(i)**2
      end do
c
c find ssr from transformed regression
c
      ssrt = 0.d0
      do i=1,nobs
        do k=1,nobs
          ssrt = ssrt + resid(i)*omega(k,i)*resid(k)
        end do
      end do
c
      return
      end
      subroutine cholx(amat,m,n,kxx)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon, 1993
c This routine uses the cholesky decomposition to invert a real
c symmetric matrix.
c
      real*8 amat(m,m)
      kxx = 0
      do 8 i=1,n
      kl = i - 1
      do 7 j=i,n
      if (i.gt.1) then
        do 3 k=1,kl
    3   amat(i,j) = amat(i,j) - amat(k,i)*amat(k,j)
      else
        if (amat(i,i).le.0.d0) then
        kxx = i
        go to 20
      end if
      end if
      if (i.eq.j) then
        amat(i,i) = dsqrt(amat(i,i))
      else
        if (j.eq.i+1) ooa = 1.d0/amat(i,i)
        amat(i,j) = amat(i,j)*ooa
      end if
    7 continue
    8 continue
      do 13 i=1,n
      do 12 j=i,n
      ooa = 1.d0/amat(j,j)
      if (i.ge.j) then
        t = 1.d0
        go to 12
      end if
      kl = j - 1
      t = 0.d0
      do 11 k=i,kl
   11 t = t - amat(i,k)*amat(k,j)
   12 amat(i,j) = t*ooa
   13 continue
      do 16 i=1,n
      do 15 j=i,n
      t = 0.d0
      do 14 k=j,n
   14 t = t + amat(i,k)*amat(j,k)
      amat(i,j) = t
   19 amat(j,i) = t
   15 continue
   16 continue
   20 return
      end
      subroutine ddnor(ystar,gauss)
      implicit real*8(a-h,o-z)
c
c Copyright (c) James G. MacKinnon, 1993
c Routine to evaluate cumulative normal distribution
c Written originally in late 1970's
c Modified 1993 to avoid changing the argument
c
c This subroutine uses Cody's method to evaluate the cumulative
c normal distribution. It is probably accurate to 19 or 20
c significant digits. It was written in 1977, based on the Cody
c article referred to in the documentation for IMSL subroutine mdnor.
c
      real*8 p(6), q(5), a(9), b(8), c(5), d(4)
      data p(1)/-6.58749161529837803157d-04/,
     1     p(2)/-1.60837851487422766278d-02/,
     2     p(3)/-1.25781726111229246204d-01/,
     3     p(4)/-3.60344899949804439429d-01/,
     4     p(5)/-3.05326634961232344035d-01/,
     5     p(6)/-1.63153871373020978498d-02/
      data q(1)/2.33520497626869185443d-03/,
     1     q(2)/6.05183413124413191178d-02/,
     2     q(3)/5.27905102951428412248d-01/,
     3     q(4)/1.87295284992346047209d00/,
     4     q(5)/2.56852019228982242072d00/
      data a(1)/1.23033935479799725272d03/,
     1     a(2)/2.05107837782607146532d03/,
     2     a(3)/1.71204761263407058314d03/,
     3     a(4)/8.81952221241769090411d02/,
     4     a(5)/2.98635138197400131132d02/,
     5     a(6)/6.61191906371416294775d01/,
     6     a(7)/8.88314979438837594118d00/,
     7     a(8)/5.64188496988670089180d-01/,
     8     a(9)/2.15311535474403846343d-08/
      data b(1)/1.23033935480374942043d03/,
     1     b(2)/3.43936767414372163696d03/,
     2     b(3)/4.36261909014324715820d03/,
     3     b(4)/3.29079923573345962678d03/,
     4     b(5)/1.62138957456669018874d03/,
     5     b(6)/5.37181101862009857509d02/,
     6     b(7)/1.17693950891312499305d02/,
     7     b(8)/1.57449261107098347253d01/
      data c(1)/3.209377589138469472562d03/,
     1     c(2)/3.774852376853020208137d02/,
     2     c(3)/1.138641541510501556495d02/,
     3     c(4)/3.161123743870565596947d00/,
     4     c(5)/1.857777061846031526730d-01/
      data d(1)/2.844236833439170622273d03/,
     1     d(2)/1.282616526077372275645d03/,
     2     d(3)/2.440246379344441733056d02/,
     3     d(4)/2.360129095234412093499d01/
      data orpi/.5641895835477562869483d0/,
     1   root2/.70710678118654752440083d0/
c
      isw = 1
      y = ystar
      if (ystar.lt.-16.d0) y = -16.d0
      if (ystar.gt.16.d0) y = 16.d0
      x = -y*root2
      if(x.gt.0.d0) go to 1
      if(x.lt.0.d0) go to 2
      gauss = .5d0
      return
    2 continue
      x = - x
      isw = -1
    1 continue
      if(x.lt..477d0) go to 10
      if(x.le.4.d0) go to 20
c
c  evaluate erfc for x.gt.4.0
c
      x2 = x*x
      xm2 = 1.d0/x2
      xm4 = xm2*xm2
      xm6 = xm4*xm2
      xm8 = xm4*xm4
      xm10 = xm6*xm4
      top = p(1) + p(2)*xm2 + p(3)*xm4 + p(4)*xm6 + p(5)*xm8 + p(6)*xm10
      bot = q(1) + q(2)*xm2 + q(3)*xm4 + q(4)*xm6 + q(5)*xm8 + xm10
      crap = orpi + top/(bot*x2)
      erfc = dexp(-x2)*crap/x
c
      if(isw.eq.-1) erfc = 2.d0 - erfc
      gauss = erfc*.5d0
      return
   20 continue
c
c  evaluate erfc for .477.lt.x.le.4.0
c
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x3*x2
      x6 = x3*x3
      x7 = x3*x4
      x8 = x4*x4
      top = a(1) + a(2)*x + a(3)*x2 + a(4)*x3 + a(5)*x4 + a(6)*x5 +
     &  a(7)*x6 + a(8)*x7 + a(9)*x8
      bot = b(1) + b(2)*x + b(3)*x2 + b(4)*x3 + b(5)*x4 + b(6)*x5 +
     &  b(7)*x6 + b(8)*x7 + x8
      erfc = dexp(-x2)*top/bot
c
      if(isw.eq.-1) erfc = 2.d0 - erfc
      gauss = erfc*.5d0
      return
   10 continue
c
c  evaluate erf for x.lt..477
c
      x2 = x*x
      x4 = x2*x2
      x6 = x4*x2
      x8 = x4*x4
      top = c(1) + c(2)*x2 + c(3)*x4 + c(4)*x6 + c(5)*x8
      bot = d(1) + d(2)*x2 + d(3)*x4 + d(4)*x6 + x8
      erf = x*top/bot
c
      erf = erf*isw
      erfc = 1.d0 - erf
      gauss = erfc*.5d0
      return
      end
      subroutine innorz(prob,anorm)
      implicit real*8 (a-h,o-z)
c
c Copyright (c) James G. MacKinnon, 1995
c Inverse normal routine that adjusts crude result twice.
c It seems to be accurate to about 14 digits.
c Crude result is taken from Abramowitz & Stegun (1968)
c It should have abs. error < 4.5 * 10^-4
c
      data c0/2.515517d0/, d1/1.432788d0/, c1/0.802853d0/
      data c2/0.010328d0/, d3/0.001308d0/, d2/0.189269d0/
      data const/.398942280401432678d0/
      if (prob.lt.0.d0.or.prob.gt.1.d0) then
         write(6,*) 'Attempt to find inverse normal of ', prob
         stop
      end if
      pr = prob
      if (prob.gt.0.5d0) pr = 1.d0 - prob
      arg = 1/pr**2
      t = sqrt(log(arg))
      anorm = t - (c0 + c1*t + c2*t**2)/
     &  (1 + d1*t + d2*t**2 + d3*t**3)
c
c now correct crude result by direct method
c
      call ddnor(anorm,prob2)
      pr2 = 1.d0 - prob2
      arg = 1/pr2**2
      t = sqrt(log(arg))
      anorm2 = t - (c0 + c1*t + c2*t**2)/
     &  (1 + d1*t + d2*t**2 + d3*t**3)
      anorm = anorm + anorm - anorm2
      if (prob.lt.0.5d0) anorm = -anorm
c
c now correct better result, using Taylor series approximation
c
      call ddnor(anorm,prob2)
      error = prob2 - prob
      dens = const*dexp(-.5d0*anorm**2)
      anorm = anorm - error/dens
      return
      end
