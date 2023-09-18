! To compile use the command:
!
! Linux or Mac:
! gfortran -o g3r_stab g3r_stab.f95 -Ofast
!
! Windows:
! gfortran -o g3r_stab.exe g3r_stab.f95 -Ofast
!
! To run program use file g3r_input.txt with input parameters
!
! Linux or Mac:
! ./g3r_stab < g3r_input.txt
!
! Windows:
! g3r_stab.exe < g3r_input.txt
!
! For output to file use:
! g3r_stab.exe < g3r_input.txt > outpitfilename
!
PROGRAM main
    IMPLICIT REAL*8 (a-h,m,o-z)
    COMMON /aaa/    a1,a2,ene0,prd1
    COMMON /diagno/ energy
    COMMON /paskel/ hs
    COMMON /jeet/   eej(3)
    COMMON /vlight/ clight
    COMMON /problem/iproblem
    REAL*8 inc1, inc2, m(3)
    REAL*8 Q, Qmax

    READ(5,*) iwr,ns,ist,max_outer,A1,A2,clight
    READ(5,*) mass1, mass2, mass3
    READ(5,*) ax1,ecc1,inc1,m1,ome1,peri1
    READ(5,*) Qmax,ecc2,inc2,m2,ome2,peri2
    total_time = SECOND()
    m = (/mass1,mass2,mass3/)
    Q = 1.0D0
    DO WHILE (Q <= Qmax)
        duration = SECOND()
        ax2 = Q*ax1/(1D0 - ecc2)
        irev = itest (Q,iwr,ns,ist,max_outer,m,ax1,ecc1,inc1,m1,ome1,peri1,ax2,ecc2,inc2,m2,ome2,peri2)
        IF (irev >= max_outer) THEN
            WRITE (6, '("Q = ",F4.1," Stable on    ",I10," period of outer binary. Duration: ", F10.2, " s")') &
                Q, irev, SECOND() - duration
            EXIT
        END IF
        CALL FLUSH(6)
        Q = Q + 0.1D0
    END DO
    WRITE (6,'("Total time:",F10.2," s")') SECOND() - total_time
END PROGRAM MAIN

FUNCTION itest (Q,iwr,ns,ist,max_outer,m,ax1,ecc1,inc1,m1,ome1,peri1,ax2,ecc2,inc2,m2,ome2,peri2)
    IMPLICIT REAL*8 (a-h,m,o-z)
    COMMON /aaa/    a1, a2, ene0, prd1
    COMMON /diagno/ energy
    COMMON /jeet/   eej(3)
    COMMON /paskel/ hs
    COMMON /vlight/ clight
    COMMON /problem/iproblem
    PARAMETER(pi    = 3.141592653589793D0, &
              twopi = 6.283185307179586D0)
    REAL*8 ej1(3),ej2(3)
    REAL*8 xw(6), vw(6)
    REAL*8 m(3), x(6), v(6)
    REAL*8 xj(6), vj(6)
    REAL*8 inc1, inc2, inc12, inc3 
    deg=ATAN(1.D0)/45
    iproblem = 0
    duration = SECOND()
    period_ratio = DSQRT(ax2**3/ax1**3)
    tmax = max_outer*period_ratio/twopi
    m12  = m(1) + m(2)
    m3 = m12 + m(3)
    ttime= 0
    IF (ns == 0) THEN
        itest = 0
        RETURN
    END IF
    CALL inite(x(1), v(1), ax1, ecc1, inc1, m1, ome1, peri1, m12)
    CALL inite(x(4), v(4), ax2, ecc2, inc2, m2, ome2, peri2, m3)
    xj = x
    vj = v
    xw = xj(1:6)
    vw = vj(1:6)
    CALL energytest(xw, vw, m, ene0, g, a1, a2, ns, hs)
    prd0  = prd1
    DO ITE=1,10
        xw = xj(1:6)
        vw = vj(1:6)
        CALL cor(hs, xw, vw, m)
        xj(1:6) = xj(1:6)+x(1:6)-xw
        vj(1:6) = vj(1:6)+v(1:6)-vw
    END DO
    xw = xj(1:6)
    vw = vj(1:6)
    CALL cor(hs*2, xw, vw, m)
    CALL energytest(xw, vw, m, ene1, g, dum, dum, ns, dum)
    CALL energytest(xj, vj, m, ene0, g, dum, dum, ns, dum)
    prd0 = prd1
    ene0 = ene1
    CALL extr(xj, vj, ene0, hs, m)
    is     = 0
    DO
        ytime = ttime/twopi
        IF (ytime > tmax*prd0) THEN
            itest = INT(ytime/prd0*twopi/period_ratio)
            RETURN
        END IF
        IF (is == ist) THEN
            is  = 0
            xw = xj(1:6)
            vw = vj(1:6)
            CALL cor(hs*2, xw, vw, m)
            CALL energytest(xw, vw, m, enet ,g, dum, dum, ns, dum)
            dene = enet/ene0 - 1
            CALL elmnts (xw(1),vw(1),m12,a12,e12,mo12,inc12,om12,oo12,alfa12,q12,tq12)
            DO k=1,3
                ej1(k)=eej(k)
            END DO
            
            CALL elmnts (xw(4),vw(4),m3, a3, e3, mo3, inc3, om3, oo3, alfa3, q3, tq3)
            DO k=1,3
                ej2(k)=eej(k)
            END DO
            IF((a12<0).or.(a3<0).or.(a12>a3).or.(e12>1).or.(e3>1)) THEN
                WRITE (6, '("Q = ",F4.1," Instabilty on", I10, " period of outer binary. Duration: ", F10.2, " s")')&
                    Q, INT(ytime/prd0*twopi/period_ratio), SECOND() - duration
                itest = INT(ytime/prd0*twopi/period_ratio)
                RETURN
            END IF
        END IF
        is = is + 1
        CALL si3(hs, xj, vj, m, ttime, 0)
        CALL si3(hs, xj, vj, m, ttime, 1)
        CALL extr(xj, vj, ene0, hs*2, m)
        IF (iproblem > 0) THEN
            itest = INT(ytime/prd0*twopi/period_ratio)
            RETURN
        END IF
    END DO
END FUNCTION itest


SUBROUTINE cross(a, b, c)
    REAL*8 a(3), b(3), c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
END SUBROUTINE cross


FUNCTION norm(x)
    REAL*8 x(3), norm
    norm = SQRT(x(1)**2 + x(2)**2 + x(3)**2)
END FUNCTION


FUNCTION dot(a, b)
    REAL*8 a(3), b(3), dot
    dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
END FUNCTION dot

SUBROUTINE inite(x,w,aks,ecc,inc0,m0,ome0,peri0,mass)
    IMPLICIT REAL*8 (a-h,m,o-z)
    REAL*8 inc0, inc, m0
    REAL*8 x(3), w(3), a(3), b(3)
    aste = DATAN(1D0)/45
    inc  = inc0*aste
    ome  = ome0*aste
    perih= peri0*aste
    anom = m0*aste
    sini = DSIN(inc)
    cosi = DCOS(inc)
    come = DCOS(ome)
    some = DSIN(ome)
    cperi= DCOS(perih)
    speri= DSIN(perih)
    e    = anom
    DO iter = 1,20
        eccsin = ecc*DSIN(e)
        perd   = 1D0/(1 - ecc*DCOS(e))
        de     = -(e-eccsin-anom)*perd
        d2     = eccsin*perd
        de     = de/DSQRT(1+de*d2)
        e      = e+de
        IF (DABS(de) < 1D-15) EXIT
    END DO
    sine = DSIN(e)
    cose = DCOS(e)
    r    = aks*(1-ecc*cose)
    edot = 1D0/(r*DSQRT(aks ))
    a(1) = aks*(cperi*come-speri*some*cosi)
    a(2) = aks*(cperi*some+speri*come*cosi)
    a(3) = aks*speri*sini
    acfii= aks*DSQRT(DABS(1-ecc*ecc))
    b(1) = acfii*(-speri*come-cperi*some*cosi)
    b(2) = acfii*(-speri*some+cperi*come*cosi)
    b(3) = acfii*cperi*sini
    ac   = cose-ecc
    bc   = sine
    ad   = -sine*edot
    bd   = cose*edot
    sqrt_mass = DSQRT(mass)
    x    =  a*ac + b*bc
    w    = (a*ad + b*bd)*sqrt_mass
    RETURN
END SUBROUTINE inite
                
        
SUBROUTINE extr(x, v, ener0, h, m)
    IMPLICIT REAL*8 (a-h,m,o-z)
    REAL*8 x(6), v(6), dvex(6), m(3), va(6)
    COMMON /mass/   mp1, mp2, my1, my2, af0, af1, af2, ak1, ak2
    COMMON /aaa/    a1, a2, ene0, prd1
    COMMON /vlight/ clight
    IF(clight == 0) RETURN
    r1  = DSQRT(cdot(x(1),x(1)))
    r2  = DSQRT(cdot(x(4),x(4)))
    g   = r1*r2/(a1*r2 + a2*r1)
    dt  = h*g
!-----------V-JUMP DUE TO EXTERNAL PERTURBATIONS--------
    m1  = m(1)
    m2  = m(2)
    m12 = m1+m2
    c   = clight    ! 10 000.D0 ! SPEED OF LIGHT
!                   SHOULD BE ITERATED (IMPL. MIDPOINT)
    va  = v
    vv  = SUM(va**2)
    emx = 1D-20*vv
    DO iteration=1,6 ! Iteration for average v(=va)
        eva     = 0
        rdotv1  = cdot(x,va)
        rdotv2  = cdot(x(4),va(4))
        vv1     = cdot(va,va) 
        vv2     = cdot(va(4),va(4))
        CALL relaterms(x(1), va(1), r1, rdotv1, vv1, m1, m2, c, dvex)
        CALL relaterms(x(4), va(4), r2, rdotv2, vv2, m12, m(3), c, dvex(4))
        dvex    = dvex*dt
        DO k=1,6
            va_old  = va(k)
            va(k)   = v(k)+dvex(k)/2
            eva     = eva +(va(k)-va_old)**2
        END DO
    !                     iteratx=iteration
        IF(eva < emx) EXIT
    END DO          ! End of Iteration
    !-----------v-jump found--------------------------------
    dw12    = 2*cdot(v(1), dvex(1)) + cdot(dvex(1), dvex(1))
    dw22    = 2*cdot(v(4), dvex(4)) + cdot(dvex(4), dvex(4))
    v       = v + dvex
    denergy = (my1*dw12 + my2*dw22)/2
    ener0   = ener0 + denergy
    RETURN
END SUBROUTINE extr


SUBROUTINE si3(h, x, v, m, ttime, ifi)
    IMPLICIT REAL*8 (a-h,m,o-z)
    REAL*8  x(6), v(6), m(3), ga(0:5), gb(0:5), gn(5), dv1(3), dv2(3)
    REAL*8  xa(3), xb(3), dg1(3), dg2(3), dr1(3), dr2(3), dp1(3), dp2(3)
    REAL*8  oy2(5) !,dvex(6)!,ot(5)
    COMMON /aaa/    a1, a2, ener0, prd1
    COMMON /diagno/ energy
    COMMON /mass/   mp1, mp2, my1, my2, af0, af1, af2, ak1, ak2
    COMMON /problem/iproblem
    SAVE
    DATA    td, y1, y2, rest/4*0.0/, dv1, dv2/6*0.0/
    !       auxiliary quantities
    twopi = 8*atan(1D0)
    ar  = a2/a1
    mj1 = m(1) + m(2)
    mj2 = mj1 + m(3)
    r10 = DSQRT(cdot(x,x))
    r20 = DSQRT(cdot(x(4),x(4)))
    !        perturbations
    DO k=1,3
        xa(k)  = x(3+k) + ak1*x(k)
        xb(k)  = x(3+k) + ak2*x(k)
    END DO
    ra   = DSQRT(cdot(xa,xa))
    rb   = DSQRT(cdot(xb,xb))
    rest = af0/r20+af1/ra+af2/rb
    g    = r10*r20/(a1*r20 + a2*r10)
    yr13 = 1/r10**3
    yr23 = 1/r20**3
    yra3 = 1/ra**3
    yrb3 = 1/rb**3
    IF(ifi > 0) THEN
        DO k=1,3
            dg1(k)  = yr13*a1*x(k)*g**2
            dg2(k)  = yr23*a2*x(k+3)*g**2
            dr1(k)  = -af1*ak1*yra3*xa(k) - af2*ak2*yrb3*xb(k)
            dr2(k)  = -af1*yra3*xa(k) - af2*yrb3*xb(k) - af0*x(3+k)*yr23
            dp1(k)  = -dg1(k)*rest - g*dr1(k)
            dp2(k)  = -dg2(k)*rest - g*dr2(k)
            dv1(k)  = h*dp1(k)/my1*2 
            dv2(k)  = h*dp2(k)/my2*2
            v(k  )  = v(k)   + dv1(k)
            v(3+k)  = v(3+k) + dv2(k)
        END DO
    END IF
    !       two-body motions
    w12     = cdot(v(1),v(1))
    w22     = cdot(v(4),v(4))
    beta1   = 2*mj1/r10 - w12
    IF (beta1 < 0) negax = 1
    beta2   = 2*mj2/r20 - w22
    vanha   = ener2
    ener2   = -.5D0*my1*beta1 - .5D0*my2*beta2
    !       energy evaluation
    energy  = ener2 + rest
    eps     = g*(ener2 - ener0)
    eps1    = eps*a1/mp1
    eps2    = eps*a2/mp2
    !       modified masses
    mj1     = mj1*(1 + eps1)
    mj2     = mj2*(1 + eps2)
    !       modified two-body motions
    beta1   = 2*mj1/r10 - w12
    beta2   = 2*mj2/r20 - w22
    eta1    = cdot(x(1),v(1))
    eta2    = cdot(x(4),v(4))
    zeta1   = mj1 - beta1*r10
    zeta2   = mj2 - beta2*r20
    kie = kie + 1 
    DO i=5,2,-1
        oy2(i) = oy2(i-1)
    END DO
    oy2(1) = y2van
    IF (kie > 5) THEN
        y2    = 4*oy2(1) - 6*oy2(2) + 4*oy2(3) - oy2(4)
        y2van = y2
    END IF
    vi1 = 0
    vi2 = h/a2
    p1  = -1
    p2  = 1
    IF ((y2 < 0).or.(y2 > vi2)) y2 = .5D0*(vi2)
    DO ite=1,100
        y1  = (h - a2*y2)/a1
        z1  = beta1*y1*y1
        z2  = beta2*y2*y2
        CALL gfun(beta1, y1, ga)
        IF (iproblem > 0) RETURN
        CALL gfun(beta2, y2, gb)
        IF (iproblem > 0) RETURN
        r1  = r10 + eta1*ga(1) + zeta1*ga(2)
        r2  = r20 + eta2*gb(1) + zeta2*gb(2)
        t1  = r10*y1 + eta1*ga(2) + zeta1*ga(3)
        t2  = r20*y2 + eta2*gb(2) + zeta2*gb(3)
        ft  = t2 - t1
        f1  = r2 + ar*r1
        IF (ft < 0) THEN
            ak  = vi2
            pk  = p2
            vi1 = y2
            p1  = ft
        ELSE
            ak  = vi1
            pk  = p1
            vi2 = y2
            p2  = ft
        END IF
        et1 = eta1*ga(0) + zeta1*ga(1)
        et2 = eta2*gb(0) + zeta2*gb(1)
        ze1 = zeta1*ga(0) - beta1*eta1*ga(1)
        ze2 = zeta2*gb(0) - beta2*eta2*gb(1)
        f2  = et1 - et2*ar**2
        f3  = ze1 + ze2*ar**3
        dy2 = -ft/f1
        dy2 = -ft/(f1 + .5D0*dy2*f2)
        dy2 = -ft/(f1 + dy2*(.5D0*f2 + .1666666666666666667D0*f3*dy2))
        dy1 = -ar*dy2
        test= a1*DABS(dy1) + a2*DABS(dy2)
        IF (test < 1D-15*h) EXIT
        y2  = y2 + dy2
        IF ((y2 - vi1) < 0) y2 = .5D0*(vi1 + vi2)
        IF ((vi2 - y2) < 0) y2 = .5D0*(vi1 + vi2)
        IF (ite == ite/7*7) y2 = .5D0*(vi1 + vi2)
        y1  = (h - a2*y2)/a1
    END DO
    IF (ite >= 100) THEN
        iproblem = 1
        RETURN
    END IF
    y2van   = y2
    ga(5)   = ga(5) + ga(4)*dy1 + .5D0*ga(3)*dy1**2
    ga(4)   = ga(4) + ga(3)*dy1 + .5D0*ga(2)*dy1**2
    gb(5)   = gb(5) + gb(4)*dy2 + .5D0*gb(3)*dy2**2
    gb(4)   = gb(4) + gb(3)*dy2 + .5D0*gb(2)*dy2**2
    y1      = y1 + dy1
    y2      = y2 + dy2
    ga(3)   = .1666666666666666667D0*y1**3 - beta1*ga(5)
    gb(3)   = .1666666666666666667D0*y2**3 - beta2*gb(5)
    ga(2)   = .5D0*y1**2 - beta1*ga(4)
    gb(2)   = .5D0*y2**2 - beta2*gb(4)
    ga(1)   = y1 - beta1*ga(3)
    gb(1)   = y2 - beta2*gb(3)
    r1      = r10 + eta1*ga(1) + zeta1*ga(2)
    r2      = r20 + eta2*gb(1) + zeta2*gb(2)
    td      = r20*y2 + eta2*gb(2) + zeta2*gb(3)
    ttime   = ttime + td
    DO k=1,2
        IF (k ==1) THEN
            DO i=1,5
                gn(i) = ga(i)
            END DO
            r0  = r10
            rx  = r1
            m12 = mj1
        ELSE
            DO i=1,5
                gn(i) = gb(i)
            END DO
            r0  = r20
            rx  = r2
            m12 = mj2
        END IF
        f   = 1 - m12*gn(2)/r0
        g   = (td - m12*gn(3))
        df  = -m12*gn(1)/(r0*rx)
        dg  = 1 - m12*gn(2)/rx
        i1  = 1 + (k-1)*3
        i3  = i1 + 2
        DO i=i1,i3
            w0i = v(i)
            v(i)= df*x(i) + dg*v(i)
            x(i)= f*x(i) + g*w0i
        END DO
    END DO
    RETURN
END SUBROUTINE si3


FUNCTION cdot(a, b)
    IMPLICIT REAL*8 (a-h,o-z)
    REAL*8 a(3), b(3), cdot
    cdot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
    RETURN
END FUNCTION cdot


SUBROUTINE gfun(beta, y, g)
    IMPLICIT REAL*8 (a-h,o-z)
    COMMON /problem/iproblem
    REAL*8 g(0:5), c(5)
    z   = beta*y*y
    CALL cfun(z, c)
    IF (iproblem > 0) RETURN
    g(0)= 1 - z*c(2)
    s   = y
    DO I=1,5
        g(i) = s*c(i)
        s    = s*y
    END DO
    RETURN
END


SUBROUTINE cfun(z,c)
    IMPLICIT REAL*8 (a-h,o-z)
    COMMON /problem/iproblem
    PARAMETER(c6  = 1D0/6,   c132 = 1D0/132, c56  = 1D0/56,             &
              c30 = 1D0/30,  c24  = 1D0/24,  c156 = 1D0/156,            &
              c90 = 1D0/90,  c110 = 1D0/110, c16  = 1D0/16,  c8 = 1D0/8,&
              c72 = 1D0/72,  c42  = 1D0/42,  c120 = 1D0/120, u  = 1)
    REAL*8 c(5)
    DATA large/0/
    h = z
    DO k=0,10
        IF(DABS(h) < 0.1D0) GOTO 2
        h = 0.25D0*h
    END DO
    IF (large >= 10) THEN
        iproblem = large
        RETURN
    END IF
    IF (large < 10)  THEN
        large = large + 1
    END IF
2   CONTINUE
    c(4)=(u-h*(u-h*(u-h*c90/(u+h*c132))*c56)*c30)*c24
    c(5)=(u-h*(u-h*(u-h*c110/(u+h*c156))*c72)*c42)*c120
    DO i=1,k
        c3  = c6 - h*c(5)
        c2  = .5D0 - h*c(4)
        c(5)= (c(5) + c(4) + c2*c3)*c16
        c(4)= c3*(2 - h*c3)*c8
        h   = 4*h
    END DO
    c(3) = c6 - z*c(5)
    c(2) = .5D0 - z*c(4)
    c(1) = 1 - z*c(3)
    RETURN
END SUBROUTINE cfun
       
       
SUBROUTINE energytest(x, v, m, energy, g, an1, an2, ns, hs)
    IMPLICIT REAL*8 (a-h, m, o-z)
    REAL*8 x(6),v(6),m(3)
    REAL*8 xa(3),xb(3)
    COMMON /aaa/    a1, a2, dummy, prd1
    COMMON /mass/   mp1, mp2, my1, my2, af0, af1, af2, ak1, ak2
    DATA    i0 /0/
    SAVE
    !i0=0
    twopi   = 8*DATAN(1D0)
    !       auxiliary quantities
    mj1 = m(1)+m(2)
    mj2 = mj1+m(3)
    mp1 = m(1)*m(2)
    mp2 = mj1*m(3)
    my1 = mp1/mj1
    my2 = mp2/mj2
    af0 = (m(1) + m(2))*m(3)
    af1 = -m(1)*m(3)
    af2 = -m(2)*m(3)
    ak1 = m(2)/(m(1) + m(2))
    ak2 = ak1 - 1
    r10 = DSQRT(cdot(x,x))
    r20 = DSQRT(cdot(x(4),x(4)))
    !        perturbations
    !-----------------------------------------------------------v
    DO k=1,3
        xa(k) = x(3+k) + ak1*x(k)
        xb(k) = x(3+k) + ak2*x(k)
    END DO
    ra      = DSQRT(cdot(xa,xa))
    rb      = DSQRT(cdot(xb,xb))
    rest    = af0/r20 + af1/ra + af2/rb
    w12     = cdot(v(1),v(1))
    w22     = cdot(v(4),v(4))
    beta1   = 2*mj1/r10 - w12
    beta2   = 2*mj2/r20 - w22
    alf1    = beta1/mj1
    alf2    = beta2/mj2
    ax1     = 1/alf1
    prd1    = twopi*ax1*DSQRT(ax1/mj1)
    hs      = .5D0*prd1/ns 
    ener2   = -.5D0*my1*beta1 - .5D0*my2*beta2
    !          energy evaluation
    energy  = ener2 + rest
    IF (i0 == 0) THEN
        ga  = 1/(a1*alf1 + a2*MAX(0D0, alf2))
        an1 = a1*ga ! normalize a1 a2 such that new ga=1
        an2 = a2*ga
        i0  = 1
    end if
    g = 1/(a1/r10 + a2/r20)
    RETURN
END SUBROUTINE energytest
       
       
SUBROUTINE cor(h,xj,vj,m)
    IMPLICIT REAL*8 (a-h,m,o-z)
    REAL*8 x(6),m(3),dv1(3,3),dv2(3,3),xj(6),vj(6),gam(3)
    REAL*8 xa(3),xb(3),dg1(3),dg2(3),dr1(3),dr2(3),dp1(3),dp2(3)
    COMMON /aaa/    a1,a2,ener0,prd1
    COMMON /diagno/ energy
    COMMON /mass/   mp1,mp2,my1,my2,af0,af1,af2,ak1,ak2
    SAVE
    !       auxiliary quantities
    mj1 = m(1) + m(2)
    mj2 = mj1 + m(3)

    r10 = DSQRT(cdot(xj(1),xj(1)))
    r20 = DSQRT(cdot(xj(4),xj(4)))
    g0  = 1./(a1/r10 + a2/r20)
    ceps=.1D0*g0
    s   = ceps*h
    twos= s + s
    DO kh=1,3
        s   = -s
        if (kh == 3) s = 0
        x   = xj + s*vj
        r10 = DSQRT(cdot(x,x))
        r20 = DSQRT(cdot(x(4),x(4)))
        DO k=1,3
            xa(k) = x(3+k) + ak1*x(k)
            xb(k) = x(3+k) + ak2*x(k)
        END DO
        ra  = DSQRT(cdot(xa,xa))
        rb  = DSQRT(cdot(xb,xb))
        rest= af0/r20 + af1/ra + af2/rb
        g   = 1/(a1/r10 + a2/r20)
        gam(kh) = g*rest
        yr13= 1./r10**3
        yr23= 1./r20**3
        yra3= 1./ra**3
        yrb3= 1./rb**3
        DO k=1,3
            dg1(k)  = yr13*a1*x(k)*g**2
            dg2(k)  = yr23*a2*x(k+3)*g**2
            dr1(k)  = -af1*ak1*yra3*xa(k) - af2*ak2*yrb3*xb(k)
            dr2(k)  = -af1*yra3*xa(k) - af2*yrb3*xb(k) - af0*x(3+k)*yr23
            dp1(k)  = -dg1(k)*rest - g*dr1(k)
            dp2(k)  = -dg2(k)*rest - g*dr2(k)
            dv1(k,kh)= dp1(k)/my1
            dv2(k,kh)= dp2(k)/my2
        END DO
    END DO
    hc  = -h**2/24
    rdot= (gam(2) - gam(1))/twos

    DO k=1,3
        xj(k)   = xj(k)   + g0*hc*dv1(k,3)
        vj(k)   = vj(k)   - g0*hc*(dv1(k,2) - dv1(k,1))/twos + dg1(k)*rdot*hc/my1
        xj(k+3) = xj(k+3) + g0*hc*dv2(k,3)
        vj(k+3) = vj(k+3) - g0*hc*(dv2(k,2) - dv2(k,1))/twos + dg2(k)*rdot*hc/my2
    END DO
    RETURN
END SUBROUTINE cor


SUBROUTINE elmnts (x, v, m, a, e, mo, inc, om, oo, alfa, q, tq)
!   x, v, m -- input arrays
!   a -- semimajor axe
!   e -- eccentrisity
!   mo 
!   inc -- inclination
!   om
!   oo
!
!   Example:
!   call elmnts (xw(1),vw(1),m12,a12,e12,mo12,inc12,om12,oo12,alfa12,q12,tq12)
!   call elmnts (xw(4),vw(4),m3, a3, e3, mo3, inc3, om3, oo3, alfa3, q3, tq3)
!
!       Note: wrong results can be produced in exeptional situations
!       where some angles are undefined in terms of the expressions used.
!       This may happen in exactly planar, rectilinear... orbits.
!       Troubles can often be avoided by a very small "perturbation" of x and/or v.
    IMPLICIT REAL*8 (a-h, m, o-z)
    PARAMETER (rad = 57.29577951308232D0)
    REAL*8 x(3) ,w(3), v(3), inc, jx, jy, jz
    COMMON /jeet/ jx, jy, jz
    mu  = DSQRT(m)
    w   = v/mu
    r   = DSQRT(SUM(x**2))
    w2  = SUM(w**2)
    eta = SUM(x*w)
    alfa= 2/r - w2
    zeta= 1 - alfa*r
    !       areal velocity vector (jx,jy,jz)
    jx  = x(2)*w(3) - x(3)*w(2)
    jy  = x(3)*w(1) - x(1)*w(3)
    jz  = x(1)*w(2) - x(2)*w(1)
    d   = DSQRT(jx*jx + jy*jy + jz*jz)
    !       eccentricity vector (ex,ey,ez)
    ex  = w(2)*jz - w(3)*jy - x(1)/r
    ey  = w(3)*jx - w(1)*jz - x(2)/r
    ez  = w(1)*jy - w(2)*jx - x(3)/r
    e   = DSQRT(ex*ex + ey*ey + ez*ez)
    b   = DSQRT(jx*jx + jy*jy)
    inc = atn2(b,jz)*rad
    om  = atn2(jx,-jy)*rad
    oo  = atn2(ez*d, ey*jx - ex*jy)*rad
    a   = 1/alfa
    sqaf= DSQRT(DABS(alfa))
    q   = d*d/(1 + e) 
    too = oot(alfa, eta, zeta, q, e, sqaf)
    tq  = too/mu
    mo  = too*rad*sqaf**3
    RETURN
END SUBROUTINE elmnts


FUNCTION atn2 (s, c)
    IMPLICIT REAL*8 (a-h, o-z)
    PARAMETER(twopi = 6.283185307179586D0)
    atn2=ATAN2(s,c)
    IF(atn2 < 0) atn2 = atn2 + twopi
    RETURN
END FUNCTION atn2


FUNCTION oot (alfa, eta, zeta, q, e, sqaf) ! oot=pericentre time
    IMPLICIT REAL*8 (a-h, o-z)
    PARAMETER(tiny_ = 1D-18)
    SAVE
    IF(zeta > 0) THEN 
    !   ellipse (near peri), parabola or hyperbola.
        ecc = MAX(e, tiny_)
        x   = eta/ecc
        z   = alfa*x*x
        oot = x*(q + x*x*as3(z))
    ELSE
    !   upper half of an elliptic orbit.
        oot = (ATAN2(eta*sqaf,zeta)/sqaf - eta)/alfa
    END IF
    RETURN
END FUNCTION oot


FUNCTION as3(z)
    IMPLICIT REAL*8 (a-h,o-z)
    SAVE
    IF(z > 0.025D0) THEN ! elliptic
        x   = DSQRT(z)
        as3 = (DASIN(x) - x)/x**3
    ELSEIF(z < -0.025D0) THEN ! hyperbolic
        x   = DSQRT(-z)
        as3 = (DLOG(x + DSQRT(1 + x*x)) - x)/x/z
    ELSE ! pade approximant for small  |z|
        as3 = (1 + 6*(-19177*z/170280 + 939109*z*z/214552800))/         &
              (6*(1 - 7987*z/7095 + 54145*z*z/204336))
    END IF
    RETURN
END FUNCTION as3


SUBROUTINE relaterms (x, v, r, rdotv, vv, m1, m2, c, dv)
    IMPLICIT REAL*8 (a-h, m, n, o-z)
    REAL*8 n(3), x(3), v(3), dv(3)
!        Mora & Will, Phys. Rev. d 69, 104021
!        http://prola.aps.org/pdf/prd/v69/i10/e104021
    SAVE
    PARAMETER (pi  = 3.1415926535897932d0,                              &
               pi2 = 9.8696044010893586d0)
    c2  = 1/c**2
    c4  = 1/c**4
    c5  = 1/c**5
    c6  = 1/c**6
    c7  = 1/c**7
    vr  = rdotv/r
    n   = x/r
    m   = m1 + m2
    eta = m1*m2/m**2
    a1  = 2*(2 + eta)*(m/r) - (1 + 3*eta)*vv + 1.5D0*eta*vr**2
    a2  = -.75D0*(12 + 29*eta)*(m/r)**2 - eta*(3 - 4*eta)*vv**2         &
          -15D0/8*eta*(1 - 3*eta)*vr**4 + .5D0*eta*(13 - 4*eta)*(m/r)*vv&
          +(2 + 25*eta + 2*eta**2)*(m/r)*vr**2 + 1.5D0*eta*(3 - 4*eta)*vv*vr**2
    a2p5= 8D0/5*eta*(m/r)*vr*(17.d0/3*(m/r) + 3*vv)
    a3  = (16 + (1399D0/12 - 41D0/16*pi2)*eta + 71D0/2*eta*eta)*(m/r)**3&
          +eta*(20827D0/840 + 123D0/64*pi2 - eta**2)*(m/r)**2*vv        &
          -(1 + (22717D0/168 + 615D0/64*pi2)*eta + 11D0/8*eta**2 - 7*eta**3)*(m/r)**2*vr**2 &
          -.25D0*eta*(11 - 49*eta + 52*eta**2)*vv**3                    &
          +35D0/16*eta*(1 - 5*eta + 5*eta**2)*vr**6                     &
          -.25D0*eta*(75 + 32*eta - 40*eta**2)*(m/r)*vv**2              &
          -.5D0*eta*(158 - 69*eta - 60*eta**2)*(m/r)*vr**4              &
          +eta*(121 - 16*eta - 20*eta**2)*(m/r)*vv*vr**2                &
          +3D0/8*eta*(20 - 79*eta + 60*eta**2)*vv**2*vr**2              &
          -15D0/8*eta*(4 - 18*eta + 17*eta**2)*vv*vr**4
    a3p5= -8D0/5*eta*(m/r)*vr*(23D0/14*(43 + 14*eta)*(m/r)**2           &
          +3D0/28*(61 + 70*eta)*vv**2                                   &
          +70*vr**4 + 1D0/42*(519 - 1267*eta)*(m/r)*vv                  &
          +.25D0*(147 + 188*eta)*(m/r)*vr**2 - 15D0/4*(19 + 2*eta)*vv*vr**2)
    b1  = 2*(2-eta)*vr
    b2  = -.5D0*vr*((4 + 41*eta+8*eta**2)*(m/r)-eta*(15 + 4*eta)*vv     &
          +3*eta*(3 + 2*eta)*vr**2)
    b2p5= -8D0/5*eta*(m/r)*(3*(m/r)+vv)
    b3  = vr*((4 + (5849D0/840 + 123D0/32*pi2)*eta                      &
          -25*eta**2 - 8*eta**3)*(m/r)**2                               &
          +1D0/8*eta*(65 - 152*eta - 48*eta**2)*vv**2                   &
          +15D0/8*eta*(3 - 8*eta - 2*eta**2)*vr**4                      &
          +eta*(15+27*eta + 10*eta**2)*(m/r)*vv                         &
          -1D0/6*eta*(329 + 177*eta + 108*eta**2)*(m/r)*vr**2           &
          -.75D0*eta*(16 - 37*eta - 16*eta**2)*vv*vr**2)                    
    b3p5= 8D0/5*eta*(m/r)*(1D0/42*(1325 + 546*eta)*(m/r)**2             &
          +1D0/28*(313 + 42*eta)*vv**2 + 75*vr**4                       &
          -1D0/42*(205 + 777*eta)*(m/r)*vv                              &
          +1D0/12*(205 + 424*eta)*(m/r)*vr**2 -.75D0*(113 + 2*eta)*vv*vr**2)
     
    atot= a1*c2 + a2*c4 + a2p5*c5 + a3*c6 + a3p5*c7 !  collection of terms
    btot= b1*c2 + b2*c4 + b2p5*c5 + b3*c6 + b3p5*c7 ! here easy to neglect something (if u want)
    dv  = m/r**2*(n*atot + v*btot)
    RETURN
END SUBROUTINE relaterms
