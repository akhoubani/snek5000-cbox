!-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'
      include 'INPUT'
      include 'TSTEP'

      udiff = 0.
      utrans = 0.

      return
      end
!-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'
      include 'PARALLEL'
      include 'SOLN'
      include 'INPUT'
      include 'SFDD'

      integer ix, iy, iz, ieg, iel
      real rtmp, Pr_, enable_sfd

      Pr_ = abs(UPARAM(1))
      enable_sfd = abs(UPARAM(7))

      ! local element number
      iel = GLLEL(ieg)

      FFX = 0
      FFY = 0
      FFZ = 0
      ! forcing, put boussinesq
      if (IFPERT) then
         ip = ix + NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1)))
         rtmp = TP(ip,1,1)*Pr_
      else
         rtmp = T(ix,iy,iz,iel,1)*Pr_
      endif

      FFX = 0
      FFY = rtmp
      if (IF3D) FFZ = 0
      if (enable_sfd.ne.0.0) then
         call sfd_forcing(FFX,FFY,FFZ,ix,iy,iz,ieg)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol   = 0.
      source = 0.

      return
      end
!-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      include 'CHKPOINTD'  ! chpt_ifrst

      integer nxyz1, nxyz2, nit_pert, nit_hist

      real vtmp(lx1*ly1*lz1*lelt,ldim),ttmp(lx1*ly1*lz1*lelt) ! temporary variables
      real ptmp(lx2*ly2*lz2*lelt), enable_sfd

      logical exist

      ! for remaining_clock_time
      character*28 rct_file_name
      real rct_seconds_between_save
      real remaining_clock_time, remaining_eq_time, delta_clock_time

      rct_seconds_between_save = real(UPARAM(12))
      rct_file_name = 'remaining_clock_time.csv'

      nit_hist = abs(UPARAM(10))
      nit_pert = abs(UPARAM(9))
      enable_sfd = UPARAM(7)

      if (ISTEP.eq.0) then
         if (IFPERT) then
            TIME = 0
         endif
      ! start framework
         call frame_start

         inquire(file=rct_file_name, exist=exist)
         if (.not. exist) then
            if (nid.eq.0) then
               open(11, File=rct_file_name)
               write(11,'(a,a)')
     &           'it,equation_times,dt,delta_clock_times,',
     &           'remaining_equation_times,remaining_clock_times'
               close(11)
            endif
         endif

      endif

      ! monitor simulation
      call frame_monitor
      call chkpt_main
      if (enable_sfd.ne.0.0) then
         call sfd_main
      endif

      ! finalise framework
      if (istep.eq.nsteps.or.lastep.eq.1) then
         call frame_end
      else if (istep .eq. 0 .and. (.not. chpt_ifrst) ) then
      ! first timestep which is not restarted
         call outpost(vx, vy, vz, pr, t, '   ')  ! write initial condition
      endif

      ! perturbation field
      if (IFPERT) then
         if (mod(ISTEP,nit_pert).eq.0) then
      ! write perturbation field
             call out_pert()
         endif
      endif

      ! history points

      ! we have two formulations for pressure solver
      nxyz1 = nx1*ny1*nz1 ! velocity
      nxyz2 = nx2*ny2*nz2 ! pressure

      if (lhis.gt.1) then
         if (mod(ISTEP,nit_hist).eq.0) then
            if (.not. ifpert) then
               call hpts()
            else
               ! save base state in temporary variables
               call opcopy(vtmp(1,1),vtmp(1,2),vtmp(1,ndim),vx,vy,vz)
               call copy(ttmp,T,nxyz1*NELV)
               call copy(ptmp,PR,nxyz2*NELV)

               call opcopy(vx, vy, vz, vxp, vyp, vzp)
               call copy(T,TP,nxyz1*NELV)
               call copy(PR,PRP,nxyz2*NELV)

               call hpts()
               ! restore base state
               call opcopy(vx,vy,vz, vtmp(1,1),vtmp(1,2),vtmp(1,ndim))
               call copy(T,ttmp,nxyz1*NELV)
               call copy(PR,ptmp,nxyz2*NELV)

            endif
         endif
      endif

      if (nid .eq. 0) then
         call compute_remaining_clock_time(
     &       rct_seconds_between_save, remaining_clock_time,
     &       remaining_eq_time, delta_clock_time)
         if (istep == 0) then
            open(10, File=rct_file_name, position='append')
            write(10,'(I12,A,g14.8,A,g14.8,A,g14.8,A)')
     &          istep,',',time,',',dt,',',delta_clock_time,
     &          ',nan,nan'
            close(10)
         elseif (remaining_clock_time >= 0.0) then
            open(10, File=rct_file_name, position='append')
            write(10,'(I12,A,g14.8,A,g14.8,A,g14.8,A,g14.8,A,g14.8)')
     &          istep,',',time,',',dt,',',delta_clock_time,
     &          ',',remaining_eq_time,
     &          ',',remaining_clock_time
            close(10)
         endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'SOLN'
      include 'NEKUSE'

      real delta_T_side, delta_T_vert
      real xmax, ymax, dTs, dTv

      delta_T_side = abs(UPARAM(5))
      delta_T_vert = abs(UPARAM(6))

      xmax = u_oper_lx
      ymax = u_oper_ly
      if (if3d) then
            zmax = u_oper_lz
      endif

      dTs = delta_T_side/2.
      dTv = delta_T_vert/2.

      ! base flow
      if (JP.eq.0) then

         ux = 0.
         uy = 0.
         uz = 0.

         if (delta_T_side.ne.0.and.delta_T_vert.eq.0) then
            temp = delta_T_side * (x/xmax - 0.5)
         elseif (delta_T_vert.ne.0.and.delta_T_side.eq.0) then
            temp = delta_T_vert * (0.5- y/ymax)
         elseif (delta_T_side.ne.0.and.delta_T_vert.ne.0) then
            if (x.eq.0) then
               temp = -dTs
            elseif (x.eq.xmax) then
               temp = dTs
            elseif(y.eq.0) then
               temp = dTv
            elseif (y.eq.ymax) then
               temp = -dTv
            elseif (x.eq.0.and.y.eq.0) then
               temp = (dTv-dTs)/2.
            elseif (x.eq.0.and.y.eq.ymax) then
               temp = (-dTv-dTs)/2.
            elseif (x.eq.xmax.and.y.eq.0) then
               temp = (dTs+dTv)/2.
            elseif (x.eq.xmax.and.y.eq.ymax) then
               temp = (dTs-dTv)/2.
            endif
         endif

      ! perturbation
      else

         ux = 0.
         uy = 0.
         uz = 0.
         temp = 0.
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'
      include 'SOLN'
      include 'GEOM'
      include 'INPUT'

      integer, parameter :: seed = 1990
      integer :: clock
      real delta_T_side, delta_T_vert, amplitude
      real xmax, ymax, ran
      integer, allocatable :: seed_array(:)
      integer :: seed_size

      delta_T_side = abs(UPARAM(5))
      delta_T_vert = abs(UPARAM(6))
      amplitude = 1e-5

      xmax = u_oper_lx
      ymax = u_oper_ly
      if (if3d) then
            zmax = u_oper_lz
      endif

      ! base flow
      if (JP.eq.0) then

         ux = 0.0
         uy = 0.0
         uz = 0.0

         call random_number(temp)
         temp = amplitude * temp

         if (delta_T_vert.ne.0.and.delta_T_side.eq.0) then
            if (IFPERT) then
               temp = delta_T_vert * (0.5- y/ymax)
            else
               temp = delta_T_vert * (0.5- y/ymax) + temp
            endif
         elseif (delta_T_side.ne.0.and.delta_T_vert.eq.0) then
                temp = delta_T_side * (x/xmax - 0.5) + temp
         endif

      ! perturbation
      else

         ux = 0.0
         uy = 0.0
         uz = 0.0

         call random_seed(size = seed_size)
         allocate(seed_array(seed_size))

         call system_clock(count=clock)
         seed_array = seed + clock

         call random_seed(put = seed_array)
         deallocate(seed_array)

         call random_number(temp)
         temp = amplitude * temp


         endif

      return
      end
!-----------------------------------------------------------------------
! This routine to modify element vertices
      subroutine usrdat
      include 'SIZE'

      return
      end
!-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'SOLN'

      integer i, ntot
      real stretch_x, stretch_y, stretch_z, xx, yy, zz, twopi
      real xmax, ymax, zmax

      stretch_y = abs(UPARAM(4))

      if (stretch_y.ne.0.0) then
         ntot = nx1*ny1*nz1*nelt

         xmax = glmax(xm1,ntot)
         ymax = glmax(ym1,ntot)
         if (if3d) then
             zmax = glmax(zm1,ntot)
         endif

         twopi = 8 * atan(1.)

      ! stretch factors
         stretch_x = stretch_y*xmax
         ! stretch_x = 0
         if (if3d) then
             stretch_z = stretch_y*zmax
         endif

         do i=1,ntot
             xx = xm1(i,1,1,1)
             yy = ym1(i,1,1,1)
             xm1(i,1,1,1) = xx - (stretch_x * (sin(twopi*xx/xmax)))
             ym1(i,1,1,1) = yy - (stretch_y * (sin(twopi*yy/ymax)))

             if (if3d) then
                 zz = zm1(i,1,1,1)
                 zm1(i,1,1,1) = zz - (stretch_z * (sin(twopi*zz/zmax)))
             endif
         enddo
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine usrdat3

      return
      end
!-----------------------------------------------------------------------
      subroutine frame_usr_register
      implicit none
      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'

      real enable_sfd

      enable_sfd = abs(UPARAM(7))

      ! register modules
      call io_register
      call chkpt_register
      if (enable_sfd.ne.0.0) then
         call sfd_register
      endif

      return
      end subroutine
!-----------------------------------------------------------------------
      subroutine frame_usr_init
      implicit none
      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'

      real enable_sfd

      enable_sfd = abs(UPARAM(7))

      ! initialise modules
      call chkpt_init
      if (enable_sfd.ne.0.0) then
         call sfd_init
      endif

      return
      end subroutine
!-----------------------------------------------------------------------
      subroutine frame_usr_end
      implicit none
      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'

      real enable_sfd

      enable_sfd = abs(UPARAM(7))

      if (enable_sfd.ne.0.0) then
         call sfd_end
      endif

      return
      end subroutine

c-----------------------------------------------------------------------
      subroutine compute_remaining_clock_time(
     &    seconds_between_save, remaining_clock_time, remaining_eq_time,
     &    delta_clock_time)
      implicit none

      include 'SIZE'
      include 'TOTAL'

      real seconds_between_save
      real remaining_clock_time, remaining_eq_time

      real clock_time, delta_clock_time, delta_eq_time

      real clock_time_last, eq_time_last
      save clock_time_last, eq_time_last

      integer istep_next, istep_last
      save istep_next

      remaining_clock_time = -1.0
      delta_clock_time = 0.0
      if (istep .eq. 0) then
        eq_time_last = time
        istep_next = 1
      elseif ((istep .eq. istep_next) .or. (lastep .eq. 1)) then
        call cpu_time(clock_time)
        delta_clock_time = clock_time - clock_time_last
        clock_time_last = clock_time

        delta_eq_time = time - eq_time_last
        eq_time_last = time

        if (fintim > 0.0) then
          remaining_eq_time = fintim - time
        else
          remaining_eq_time = (nsteps - istep) * dt
        endif
        if (istep .eq. 1) then
          istep_next = 2
        else
          istep_next = istep + int((istep - istep_last) *
     &      seconds_between_save / delta_clock_time)
          if (istep_next .eq. istep) then
            istep_next = istep_next + 1
          endif
        endif
        istep_last = istep
        remaining_clock_time = remaining_eq_time / delta_eq_time *
     &      delta_clock_time
        if (remaining_clock_time < 0.0) then
          remaining_clock_time = 0.0
        endif
      endif

      end subroutine
