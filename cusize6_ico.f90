      program cusize6
      
      !use netcdf

      implicit none
      
      include "netcdf.inc"

      include 'parconst'
      include 'comvar'
      include 'comnc'

      integer f,itim,itimold,nx,nk,i,v,k,yyyy,mm,dd,temp3(cmax,vmax)
      INTEGER  STATUS, timID,xID,varid,start1(1),count1(1), START2(2), COUNT2(2),start3(2),count3(2)
      character fname_var*75, fname_grid*75
      real temp0,temp1(cmax),temp2(vmax,cmax),tempzf(kmax),tempzh(cmax,kmaxp1),xx0(3,cmax),xx1(3,cmax),xx2(3,cmax),test_area,indcentre4(cmax,4)
      real*8 time,dayfrac


      integer, external :: JD, GDATE
      real,parameter :: radius = 6371000


      !--- open nc files ---
      !fname_var = 'lasttimestep.nc'
      !fname_var = '20131220_TA_lowhalf.nc'
      !fname_var = 'dom4_time16.nc'
      !fname_var = 'HDCP2TA_NarvalDom2_NestE-P_01_HLD_DOM03_ML_20131220T120000Z.nc'
      fname_var = 'TA_dom4_qc_time41.nc'
      print *,'which input file (ql)?'
      !read(5,*) fname_var
      print *,'opening file: ',fname_var
      STATUS = nf_open(fname_var, nf_nowrite, ncidql)
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      !fname_grid = "subsubgrid.nc"
      fname_grid = "NarvalDom2_NestE-R02B14_DOM03.nc"
      !fname_grid = 'NarvalDom2_NestE-R02B13_DOM02.nc'
      !fname_grid = 'NarvalDom2_NestE-R02B12_DOM01.nc'
      print *,'which input file (grid)?'
      !read(5,*) fname_var
      print *,'opening file: ',fname_grid
      STATUS = nf_open(fname_grid, nf_nowrite, ncidgrid)
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)


      !--- get dimensions and grid ---

      ! load time dimension
      !status = nf_inq_dimid(ncidql, "time", xID)
      !if (status /= nf_noerr) call handle_err(status)
      !status = nf_inq_dimlen(ncidql, xID, ntim)
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      ntim =1

      print *,'nr of timepoints in netcdf file: ',ntim
      print *,''
      start1 = (/ 1    /)
      count1 = (/ ntim /)

      !status = NF_GET_VARA_DOUBLE(ncidql, xID, start1, count1, time)    !yyyymmdd.fraction
      !if (status /= nf_noerr) call handle_err(status)

      time = 1

      print *, 'time:',time


      ! load vertical dimension
      !status = nf_inq_dimid(ncidgrid, "height", xID)
      !if (status /= nf_noerr) call handle_err(status)
      !status = nf_inq_dimlen(ncidgrid, xID, nk)
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      !nk = 38
      nk = 150

      print *,'nr of levels in netcdf file: ',nk
      print *,'nr of levels in executable : ',kmax
      if (nk.ne.kmax) then
        print *,'WARNING!! vertical grid dimension in &
             executable and file do not match - stop'
        stop
      endif 
      print *,''


      !  load nr of grid cells
      !status = nf_inq_dimid(ncidgrid, "ncells", xID)
      !if (status /= nf_noerr) call handle_err(status)
      !status = nf_inq_dimlen(ncidgrid, xID, nx)
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      !nx = 2830120
      nx = 2846164
      !nx = 1253432

      print *,'nr of points on x-axis in netcdf file: ',nx
      print *,'nr of grid cells in executable       : ',cmax
      if (nx.ne.cmax) then
        print *,'WARNING!! horizontal grid dimension in&
            executable and file do not match - stop'
        stop
      endif 
      print *,''

           

      !  load cell midpoint locations
      start1 = (/ 1    /)
      count1 = (/ cmax /)
    
      status = NF_INQ_VARID(ncidgrid, 'clon', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL(ncidgrid, varid, start1, count1, temp1)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        clon = temp1
      endif



      status = NF_INQ_VARID(ncidgrid, 'clat', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL(ncidgrid, varid, start1, count1, temp1)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        clat = temp1
      endif



      !  load cell areas
      status = NF_INQ_VARID     (ncidgrid, 'cell_area', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL (ncidgrid, varid, start1, count1, temp1)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        cell_area = temp1
      endif


      darea = maxval(cell_area)   !search for the maximum gridbox area, to define the spacing of the size-histograms
      totalarea = SUM(cell_area)
 
      print *, 'darea:',darea
      print *, 'total area:',totalarea

      !  load heights 
      start1 = (/  1    /)
      count1 = (/  kmax /)



      !  z_mc(height, ncells)  geometric height at full level center "m"
      !    used in variables  u,v,temp,theta_v,qv,qc,qi,qr,qs,clc,rho,pres
      count2 = (/ cmax, kmax /)
      !status = NF_INQ_VARID     (ncidql, 'height', varid) 
      status = NF_INQ_VARID     (ncidgrid, 'height', varid) 
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL (ncidgrid, varid, start1, count1, tempzf)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        !z_f = tempzf
      endif


      !load cell vertex locations
      start2 = (/ 1   , 1    /)
      count2 = (/ vmax, cmax /)
     
      status = NF_INQ_VARID(ncidgrid, 'clon_vertices', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL (ncidgrid, varid, start2, count2, temp2)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        cell_vx = temp2
      endif


      status = NF_INQ_VARID(ncidgrid, 'clat_vertices', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL (ncidgrid, varid, start2, count2, temp2)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        cell_vy = temp2
      endif


      start3 = (/ 1   , 1    /)
      count3 = (/ cmax, vmax /)
 
      status = NF_INQ_VARID(ncidgrid, 'neighbor_cell_index', varid)
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_INT (ncidgrid, varid, start3, count3, temp3)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        cell_neighbor = temp3
      endif

      print *,'closing NetCDF file grid'
      print *,'ncidgrid:',ncidgrid

      status = NF_CLOSE(ncidgrid)
      if (status /= nf_noerr) call handle_err(STATUS)

      !--- user input ---
      call userinput


      !--- reset the size densities ---
      call inithistogram


      !################# main time loop ####################
      !
      !status = NF_INQ_VARID (ncidql, 'time', varid)
      !if (status /= nf_noerr) call handle_err(status)

      itimold = 1

      do f=1,ntim
      
        !--- read time --- 
        start1 = (/ f    /)
        count1 = (/ 1    /)
     
        status = NF_GET_VARA_DOUBLE(ncidql, varid, start1, count1, time)    !yyyymmdd.fraction
        if (status /= nf_noerr) call handle_err(status)

        time = 1

        dayfrac = time - floor(time)
        yyyy    = floor(  time                         /10000 ) 
        mm      = floor( (time - yyyy*10000)           /100   )
        dd      = floor(  time - yyyy*10000 - mm * 100        )
        jdate   = JD( yyyy, mm, dd ) + dayfrac  !convert (proleptic) Gregorian to Julian date 

        if (f.eq.1) then
          !jdate0 = jdate
          jdate0 = JD( yyyy, mm, dd ) + floor( dayfrac * 24. * 3600. / dtav ) * ( dtav / (24. * 3600.) )
        endif
        !print *,'reading field:',f,'  t=',jdate,' t0=',jdate0, dayfrac, yyyy, mm, dd

        seconds = (jdate - jdate0) * 3600. * 24.   !seconds since initialization

        !--- determine time index ---
        itim = 1 + floor(seconds/dtav)

        !--- time-averaging procedures ---
        if (itim.gt.itimold) then
          call finishhistogram(itimold)
          call inithistogram
        endif

        !--- process field ---
        print *, 'reading field:',f,'julian date =',jdate
        print *, 'seconds since init = ',seconds
        print *, ' itim',itim


        call nextfield(f)

        print *,'closing NetCDF file ql 2'
        status = NF_CLOSE(ncidql)
        if (status /= nf_noerr) call handle_err(STATUS)

        call slabstat

        call cluster(indcentre4)

        if (ncloud.gt.0) then
          call nnspacing(indcentre4)
        endif

        !--- update time index ---
        itimold = itim

      enddo

      !#####################################################
      

      call finishhistogram(itim)


      !--- close nc files ---
      !print *,'closing NetCDF file ql'
      !status = NF_CLOSE(ncidql)
      !if (status /= nf_noerr) call handle_err(STATUS)

      print *,'closing NetCDF file out'
      status = NF_CLOSE(ncidout)
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      


      print *,''



      end
     
     
!
!-----------------------------------------------------------------
!								  
      subroutine userinput
      
      include 'parconst'
      include 'comvar'

      integer n
      
      
      print *,'which size-def?'
      print *,'  1 - height               '
      print *,'  2 - sqrt(height-av.area) '
      print *,'  3 - sqrt(proj.area)      '
      print *,'  4 - (volume)**(1./3.)    '
      !read(5,*) isize
      isize = 3
      
      select case (isize)
      case (1)
        sizmin =      dz
        sizmax = kmax*dz
        sizname = 'cloudheight'
      case (2)
        sizmin =       darea **(0.5)
        sizmax = (cmax*darea)**(0.5)
        sizname = 'height-averaged area'
      case (3)
        sizmin =       darea **(0.5)
        sizmax = (cmax*darea)**(0.5)
        sizname = 'projected area'
      case (4)
        sizmin = (          darea*dz)**(1./3.)
        sizmax = (kmax*cmax*darea*dz)**(1./3.)
        sizname = 'volume'
      end select
      
      
      print *,'sizemin=',sizmin
      print *,'sizemax=',sizmax
      print *,'  '
      
      nbin = ceiling(sizmax/sizmin)
      print *,'  nr of bins =',nbin
      print *,' '
      print *,' '
      
      !--- compute size for bins ---
      do n=1,nbin
        hsize(n) = sizmin * n
      enddo
      
      
      print *,'vertical profiles with cloud type - options'
      print *,'which overlap style?'
      print *,'  1 - equalize cloudbases'
      print *,'  2 - equalize cloudtops'
      print *,'  3 - no shifting (vertical position preserved)'
      
      !read(5,*) ioverlap
      print *,'  '
      ioverlap = 3      

      
      print *,'cloud selection criteria'
      print *,'which clouds should be diagnosed?'
      print *,'  1 - all clouds'
      print *,'  2 - those rooted below LFC'
      print *,'  3 - those rooted below LFC and extending above LFC'
      
      !read(5,*) iselect
      print *,'  '
      iselect = 1
           
     
      return
      end
     
!
!-----------------------------------------------------------------
!								  
      subroutine inithistogram
      
      include 'parconst'
      include 'comvar'

      integer n,k
          
      print *,'initializing the histograms...'
      print *,' '

      !--- reset nr of fields ---
      nfield = 0
      
      !--- reset bin props ---
      do n=1,binmax
      
        hn(n)   = 0.
        hac(n)  = 0.
        hwc(n)  = 0.
        hmc(n)  = 0.
        hqlc(n) = 0.
        hnns_mean(n) = 0.
        hnns_stdev(n) = 0.
        D0(n)   = 0.
        D1(n)   = 0.
        
        do k=1,kmax
          hnlev(n) = 0.
          hthl(n) = 0.
          hthv(n) = 0.
          hqt(n) = 0.
          hql(n) = 0.
          hw(n) = 0.
          ha(n) = 0.
          hm(n) = 0.
       enddo
      enddo
    
 
      return
      end
      
                  
      
      
!
!-----------------------------------------------------------------
!								  
      subroutine finishhistogram(itim)
      
      implicit none

      include 'parconst'
      include 'comvar'
      include 'comnc'

      include "netcdf.inc"

      integer n,k,slen,itim
      character ext*10,sform1*20,sform2,ext4*4
      real ndum,mdum,adum,lcloudmax,nlcloudmax

      integer status,dd(3),start1(1),count1(1),start2(2),count2(2),start3(3),count3(3)
      character fname_out*75

      LOGICAL,SAVE::FIRST_CALL=.TRUE.
            
      !print *,'finishing the histograms...'

      if (nfield.gt.0) then

      ndum=0
      !mdum=0
      adum=0
     
      do n=1,nbin
        
  
        if (isnan(hn(n))) then
          hn(n) = 0.
        endif

        !--- sub-ensemble-averages (div by nr of clouds in sub-ens.) ---
        if (hn(n).gt.0) then
          hwc(n)  = hwc(n)  / hn(n)   
          hqlc(n) = hqlc(n) / hn(n) 
        endif  
        
        !--- sub-ensemble-contributions (div. by binsize) ---
        if (nfield.gt.0) then
          hn(n)   = hn(n)   / (nfield*sizmin)
          hac(n)  = hac(n)  / (nfield*sizmin)
          hmc(n)  = hmc(n)  / (nfield*sizmin)
        endif


        if (hnlev(n).gt.0) then
          
          !--- sub-ensemble averages at z----
          hthl(n) = hthl(n) / hnlev(n)
          hthv(n) = hthv(n) / hnlev(n)
          hql(n)  = hql(n)  / hnlev(n)
          hqt(n)  = hqt(n)  / hnlev(n)
          hw(n)   = hw(n)   / hnlev(n)
          hm(n)   = hm(n)   / hnlev(n)
          
          !--- nr of clouds in sub-ensemble at z----
          hnlev(n) = hnlev(n)      ! What is this for?
          
        endif
        
        if (hn(n).ne.0.) then  
          ndum = ndum + hn(n)*sizmin
          adum = adum + hac(n)*sizmin
          !mdum = mdum + hmc(n)*sizmin
        endif

      enddo

     
      print *,'------------------------------------------ '
      print *,'----------time averaging ----------------- '
      print *,'------------------------------------------ '
      print *,'--- itim:',itim, '  averaging period:',(itim-1)*dtav,'-',(itim)*dtav
      print *,'--- total nr. of fields read from file:',nfield
      print *,'--- total nr. of selected clouds:',ndum*nfield
      print *,'---   selected clouds per field:',ndum
      print *,'---   their proj. cloud fraction:',adum*100,'%'
      !print *,'---   their mass flux:',mdum,'m/s'
      print *,'------------------------------------------ '
      print *,' '
 
     
      !--- request nr of output bins ----
      lcloudmax=0.
      nlcloudmax = 1
      do n=1,nbin
        if (hn(n).gt.0) then
          lcloudmax=hsize(n)
          nlcloudmax=n
        endif
      enddo
      
      
      !--- NetCDF output ----
      fname_out = "cusize_output.nc"
      print *,'writing size histograms to NetCDF file ',fname_out

      if (FIRST_CALL) then

        print *,'   --> creating file '
        STATUS = nf_create(fname_out, nf_clobber, ncidout)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        ! global attributes
        status = nf_put_att_text  (ncidout, nf_global, 'title', 24, 'ICON cloud size analysis')
        STATUS = NF_PUT_ATT_real  (NCIDout, nf_global, 'dtav', nf_real, 1, dtav)
        STATUS = NF_PUT_ATT_double(NCIDout, nf_global, 'jdate0', nf_double, 1, jdate0)
        status = nf_put_att_text  (ncidout, nf_global, 'size_definition', 30, sizname)
        STATUS = NF_PUT_ATT_int   (NCIDout, nf_global, 'profile_overlap_mode', nf_int, 1, ioverlap)
        STATUS = NF_PUT_ATT_int   (NCIDout, nf_global, 'cloud_selection_mode', nf_int, 1, iselect)

        ! dimensions
        status = NF_DEF_DIM(ncidout,'time'    ,nf_unlimited,dimidt)  !make time dimension unlimited
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_DIM(ncidout,'size'    ,nbin,dimidl)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_DIM(ncidout,'height'  ,kmax,dimidk)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        ! coordinate variables
        dd(1)=dimidt
        status=nf_def_var(ncidout,'time'  ,nf_real,1,dd,varidt)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        dd(1)=dimidl
        status=nf_def_var(ncidout,'size'  ,nf_real,1,dd,varidl)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        dd(1)=dimidk
        status=nf_def_var(ncidout,'height',nf_real,1,dd,varidk)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        status = nf_put_att_text(ncidout, varidt, 'long_name', 25, 'Time since initialization')
        status = nf_put_att_text(ncidout, varidt, 'units', 1, 's')
        status = nf_put_att_text(ncidout, varidl, 'long_name', 4, 'Size')
        status = nf_put_att_text(ncidout, varidl, 'units', 1, 'm')
        status = nf_put_att_text(ncidout, varidk, 'long_name', 6, 'Height')
        status = nf_put_att_text(ncidout, varidk, 'units', 1, 'm')

        ! scalar statistics
        dd(1)=dimidt
        status=nf_def_var(ncidout,'jdate'  ,nf_double,1,dd,varidjd)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'seconds',nf_double,1,dd,varidsec)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'nfield' ,nf_int  ,1,dd,varidfld)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'nclouds',nf_float,1,dd,varidncld)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'atot'   ,nf_float,1,dd,varidatot)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        !status=nf_def_var(ncidout,'mtot'   ,nf_float,1,dd,varidmtot)
        !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'lmax'   ,nf_float,1,dd,varidlmax)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status=nf_def_var(ncidout,'aptot'  ,nf_float,1,dd,varidaptot)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = nf_put_att_text(ncidout, varidjd  , 'long_name', 11, 'Julian date')
        status = nf_put_att_text(ncidout, varidjd  , 'units', 16, '<daynr>.fraction')
        status = nf_put_att_text(ncidout, varidsec , 'long_name', 28, 'Seconds since initialization')
        status = nf_put_att_text(ncidout, varidsec , 'units', 1, 's')
        status = nf_put_att_text(ncidout, varidfld , 'long_name', 23, 'Nr of scenes in average')
        status = nf_put_att_text(ncidout, varidfld , 'units', 1, '#')
        status = nf_put_att_text(ncidout, varidncld, 'long_name', 20, 'Nr of sampled clouds')
        status = nf_put_att_text(ncidout, varidncld, 'units', 1, '#')
        status = nf_put_att_text(ncidout, varidatot, 'long_name', 33, 'Summated projected cloud fraction')
        status = nf_put_att_text(ncidout, varidatot, 'units', 1, '0-1')
        !status = nf_put_att_text(ncidout, varidmtot, 'long_name', 20, 'Cumulative mass flux')
        !status = nf_put_att_text(ncidout, varidmtot, 'units', 3, 'm/s')
        status = nf_put_att_text(ncidout, varidlmax, 'long_name', 18, 'Maximum cloud size')
        status = nf_put_att_text(ncidout, varidlmax, 'units', 1, 'm')
        status = nf_put_att_text(ncidout, varidaptot, 'long_name', 24, 'Projected cloud fraction')
        status = nf_put_att_text(ncidout, varidaptot, 'units', 3, '0-1')

        ! profile statistics
        dd(1)=dimidk
        dd(2)=dimidt
        status = NF_DEF_VAR(ncidout,'a_prof',nf_real,2,dd,varidap)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'ql_prof',nf_real,2,dd,varidqlp)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = nf_put_att_text(ncidout, varidap, 'long_name', 24, 'Slab-mean cloud fraction')
        status = nf_put_att_text(ncidout, varidap, 'units', 3, '0-1')
        status = nf_put_att_text(ncidout, varidqlp, 'long_name', 26, 'Slab-mean cloud condensate')
        status = nf_put_att_text(ncidout, varidqlp, 'units', 5, 'kg/kg')

        ! size densities
        dd(1)=dimidl
        dd(2)=dimidt
        status = NF_DEF_VAR(ncidout,'hn',nf_real,2,dd,varidn)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hac',nf_real,2,dd,varidac)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hwc',nf_real,2,dd,varidwc)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hmc',nf_real,2,dd,varidmc)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hqlc',nf_real,2,dd,varidqlc)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hnns_mean',nf_real,2,dd,varidmean)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hnns_stdev',nf_real,2,dd,varidstdev)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'D0',nf_real,2,dd,varidd0)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'D1',nf_real,2,dd,varidd1)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'ncloud_bin',nf_real,2,dd,varidncb)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        status = nf_put_att_text(ncidout, varidn  , 'long_name', 14, 'Number density')
        status = nf_put_att_text(ncidout, varidn  , 'units', 3, '1/m')
        status = nf_put_att_text(ncidout, varidac , 'long_name', 12, 'Area density')
        status = nf_put_att_text(ncidout, varidac , 'units', 3, '1/m')
        status = nf_put_att_text(ncidout, varidwc , 'long_name', 25, 'Average vertical velocity')
        status = nf_put_att_text(ncidout, varidwc , 'units', 3, 'm/s')
        status = nf_put_att_text(ncidout, varidmc , 'long_name', 17, 'Mass flux density')
        status = nf_put_att_text(ncidout, varidmc , 'units', 3, '1/s')
        status = nf_put_att_text(ncidout, varidqlc, 'long_name', 24, 'Average cloud condensate')
        status = nf_put_att_text(ncidout, varidqlc, 'units', 5, 'kg/kg')
        status = nf_put_att_text(ncidout, varidmean, 'long_name', 45, 'Meanshortest distance between clouds in bin')
        status = nf_put_att_text(ncidout, varidmean, 'units', 1, 'm')
        status = nf_put_att_text(ncidout, varidstdev, 'long_name', 60, 'Standard deviation of shortest distance between clouds in bin')
        status = nf_put_att_text(ncidout, varidstdev, 'units', 1, 'm')
        status = nf_put_att_text(ncidout, varidd0, 'long_name', 42, 'Geometrical mean distance of clouds in bin')
        status = nf_put_att_text(ncidout, varidd0, 'units', 1, 'm')
        status = nf_put_att_text(ncidout, varidd1, 'long_name', 42, 'Arithmetical mean distance of clouds in bin')
        status = nf_put_att_text(ncidout, varidd1, 'units', 1, 'm')

        ! size density profiles
        dd(1)=dimidl
        dd(2)=dimidk
        dd(3)=dimidt
        status = NF_DEF_VAR(ncidout,'hnlev_prof',nf_real,3,dd,varidpnlev)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hw_prof',nf_real,3,dd,varidpw)  
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hthl_prof',nf_real,3,dd,varidpthl)  
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hthv_prof',nf_real,3,dd,varidpthv)  
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hqt_prof',nf_real,3,dd,varidpqt)  
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        status = NF_DEF_VAR(ncidout,'hql_prof',nf_real,3,dd,varidpql)  
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        status = nf_put_att_text(ncidout, varidpnlev, 'long_name', 24, 'Number of sampled clouds')
        status = nf_put_att_text(ncidout, varidpnlev, 'units', 1, '#')
        status = nf_put_att_text(ncidout, varidpw   , 'long_name', 25, 'Average vertical velocity')
        status = nf_put_att_text(ncidout, varidpw   , 'units', 3, 'm/s')
        status = nf_put_att_text(ncidout, varidpthl , 'long_name', 15, 'Average theta_l')
        status = nf_put_att_text(ncidout, varidpthl , 'units', 1, 'K')
        status = nf_put_att_text(ncidout, varidpthv , 'long_name', 15, 'Average theta_v')
        status = nf_put_att_text(ncidout, varidpthv , 'units', 1, 'K')
        status = nf_put_att_text(ncidout, varidpqt  , 'long_name', 11, 'Average q_t')
        status = nf_put_att_text(ncidout, varidpqt  , 'units', 5, 'kg/kg')
        status = nf_put_att_text(ncidout, varidpql  , 'long_name', 11, 'Average q_l')
        status = nf_put_att_text(ncidout, varidpql  , 'units', 5, 'kg/kg')

        ! Change mode of netCDF operation
        status = NF_ENDDEF(ncidout)
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        ! fill coordinate variables
        start1(1)=1
        count1(1)=nbin
        status = NF_PUT_VARA_REAL(ncidout,varidl,start1,count1,hsize(1:nbin) )
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)
        start1(1)=1  
        count1(1)=kmax
        status = NF_PUT_VARA_REAL(ncidout,varidk,start1,count1,z_f(1,1:kmax) )
        if (STATUS .ne. nf_noerr) call handle_err(STATUS)

        FIRST_CALL = .FALSE.

      ENDIF

      print *,'   --> adding values '

      ! add time value
      !start1(1)=itim
      start1(1)=1
      count1(1)=1


      status = NF_PUT_VARA_REAL(ncidout,varidt,start1,count1, (itim-0.5)*dtav )
      !status = NF_PUT_VARA_REAL(ncidout,varidt,start1,count1, 1*dtav )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)


      ! add scalar values
      status = NF_PUT_VARA_DOUBLE(ncidout,varidjd ,start1,count1, jdate )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_DOUBLE(ncidout,varidsec,start1,count1, seconds )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_INT (ncidout,varidfld  ,start1,count1, nfield )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidncld ,start1,count1, ndum*nfield )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidatot ,start1,count1, adum )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidmtot ,start1,count1, mdum )
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidlmax ,start1,count1, lcloudmax )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidaptot,start1,count1, aptot )
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      ! add profiles
      start2(1)=1
      !start2(2)=itim
      start2(2)=1
      count2(1)=kmax
      count2(2)=1

      
      status = NF_PUT_VARA_REAL(ncidout,varidap,start2,count2,pa)
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidqlp,start2,count2,pql)
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      ! add 1D size densities
      start2(1)=1
      !start2(2)=itim
      start2(2)=1
      count2(1)=nbin
      count2(2)=1


      status = NF_PUT_VARA_REAL(ncidout,varidn,start2,count2,hn(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidac,start2,count2,hac(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidwc,start2,count2,hwc(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidmc,start2,count2,hmc(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidqlc,start2,count2,hqlc(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      status = NF_PUT_VARA_REAL(ncidout,varidmean,start2,count2,hnns_mean(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidstdev,start2,count2,hnns_stdev(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidd0,start2,count2,D0(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidd1,start2,count2,D1(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      status = NF_PUT_VARA_REAL(ncidout,varidncb,start2,count2,ncloud_bin(1:nbin))
      if (STATUS .ne. nf_noerr) call handle_err(STATUS)

      ! add 2D size densities
      start3(1)=1
      start3(2)=1
      !start3(3)=itim
      start3(3)=1
      count3(1)=nbin
      !count3(2)=kmax
      !count3(3)=1
      count3(2)=1
      !status = NF_PUT_VARA_REAL(ncidout,varidpnlev,start3,count3,hnlev(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidpw,start3,count3,hw(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidpthl,start3,count3,hthl(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidpthv,start3,count3,hthv(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidpqt,start3,count3,hqt(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)
      !status = NF_PUT_VARA_REAL(ncidout,varidpql,start3,count3,hql(1:nbin))
      !if (STATUS .ne. nf_noerr) call handle_err(STATUS)



      print *,' '
     
     
      endif  !nfield>0

      
      return
      end
      
      
      
!
!-----------------------------------------------------------------
!								  
      subroutine nextfield(ff)
      
      implicit none

      include "netcdf.inc"
      include 'parconst'
      include 'comvar'
      include 'comnc'

      integer STATUS, VARID, START1, COUNT1, START3(3), COUNT3(3)

      integer i,ff,v
      integer dql,dw,dqt,dthl,dthv 

      real temp0,temp3(cmax,kmax),dist


      !--- reset ---
      dat (:)   = 0
      mark(:)   = -1
      qt  (:,:) = 0
      thv (:,:) = 0
      thl (:,:) = 0
      ql  (:,:) = 0
      w   (:,:) = 0

      
      !--- read 3D fields from nc files ---
      start3 = (/ 1   , 1   , ff  /)
    
      !  cloud liquid water
      count3 = (/ cmax, kmax, 1   /)
      status = NF_INQ_VARID     (ncidql, 'qc', varid)    !qc(time, height, ncells)   specific_cloud_water_content    "kg kg-1"
      if (status /= nf_noerr) call handle_err(status)
      status = NF_GET_VARA_REAL (ncidql, varid, start3, count3, temp3)
      if (status /= nf_noerr) then
        call handle_err(status)
      else
        ql = temp3
      endif
      

      !--- identify cloudy gridboxes & set cloud id (mask) to 0 ---
      do i=1,cmax

        if (SUM(ql(i,:)).gt.1.E-8) then     !criterion: condensate>0 kg/kg
          dat (i) = 1
          mark(i) = 0
        end if

      enddo


      !--- increase nr of fields for averaging ---
      nfield = nfield + 1
      
      
      return
      end
      
      
      
!
!-----------------------------------------------------------------
!								  
      subroutine slabstat
      
      include 'parconst'
      include 'comvar'

      integer i
      real ap(cmax)
      
      print *,'calculating slab-average profiles ...'
      
      pql=0
      pa=0
      
      do i=1,cmax
        ap(i) = 0
      enddo
      
      do i=1,cmax
        pql = pql + SUM(ql(i,:))
        pa  = pa  + dat(i)
        if (dat(i).eq.1) then
          ap(i) = 1
        endif
      enddo
      
      aptot = 0.
      do i=1,cmax
        aptot = aptot + ap(i)
      enddo
      aptot = aptot / (cmax)
      
      pql = pql / (cmax)
      pa  = pa  / (cmax)
      
      !print *,'   maximum cloud fraction = ',maxval(pa)
      print *,'   projected cloud fraction = ',aptot*100.,'%'
      
      return
      end


!
!-----------------------------------------------------------------
!								  
      subroutine cluster(indcentre3)
      
      include 'parconst'
      include 'comvar'

      integer i,n
      real indcentre3(cmax,4)
 
      print *,'clustering ...'
      
      ncloud = 0
      do i=1,cmax
      
        if (dat(i).gt.0.and.mark(i).eq.0) then
          ncloud = ncloud + 1
          mark(i) = ncloud
          call constructcloud(i,indcentre3)
        endif   
        
      enddo
      
      print *,'  nr of clouds in domain =',ncloud
      print *,'  '
      
      return
      end
   
                
      
!
!------ find all points connected to cloudy point i1,k1 ---------
!								  

      subroutine constructcloud(i1,indcentre2)
      
      include 'parconst'
      include 'comvar'

      integer i,n,m,i1,i2,k2,pos,hbot,htop
      integer nr2,nrtemp,start3(2),count3(2)
      real indcentre2(cmax,4)

      start3 = (/ 1   , 1    /)
      count3 = (/ cmax, vmax /)
 
      do n=1,cldmax      !maximum possible amount of clouds (cmax*kmax)
        ind(n,2)=-1      !checker status. -1: untreated. 0: part of cloud, but unchecked for neighbors. Other: part of cloud and checked
      enddo
      
      nr2=1
      ind(nr2,1)=i1
      ind(nr2,2)=0
      
      m=mark(i1)
      print *,'    cloud found, label ',m,i1

      !print *,'    searching neighbouring points...'
      nrtemp=0
      do while(nr2.gt.nrtemp)
      
       nrtemp = nr2    !-- new check --
       do n=1,nrtemp
        
        if (ind(n,2).eq.0) then

          !--belongs to cloud, but not yet checked for neighbours ---
          
          ind(n,2)=m 
          i2=ind(n,1)
          
          !check horizontally adjacent cells 
          do i=1,vmax
            pos = cell_neighbor(ind(n,1),i)
           
            !print *, '        checking:',ind(n,1),i,cell_neighbor(ind(n,1),i)
            if (pos.gt.0) then
              call neighbour(pos,m,nr2)
            endif
          enddo
       
   
          !check vertically adjacent cells
          !call mod(k2+1,1,kmax,pos)
          !call neighbour(i2,pos,m,nr2)

          !call mod(k2-1,1,kmax,pos)
          !call neighbour(i2,pos,m,nr2)
          
        endif  
          
       enddo  
       
      enddo  !-- no new neighbours --

      print *,'      nr of cloudpoints:',nr2

      call cloudsizes(nr2,m,indcentre2)

      
      select case (iselect)
      case (1)
        call addtohist(nr2,m,indcentre2)
      case (2)
        if (kbot(m).le.kselect) then
          call addtohist(nr2,m,indcentre2)
        endif  
      case (3)
        if (kbot(m).le.kselect.and.ktop(m).gt.kselect) then
!          print *,kbot(m),kselect,ktop(m)
          call addtohist(nr2,m,indcentre2)
        endif  
      end select  
      
            
      return
      end
        


!
!------ calc cloud sizes ---------
!	
							  
      subroutine cloudsizes(nrr,mm,indcentre)
      
      include 'parconst'
      include 'comvar'

      integer ii,n,mm,prof(kmax),surf(cmax),nrr,v
      real max_lat,max_lon,min_lat,min_lon
      real indcentre(cmax,4)
      

!      !--- cloud size: height ---  
!      do kk=1,kmax 
!        prof(kk)=0
!      enddo
!      
!      do n=1,nrr
!        prof(ind(n,3))=1
!      enddo
!      
!      n=0
!      ii=0
!      do kk=1,kmax 
!        n=n+prof(kk)
!        if (prof(kk).eq.1.and.ii.eq.0) then
!          kbot(mm)=kk
!          ii=1
!        endif
!        if (prof(kk).eq.0.and.ii.eq.1) then
!          ktop(mm)=kk-1
!          ii=0
!        endif
!      enddo
!      
!      l(1,mm) = n*dz
!      
!      
!      
!      !--- cloud size: sqrt(height-averaged area) ---  
!      l(2,mm) = ((nrr*darea)/n)**0.5
      
      
      
      !--- cloud size: sqrt(vert.projected area) ---  
      do ii=1,cmax
        surf(ii)=0
      enddo
      
      do n=1,nrr
        surf(ind(n,1))=1
      enddo
      
      n=0
      do ii=1,cmax
        n=n+surf(ii)
      enddo
         
 
      l(3,mm) = (n*darea)**0.5

      print*, 'projected cloud size:',l(3,mm)

      ! ---- coordinates centre points-----------

      max_lon = -1.
      max_lat = 0.
      min_lon = 1.
      min_lat = 1.

      do n=1,nrr ! this is wrong! (but why again?)
        do v=1,3
          if (cell_vx(v,ind(n,1)).gt.max_lon) then
            max_lon = cell_vx(v,ind(n,1))
          endif
          if (cell_vy(v,ind(n,1)).gt.max_lat) then      
            max_lat = cell_vy(v,ind(n,1))
          endif
          if (cell_vx(v,ind(n,1)).lt.min_lon) then
            min_lon = cell_vx(v,ind(n,1))
          endif
          if (cell_vy(v,ind(n,1)).lt.min_lat) then
            min_lat = cell_vy(v,ind(n,1))
          endif
        enddo
      enddo
     
      indcentre(mm,1) = (max_lon+min_lon)/2
      indcentre(mm,2) = (max_lat+min_lat)/2

      !--- cloud size: (volume)**(1./3.) ---  
      l(4,mm) = (nrr*darea*dz)**(1./3.)    
      


      !print *, '       sizes:',(l(n,mm),n=1,4)

      
      return
      end
        
      
!
!-----------------------------------------------------------------
!								  
      subroutine addtohist(nrr,mm,indcentre)
      
      include 'parconst'
      include 'comvar'

      integer nn,binnr,ii,jj,kb,nrr,mm,nrcldlev
      
      real dthl,dqt,dw,dql,dthv,da,dm
      real dcw,dcql,dca,dcm,indcentre(cmax,4)
              
        !--- sort the cloud by its size ---
        binnr = floor(l(isize,mm)/sizmin)
        indcentre(mm,4) = int(binnr)       
        
        dthl=0.
        dthv=0.
        dqt=0.
        dql=0.
        dw=0.
        da=0.
        
        dcql=0.
        dcw=0.
        dcm=0.
        dca=0.

      
!        !--- compute hor. cloud-averages ---
!        do nn=1,nrr
!          ii=ind(nn,1)
!          kk=ind(nn,3)
!          
!          
!          !--- cloud overlay options ---
!          select case (ioverlap)
!          case (1)
!            kb = kk - kbot(mm) + 1 
!          case (2)
!            kb = kmax + kk - ktop(mm)
!          case (3)
!            kb = kk
!          end select
!          
!          da(kb)   = da(kb) + 1.
!          dthl(kb) = dthl(kb) + thl(ii,kk)
!          dthv(kb) = dthv(kb) + thv(ii,kk)
!          dqt(kb)  = dqt(kb) + qt(ii,kk)
!          dql(kb)  = dql(kb) + ql(ii,kk)
!          dw(kb)   = dw(kb) + w(ii,kk)
!
!        enddo
!        
!
!        do kk=1,kmax 
!        
!        if (da(kk).gt.0) then
!          
!          nrcldlev = da(kk)
!          
!          !--- fractional stuff ---
!          da(kk) = da(kk) / (cmax)
!          dm(kk) = dw(kk) / (cmax)
!        
!          !--- add to sub-ensemble ---
!          dcql = dcql + dql(kk)
!          dcw  = dcw + dw(kk)
!          dcm  = dcm + dm(kk)
!          
!          !--- level-averages ---
!          dthl(kk) = dthl(kk) / nrcldlev 
!          dthv(kk) = dthv(kk) / nrcldlev 
!          dqt(kk)  = dqt(kk) / nrcldlev 
!          dql(kk)  = dql(kk) / nrcldlev 
!          dw(kk)   = dw(kk)  / nrcldlev 
!          
!!          print *,kk,nrcldlev,dw(kk) 
!
!          !--- add to sub-ensemble-profile stack ---
!          hnlev(binnr,kk) = hnlev(binnr,kk) + 1.
!          ha(binnr,kk)    = ha(binnr,kk) + da(kk)
!          hthl(binnr,kk)  = hthl(binnr,kk) + dthl(kk)
!          hthv(binnr,kk)  = hthv(binnr,kk) + dthv(kk)
!          hql(binnr,kk)   = hql(binnr,kk) + dql(kk)
!          hqt(binnr,kk)   = hqt(binnr,kk) + dqt(kk)
!          hw(binnr,kk)    = hw(binnr,kk) + dw(kk)
!          hm(binnr,kk)    = hm(binnr,kk) + dm(kk)
!                  
!        endif
!        
!        enddo
        
      !--- normalize sub-ensemble ---
      dcql = dcql / nrr          !volume average
      dcw  = dcw  / nrr          !volume average
      dcm  = dcm  / (l(1,mm)/dz) !height-averaged mass flux
      dca  = ( l(3,mm)**2. )/(darea*cmax)  !vertically projected area


      !--- add to sub-ensemble stack ---
      hn(binnr)   = hn(binnr) + 1
      hac(binnr)  = hac(binnr) + dca
      hwc(binnr)  = hwc(binnr) + dcw
      hmc(binnr)  = hmc(binnr) + dcm
      hqlc(binnr) = hqlc(binnr) + dcql
     
      return
      end
      
      
      
!
!------ label the neighbours ---------------
!								  
        
      subroutine neighbour(i3,mm,nrr)
      
      include 'parconst'
      include 'comvar'

      integer i3,mm,nrr
     
      if ( dat(i3).gt.0. .and. mark(i3).eq.0 ) then
        !print *, '        found:',i3
        mark(i3)=mm
        nrr=nrr+1

        ind(nrr,1)=i3
        ind(nrr,2)=0
      endif
 
      return
      end


      !-----calculate nearest-neighbour distances-------------
      
      subroutine nnspacing(indcentre)
 
      include 'parconst'
      include 'comvar'

      integer bb,cc,aa,ff,ss,pairs
      real indcentre(cmax,4),firstx,firsty,secondx,secondy,mean,std
      real dist1,dist(ncloud*(ncloud-1)/2),binclouds(ncloud,2)
      real distnew,distold,disttotal,dist_all(ncloud*(ncloud-1)/2)
      real D1_sum
      double precision D0_product
      
      Nmax = cmax/2

      do bb=1,nbin
        binclouds = 0
        cc = 0
        D0_product = 1. ! D0 from Tobin et al (2012), geometrical mean
        D1_sum = 0. ! D1 from Tobin et al (2012), arithmetical mean
        dist_all(:) = 0.
        do aa = 1,ncloud !total nr of clouds
          if (indcentre(aa,4)==bb) then
            cc = cc+1
            binclouds(cc,1) = indcentre(aa,1)
            binclouds(cc,2) = indcentre(aa,2)
          endif
        enddo
        ! calculate the nns between the x- and y-coordinates (lon and lat)
        pairs = 0
        dist = 0
        do ff = 1,cc
          distold = 400000.
          disttotal = 0.
          firstx = binclouds(ff,1)
          firsty = binclouds(ff,2)
          do ss = 2,cc
            if (ss.gt.ff) then
              secondx = binclouds(ss,1)
              secondy = binclouds(ss,2)
              distnew = haversine(firsty,firstx,secondy,secondx)
              pairs   = pairs+1
              D0_product = D0_product * distnew * 1E-5
        print *, 'D0_product:',D0_product
              D1_sum  = D1_sum + distnew
              dist_all(pairs) = distnew
              if (distnew.le.distold) then
                distold = distnew
              endif
            endif
            enddo
          enddo 
          disttotal       = disttotal+distnew
          dist(cc)        = distnew
          ncloud_bin(bb)  = cc
          bin_density(bb) = real(cc)/Nmax
          mean            = disttotal/real(cc)
          !mean_all        = SUM(dist_all)/real(pairs)
          std             = sqrt(SUM((dist(1:cc)-mean)**2)/(real(cc)-1))
          !std_all         = sqrt(SUM((dist_all(1:pairs)-mean_all)**2)/(real(pairs)-1))
          if (cc.ge.2) then
            !mean = SUM(dist/real(pairs))
            hnns_mean(bb) = mean  
          else
            hnns_mean(bb) = 0.
          endif
          if (cc.ge.3) then
            !std = sqrt(SUM((dist(1:pairs)-mean)**2)/(real(pairs)-1))
            hnns_stdev(bb) = std  
          else
            hnns_stdev(bb) = 0.
          endif 
          if (pairs.ge.1) then
            D0(bb) = D0_product**(1/real(pairs)) *1E5
            D1(bb) = D1_sum*(1/real(pairs))
          else
            D0(bb) = 0.
            D1(bb) = 1.
          endif

        !enddo
      enddo

        print *, 'ncloud_bin:',ncloud_bin(1:6)
        print *, 'hnns_mean:',hnns_mean(1:6)
        print *, 'hnns_stdev:',hnns_stdev(1:6)
        print *, 'bin_density:',bin_density(1:6)
        print *, 'D0:',D0(1:6)
        print *, 'D1:',D1(1:6)

      return
      end


      !-----haversine formula to calculate distance between two points----
      
      function haversine(deglat1,deglon1,deglat2,deglon2) result (distance)
        
      real,intent(in) :: deglat1,deglon1,deglat2,deglon2
      real :: a,c,distance,dlat,dlon,lat1,lat2
      real,parameter :: radius = 6371000

      dlat = deglat2-deglat1
      dlon = deglon2-deglon1
      lat1 = deglat1
      lat2 = deglat2
      a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
      c = 2*asin(sqrt(a))
      distance = radius*c
     
      end function haversine

!
!------ check boundaries ---------------
!								  
        
      subroutine mod(x,lx,ux,ps)
      
      integer x,lx,ux,ps
      
      ps = x  
      if (x.lt.lx) then
        ps = ux
        !ps = lx
      else if(x.gt.ux) then
        ps = lx
        !ps = ux
      endif
      
      return
      end



!
!   Conversion from a Gregorian calendar date to a Julian date. Valid for any Gregorian calendar date 
!   producing a Julian date greater than zero:
!
!   http://aa.usno.navy.mil/faq/docs/JD_Formula.php
!
      INTEGER FUNCTION JD (YEAR,MONTH,DAY)

!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY).

      INTEGER YEAR,MONTH,DAY,I,J,K

      I= YEAR
      J= MONTH
      K= DAY

      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)  /12-3*((I+4900+(J-14)/12)/100)/4

      RETURN
      END



!
!   Conversion from a Julian date to a Gregorian calendar date.
!
!   http://aa.usno.navy.mil/faq/docs/JD_Formula.php
!
      SUBROUTINE GDATE (JD, YEAR,MONTH,DAY)

!---COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!   GIVEN THE JULIAN DATE (JD).

      INTEGER JD,YEAR,MONTH,DAY,I,J,K
  
      L= JD+68569
      N= 4*L/146097
      L= L-(146097*N+3)/4
      I= 4000*(L+1)/1461001
      L= L-1461*I/4+31
      J= 80*L/2447
      K= L-2447*J/80
      L= J/11
      J= J+2-12*L
      I= 100*(N-49)+I+L
  
      YEAR= I
      MONTH= J
      DAY= K
  
      RETURN
      END


      subroutine handle_err(errcode)

      !use netcdf

      implicit none

      include "netcdf.inc"
      
      integer errcode
     
      write(6,*) 'Error: ', nf_strerror(errcode)
      !write(6,*) 'Error: ', nf90_strerror(errcode)

      !stop 2
      
      end subroutine handle_err
           

       !ELEMENTAL FUNCTION vector_product (x0, x1) result(x2)
       !REAL FUNCTION vector_product (xx0, xx1) result(xx2)
    
       subroutine vector_product(xx0,xx1,xx2)

        ! TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1
    
        ! TYPE(t_cartesian_coordinates) :: x2

        real xx0(3),xx1(3),xx2(3)        
    
         !------------------------------------------------------
         xx2(1) = xx0(2)*xx1(3) - xx0(3)*xx1(2)
         xx2(2) = xx0(3)*xx1(1) - xx0(1)*xx1(3)
         xx2(3) = xx0(1)*xx1(2) - xx0(2)*xx1(1)

         !x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
         !x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
         !x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)
    
       !END FUNCTION vector_product

       end

