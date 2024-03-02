  implicit real *8 (a-h,o-z)
  character *100, fstr
  
  fstr = 'a380.msh'
  call get_filetype(fstr, ifiletype, ier)
  print *, ".msh file"
  print *, "ifiletype=,", ifiletype, "ier=", ier

  fstr = 'prism_50.gidmsh'
  call get_filetype(fstr, ifiletype, ier)
  print *, ".gidmsh file"
  print *, "ifiletype=,", ifiletype, "ier=", ier

  fstr = 'cuboid_a1_b2_c1p3.tri'
  call get_filetype(fstr, ifiletype, ier)
  print *, ".tri file"
  print *, "ifiletype=,", ifiletype, "ier=", ier

  fstr = 'lens_r00.msh'
  call get_filetype(fstr, ifiletype, ier)
  print *, "gmsh file"
  print *, "ifiletype=,", ifiletype, "ier=", ier

  stop
  end





subroutine get_filetype(fstr, ifiletype, ier)
!
!  This subroutine determines the mesh type of a given
!  file. Currently supported formats include
!
!  * .msh
!  * .tri
!  * .gidmsh
!  * .msh (from gmsh v2 or v4)
!
!  Input arguments:
!    - fstr: string
!         file name
!  Output arguments
!    - ifiletype: integer
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .gmsh
!    - ier: integer
!        Error code
!        * ier = 0, successful execution
!        * ier = 8, file format not recognized
!
!
  implicit real *8 (a-h,o-z)
  integer i1,i2, io
  character (len=*) fstr
  character *100, fstr_gid_tritest, fstr_gmshtest

  ier = 0
  fstr_gid_tritest = 'MESH dimension 3 ElemType Triangle Nnode '
  fstr_gmshtest = '$MeshFormat'

! 
  open(unit=33, file=trim(fstr), status='old')

!  Check if it is a .msh file
  io = 0   
  read(33,*, iostat = io) i1,i2, i3, i4, i5
  if(io.eq.0) then
    ifiletype = 1
    return
  endif
      
!   Check if it is a .tri file now      
  rewind(33)
  read(33,*, iostat = io) i1,i2
  if(io.eq.0) then
    ifiletype = 2
    return
  endif

!
!  Now check if it is a gidmsh
!
  rewind(33)
  read(33,'(a)', iostat = io) fstr
      
  len1 = len(trim(fstr_gid_tritest))
  len2 = len(trim(fstr))
  if(trim(fstr(1:len1)) .eq. trim(fstr_gid_tritest)) then
     read(fstr((len1+1):len2),*,iostat=io) iind
     ifiletype = 3
     return
  endif

  len1 = len(trim(fstr_gmshtest))
  if(trim(fstr(1:len1)) .eq. trim(fstr_gmshtest)) then
     read(33,*,iostat=io) tmp, tmp2, tmp3
     if(abs(tmp-2).le.1) then
       ifiletype = 4
     elseif(abs(tmp-4).le.1) then
       ifiletype = 5
     else
       print *, "Invalid gmsh file type"
       ier = 8
     endif
     return
  endif

  print *, "Invalid file type"
  ier = 8

  close(33)

  return
  end
