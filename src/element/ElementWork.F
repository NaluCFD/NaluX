      subroutine hex_scs_det( nelem, npe, nscs,
     .  cordel, area_vec )
c
c***********************************************************************
c***********************************************************************
c
c description:
c     This  routine returns the area vectors for each of the 12
c     subcontrol surfaces of the hex.
c
c formal parameters - input:
c     nelem         int   number of elements in the workset
c     npe           int   number of nodes per element
c     nscs          int   number of sub control surfaces
c     cordel        real  element local coordinates
c
c formal parameters - output:
c     area_vec      real  area vector: A = Ax i + Ay j + Az k
c                                        = Area*(i + j + k)
c
c***********************************************************************
c
      implicit none
      integer nelem, npe, nscs
      double precision cordel, area_vec

      dimension cordel(3, npe, nelem)
      dimension area_vec(3, nelem, nscs)

      double precision coords, scscoords

      integer ielem, j, k, ics, inode, itrianglenode, npf
      ! coordinates of the vertices needed to define the scs
      dimension coords(3,27)
      ! coordinates of the vertices of the scs
      dimension scscoords(3,4)

      double precision half, one4th, one8th

      ! this table defines the scs
      integer HexEdgeFacetTable(4,12)
      data HexEdgeFacetTable /
     .  21, 9, 13, 27,   ! sc face 1 -- points from 1 -> 2
     .  25, 10, 13, 27,  ! sc face 2 -- points from 2 -> 3
     .  11, 13, 27, 24,  ! sc face 3 -- points from 3 -> 4
     .  12, 26, 27, 13,  ! sc face 4 -- points from 1 -> 4
     .  14, 21, 27, 18,  ! sc face 5 -- points from 5 -> 6
     .  18, 15, 25, 27,  ! sc face 6 -- points from 6 -> 7
     .  18, 16, 24, 27,  ! sc face 7 -- points from 7 -> 8
     .  17, 18, 27, 26,  ! sc face 8 -- points from 5 -> 8
     .  20, 21, 27, 26,  ! sc face 9 -- points from 1 -> 5
     .  21, 19, 25, 27,  ! sc face 10 -- points from 2 -> 6
     .  23, 24, 27, 25,  ! sc face 11 -- points from 3 -> 7
     .  22, 26, 27, 24 / ! sc face 12 -- points from 4 -> 8

      half = 0.5d0
      one4th = 0.25d0
      one8th = 0.125d0

      ! loop over elements
      do ielem = 1, nelem
        ! element vertices
        do j = 1,8
          do k = 1,3
            coords(k,j) = cordel(k,j,ielem)
          end do
        end do

        ! face 1 (front)
        ! 4++++11+++3
        ! +         +
        ! +         +
        ! 12   13   10
        ! +         +
        ! +         +
        ! 1++++9++++2

        ! edge midpoints
        do k = 1,3
          coords(k,9) = half*(cordel(k,1,ielem) + cordel(k,2,ielem))
        end do
        do k = 1,3
          coords(k,10) = half*(cordel(k,2,ielem) + cordel(k,3,ielem))
        end do
        do k = 1,3
          coords(k,11) = half*(cordel(k,3,ielem) + cordel(k,4,ielem))
        end do
        do k = 1,3
          coords(k,12) = half*(cordel(k,4,ielem) + cordel(k,1,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,13) = one4th*(cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,3,ielem) + cordel(k,4,ielem))
        end do

        ! face 2 (back)
        ! 8++++16+++7
        ! +         +
        ! +         +
        ! 17   18   15
        ! +         +
        ! +         +
        ! 5++++14+++6

        ! edge midpoints
        do k = 1,3
          coords(k,14) = half*(cordel(k,5,ielem) + cordel(k,6,ielem))
        end do
        do k = 1,3
          coords(k,15) = half*(cordel(k,6,ielem) + cordel(k,7,ielem))
        end do
        do k = 1,3
          coords(k,16) = half*(cordel(k,7,ielem) + cordel(k,8,ielem))
        end do
        do k = 1,3
          coords(k,17) = half*(cordel(k,8,ielem) + cordel(k,5,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,18) = one4th*(cordel(k,5,ielem) + cordel(k,6,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem))
        end do

        ! face 3 (bottom)
        ! 5++++14+++6
        ! +         +
        ! +         +
        ! 20   21   19
        ! +         +
        ! +         +
        ! 1++++9++++2

        ! edge midpoints
        do k = 1,3
          coords(k,19) = half*(cordel(k,2,ielem) + cordel(k,6,ielem))
        end do
        do k = 1,3
          coords(k,20) = half*(cordel(k,1,ielem) + cordel(k,5,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,21) = one4th*(cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,5,ielem) + cordel(k,6,ielem))
        end do

        ! face 4 (top)
        ! 8++++16+++7
        ! +         +
        ! +         +
        ! 22   24   23
        ! +         +
        ! +         +
        ! 4++++11+++3

        ! edge mipdoints
        do k = 1,3
          coords(k,22) = half*(cordel(k,4,ielem) + cordel(k,8,ielem))
        end do
        do k = 1,3
          coords(k,23) = half*(cordel(k,3,ielem) + cordel(k,7,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,24) = one4th*(cordel(k,3,ielem) + cordel(k,4,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem))
        end do

        ! face 5 (right)
        ! 3++++23+++7
        ! +         +
        ! +         +
        ! 10   25   15
        ! +         +
        ! +         +
        ! 2++++19+++6

        ! face midpoint
        do k = 1,3
          coords(k,25) = one4th*(cordel(k,2,ielem) + cordel(k,3,ielem)
     .      + cordel(k,6,ielem) + cordel(k,7,ielem))
        end do

        ! face 6 (left)
        ! 4++++22+++8
        ! +         +
        ! +         +
        ! 12   26   18
        ! +         +
        ! +         +
        ! 1++++20+++5

        ! face midpoint
        do k = 1,3
          coords(k,26) = one4th*(cordel(k,1,ielem) + cordel(k,4,ielem)
     .      + cordel(k,5,ielem) + cordel(k,8,ielem))
        end do

        ! volume centroid
        do k = 1,3
          coords(k,27) = one8th*(
     .        cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,3,ielem) + cordel(k,4,ielem)
     .      + cordel(k,5,ielem) + cordel(k,6,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem) )
        end do

        npf = 4
        ! loop over subcontrol surfaces
        do ics = 1, nscs
          ! loop over vertices of scs
          do inode = 1, npf
            itrianglenode = HexEdgeFacetTable(inode,ics)
            ! set scs coordinates using node table
            do k = 1,3
              scscoords(k,inode) = coords(k,itrianglenode)
            end do
          end do
          ! compute quad area vector using triangle decomposition
          call quadAreaByTriangleFacets( scscoords(1,1),
     .      area_vec(1,ielem,ics) )
        end do

      end do

      return
      end

      subroutine hex_scv_det( nelem, npe, nscv,
     .      cordel, vol, err, nerr )
c
c***********************************************************************
c***********************************************************************
c
c description:
c     This  routine returns the volume of each of the 8 subcontrol
c     volumes of the hex element for each element.
c
c formal parameters - input:
c     nelem         int   number of elements in the workset
c     npe           int   number of nodes per element
c     nscv          int   number of sub control volumes
c     cordel        real  element local coordinates
c
c formal parameters - output:
c     vol           real  volume of each sub control volume
c     err           real  positive volume check (0. = no error, 1. = error)
c     nerr          int   element number which fails positive volume check
c
c***********************************************************************
c
      implicit none
      integer nelem, npe, nscv, nerr
      double precision cordel, vol, err
      dimension cordel(3,npe,nelem),
     .          vol(nelem,nscv),
     .          err(nelem)
c
      double precision coords, scvcoords, voltmp
      ! coordinates of each vertex defining the subcontrol volumes
      dimension coords(3,27)
      ! coordinates of the subcontrol volume
      dimension scvcoords(3,8)

      integer ielem, j, k, icv, inode
      ! this table defines the vertices that compose each scv
      integer HexSubcontrolNodeTable(8,8)
      data HexSubcontrolNodeTable /
     .  0, 8, 12, 11, 19, 20, 26, 25,
     .  8, 1, 9, 12, 20, 18, 24, 26,
     .  12, 9, 2, 10, 26, 24, 22, 23,
     .  11, 12, 10, 3, 25, 26, 23, 21,
     .  19, 20, 26, 25, 4, 13, 17, 16,
     .  20, 18, 24, 26, 13, 5, 14, 17,
     .  26, 24, 22, 23, 17, 14, 6, 15,
     .  25, 26, 23, 21, 16, 17, 15, 7 /

      ! loop over elements
      do ielem = 1,nelem
        ! element vertices
        do j = 1,8
          do k = 1,3
            coords(k,j) = cordel(k,j,ielem)
          end do
        end do

        ! face 1 (front)
        ! 4++++11+++3
        ! +         +
        ! +         +
        ! 12   13   10
        ! +         +
        ! +         +
        ! 1++++9++++2

        ! edge midpoints
        do k = 1,3
          coords(k,9) = 0.5d0*(cordel(k,1,ielem) + cordel(k,2,ielem))
        end do
        do k = 1,3
          coords(k,10) = 0.5d0*(cordel(k,2,ielem) + cordel(k,3,ielem))
        end do
        do k = 1,3
          coords(k,11) = 0.5d0*(cordel(k,3,ielem) + cordel(k,4,ielem))
        end do
        do k = 1,3
          coords(k,12) = 0.5d0*(cordel(k,4,ielem) + cordel(k,1,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,13) = 0.25d0*(cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,3,ielem) + cordel(k,4,ielem))
        end do

        ! face 2 (back)
        ! 8++++16+++7
        ! +         +
        ! +         +
        ! 17   18   15
        ! +         +
        ! +         +
        ! 5++++14+++6

        ! edge midpoints
        do k = 1,3
          coords(k,14) = 0.5d0*(cordel(k,5,ielem) + cordel(k,6,ielem))
        end do
        do k = 1,3
          coords(k,15) = 0.5d0*(cordel(k,6,ielem) + cordel(k,7,ielem))
        end do
        do k = 1,3
          coords(k,16) = 0.5d0*(cordel(k,7,ielem) + cordel(k,8,ielem))
        end do
        do k = 1,3
          coords(k,17) = 0.5d0*(cordel(k,8,ielem) + cordel(k,5,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,18) = 0.25d0*(cordel(k,5,ielem) + cordel(k,6,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem))
        end do

        ! face 3 (bottom)
        ! 5++++14+++6
        ! +         +
        ! +         +
        ! 20   21   19
        ! +         +
        ! +         +
        ! 1++++9++++2

        ! edge midpoints
        do k = 1,3
          coords(k,19) = 0.5d0*(cordel(k,2,ielem) + cordel(k,6,ielem))
        end do
        do k = 1,3
          coords(k,20) = 0.5d0*(cordel(k,1,ielem) + cordel(k,5,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,21) = 0.25d0*(cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,5,ielem) + cordel(k,6,ielem))
        end do

        ! face 4 (top)
        ! 8++++16+++7
        ! +         +
        ! +         +
        ! 22   24   23
        ! +         +
        ! +         +
        ! 4++++11+++3

        ! edge mipdoints
        do k = 1,3
          coords(k,22) = 0.5d0*(cordel(k,4,ielem) + cordel(k,8,ielem))
        end do
        do k = 1,3
          coords(k,23) = 0.5d0*(cordel(k,3,ielem) + cordel(k,7,ielem))
        end do

        ! face midpoint
        do k = 1,3
          coords(k,24) = 0.25d0*(cordel(k,3,ielem) + cordel(k,4,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem))
        end do

        ! face 5 (right)
        ! 3++++23+++7
        ! +         +
        ! +         +
        ! 10   25   15
        ! +         +
        ! +         +
        ! 2++++19+++6

        ! face midpoint
        do k = 1,3
          coords(k,25) = 0.25d0*(cordel(k,2,ielem) + cordel(k,3,ielem)
     .      + cordel(k,6,ielem) + cordel(k,7,ielem))
        end do

        ! face 6 (left)
        ! 4++++22+++8
        ! +         +
        ! +         +
        ! 12   26   18
        ! +         +
        ! +         +
        ! 1++++20+++5

        ! face midpoint
        do k = 1,3
          coords(k,26) = 0.25d0*(cordel(k,1,ielem) + cordel(k,4,ielem)
     .      + cordel(k,5,ielem) + cordel(k,8,ielem))
        end do

        ! volume centroid
        do k = 1,3
          coords(k,27) = 0.125d0*(
     .        cordel(k,1,ielem) + cordel(k,2,ielem)
     .      + cordel(k,3,ielem) + cordel(k,4,ielem)
     .      + cordel(k,5,ielem) + cordel(k,6,ielem)
     .      + cordel(k,7,ielem) + cordel(k,8,ielem) )
        end do

        ! loop over subcontrol volumes
        do icv = 1, nscv
          ! loop over each vertex of scv
          do inode = 1,8
            ! set coordinates of scv vertex using the node table
            do k = 1,3
              scvcoords(k,inode) = coords(k,
     .          HexSubcontrolNodeTable(inode,icv)+1)
            end do
          end do
          ! compute the volume using an equivalent polyhedron
          call hexVolumeByTriangleFacets(scvcoords(1,1),voltmp)
          vol(ielem,icv) = voltmp
          ! check for negative volume
          if ( vol(ielem,icv) < 0.0d0) then
            err(nelem) = 1.0d0
            nerr = ielem
          end if
        end do

      end do

      return
      end

      subroutine quad3d_scs_det( nelem, cordel, area )
c
c***********************************************************************
c***********************************************************************
c
c description:
c     This  routine returns the area vectors for each of the 4
c     subcontrol surfaces of the quad.
c
c formal parameters - input:
c     nelem         int   number of elements in the workset
c     cordel        real  element local coordinates
c
c formal parameters - output:
c     area          real  area vector: A = Ax i + Ay j + Az k
c                                        = Area*(i + j + k)
c    
c***********************************************************************
c

      implicit none
      
c
      integer nelem
      integer k,i, idim
      double precision cordel, area
      double precision p, e, c, half, one4th
      double precision dx13, dx24, dy13, dy24, dz13, dz24
c
      dimension cordel(3,4,nelem)
      dimension area(3,nelem,4)
c
      dimension p(3,4), e(3,4), c(3)
c
c----------- first executable line of quad3d_scs_det ---------
c
      half = 0.5d0
      one4th = 1.d0/4.d0
c
c ... INITIALIZE FACE DEFINITION
c
      do 100 k=1,nelem
c
        do 20 i = 1, 4
          do 10 idim = 1, 3
            p(idim, i) = cordel(idim, i, k)
 10       continue
 20     continue
c
        do 30 idim = 1, 3
c
          e(idim, 1) = ( p(idim, 1) + p(idim, 2) ) * half
          e(idim, 2) = ( p(idim, 2) + p(idim, 3) ) * half
          e(idim, 3) = ( p(idim, 3) + p(idim, 4) ) * half
          e(idim, 4) = ( p(idim, 4) + p(idim, 1) ) * half
c
          c(idim) = ( p(idim, 1) + p(idim, 2) + 
     *                p(idim, 3) + p(idim, 4)   ) * one4th
c
 30     continue
c
c ... CALCULATE SUBCONTROL VOLUME FACE AREAS ...
c     ... subcontrol volume face 1
c
        dx13 = c(1)    - p(1, 1)
        dx24 = e(1, 4) - e(1, 1)
        dy13 = c(2)    - p(2, 1)
        dy24 = e(2, 4) - e(2, 1)
        dz13 = c(3)    - p(3, 1)
        dz24 = e(3, 4) - e(3, 1)
c
        area(1, k, 1) = half * ( dz24 * dy13 - dz13 * dy24 )
        area(2, k, 1) = half * ( dx24 * dz13 - dx13 * dz24 )
        area(3, k, 1) = half * ( dy24 * dx13 - dy13 * dx24 )
c
c     ... subcontrol volume face 2
c
        dx13 = c(1)    - p(1, 2)
        dx24 = e(1, 1) - e(1, 2)
        dy13 = c(2)    - p(2, 2)
        dy24 = e(2, 1) - e(2, 2)
        dz13 = c(3)    - p(3, 2)
        dz24 = e(3, 1) - e(3, 2)
c
        area(1, k, 2) = half * ( dz24 * dy13 - dz13 * dy24 )
        area(2, k, 2) = half * ( dx24 * dz13 - dx13 * dz24 )
        area(3, k, 2) = half * ( dy24 * dx13 - dy13 * dx24 )
c
c     ... subcontrol volume face 3
c
        dx13 = c(1)    - p(1, 3)
        dx24 = e(1, 2) - e(1, 3)
        dy13 = c(2)    - p(2, 3)
        dy24 = e(2, 2) - e(2, 3)
        dz13 = c(3)    - p(3, 3)
        dz24 = e(3, 2) - e(3, 3)
c
        area(1, k, 3) = half * ( dz24 * dy13 - dz13 * dy24 )
        area(2, k, 3) = half * ( dx24 * dz13 - dx13 * dz24 )
        area(3, k, 3) = half * ( dy24 * dx13 - dy13 * dx24 )
c
c     ... subcontrol volume face 4
c
        dx13 = c(1)    - p(1, 4)
        dx24 = e(1, 3) - e(1, 4)
        dy13 = c(2)    - p(2, 4)
        dy24 = e(2, 3) - e(2, 4)
        dz13 = c(3)    - p(3, 4)
        dz24 = e(3, 3) - e(3, 4)
c
        area(1, k, 4) = half * ( dz24 * dy13 - dz13 * dy24 )
        area(2, k, 4) = half * ( dx24 * dz13 - dx13 * dz24 )
        area(3, k, 4) = half * ( dy24 * dx13 - dy13 * dx24 )
c
  100 continue
c
c---------- last executable line of quad3d_scs_det ----------
c
      return
      end

      subroutine hex_shape_fcn( npts, par_coord, shape_fcn )
c
c***********************************************************************
c***********************************************************************
c
c formal parameters - input:
c     npts          int   number of points to evaluate (usually
c                         the number of Gauss Points)
c     par_coord     real  parametric coordinates of the points to be
c                         evaluated (typically, the gauss pts) -1:1
c
c formal parameters - output:
c     shape_fcn     real  shape functions evaluated at the evaluation
c                         points
c
c***********************************************************************
      implicit none

      integer npts

      double precision par_coord
      double precision shape_fcn

      dimension par_coord(3,npts)
      dimension shape_fcn(8,npts)

      integer j
      double precision s1, s2, s3, half, one4th, one8th
c
      half = 1.d0/2.d0
      one4th = 1.d0/4.d0
      one8th = 1.d0/8.d0
c
      do j = 1,npts
         s1 = par_coord(1,j)
         s2 = par_coord(2,j)
         s3 = par_coord(3,j)
c
         shape_fcn(1,j) = one8th + one4th*(-s1 - s2 - s3)
     *                  + half*( s2*s3 + s3*s1 + s1*s2 ) - s1*s2*s3
         shape_fcn(2,j) = one8th + one4th*( s1 - s2 - s3)
     *                  + half*( s2*s3 - s3*s1 - s1*s2 ) + s1*s2*s3
         shape_fcn(3,j) = one8th + one4th*( s1 + s2 - s3)
     *                  + half*(-s2*s3 - s3*s1 + s1*s2 ) - s1*s2*s3
         shape_fcn(4,j) = one8th + one4th*(-s1 + s2 - s3)
     *                  + half*(-s2*s3 + s3*s1 - s1*s2 ) + s1*s2*s3
         shape_fcn(5,j) = one8th + one4th*(-s1 - s2 + s3)
     *                  + half*(-s2*s3 - s3*s1 + s1*s2 ) + s1*s2*s3
         shape_fcn(6,j) = one8th + one4th*( s1 - s2 + s3)
     *                  + half*(-s2*s3 + s3*s1 - s1*s2 ) - s1*s2*s3
         shape_fcn(7,j) = one8th + one4th*( s1 + s2 + s3)
     *                  + half*( s2*s3 + s3*s1 + s1*s2 ) + s1*s2*s3
         shape_fcn(8,j) = one8th + one4th*(-s1 + s2 + s3)
     *                  + half*( s2*s3 - s3*s1 - s1*s2 ) - s1*s2*s3
c
      end do
c
      return
      end     

      subroutine hex_derivative( npts, par_coord, deriv )
c
c***********************************************************************
c***********************************************************************
c
c formal parameters - input:
c     npts          int   number of points to evaluate (usually
c                         the number of Gauss Points)
c     ncoord        int   number of parametric coordinates in element
c     par_coord     real  parametric coordinates of the points to be
c                         evaluated (typically, the gauss pts)
c
c formal parameters - output:
c     deriv         real  shape function derivatives evaluated at
c                         evaluation points.
c
c***********************************************************************
c
c
      implicit none

      integer npts
      double precision par_coord, deriv

      dimension par_coord(3,npts)
      dimension deriv(3,8,npts)
c
      integer j

      double precision s1, s2, s3

      double precision s1s2
      double precision s2s3
      double precision s1s3
c
      double precision half, one4th
      half = 1.d0/2.d0
      one4th = 1.d0/4.d0
c
      do j = 1,npts
         s1 = par_coord(1,j)
         s2 = par_coord(2,j)
         s3 = par_coord(3,j)
c
         s1s2 = s1*s2
         s2s3 = s2*s3
         s1s3 = s1*s3
c shape function derivative in the s1 direction -
         deriv(1,1,j) = half*( s3 + s2 ) - s2s3 - one4th
         deriv(1,2,j) = half*(-s3 - s2 ) + s2s3 + one4th
         deriv(1,3,j) = half*(-s3 + s2 ) - s2s3 + one4th
         deriv(1,4,j) = half*(+s3 - s2 ) + s2s3 - one4th
         deriv(1,5,j) = half*(-s3 + s2 ) + s2s3 - one4th
         deriv(1,6,j) = half*(+s3 - s2 ) - s2s3 + one4th
         deriv(1,7,j) = half*(+s3 + s2 ) + s2s3 + one4th
         deriv(1,8,j) = half*(-s3 - s2 ) - s2s3 - one4th
c
c shape function derivative in the s2 direction -
         deriv(2,1,j) = half*( s3 + s1 ) - s1s3 - one4th
         deriv(2,2,j) = half*( s3 - s1 ) + s1s3 - one4th
         deriv(2,3,j) = half*(-s3 + s1 ) - s1s3 + one4th
         deriv(2,4,j) = half*(-s3 - s1 ) + s1s3 + one4th
         deriv(2,5,j) = half*(-s3 + s1 ) + s1s3 - one4th
         deriv(2,6,j) = half*(-s3 - s1 ) - s1s3 - one4th
         deriv(2,7,j) = half*( s3 + s1 ) + s1s3 + one4th
         deriv(2,8,j) = half*( s3 - s1 ) - s1s3 + one4th
c
c shape function derivative in the s3 direction -
         deriv(3,1,j) = half*( s2 + s1 ) - s1s2 - one4th
         deriv(3,2,j) = half*( s2 - s1 ) + s1s2 - one4th
         deriv(3,3,j) = half*(-s2 - s1 ) - s1s2 - one4th
         deriv(3,4,j) = half*(-s2 + s1 ) + s1s2 - one4th
         deriv(3,5,j) = half*(-s2 - s1 ) + s1s2 + one4th
         deriv(3,6,j) = half*(-s2 + s1 ) - s1s2 + one4th
         deriv(3,7,j) = half*( s2 + s1 ) + s1s2 + one4th
         deriv(3,8,j) = half*( s2 - s1 ) - s1s2 + one4th
c
      end do
c
      return
      end

      subroutine hex_gradient_operator( nelem, npe, 
     *    nint, deriv, cordel, gradop, det_j, err, nerr )
c       
c***********************************************************************
c***********************************************************************
c
c description:
c     This  routine returns the gradient operator, determinate of 
c     the Jacobian, and error count for an element workset of 3D 
c     subcontrol surface elements The gradient operator and the 
c     determinate of the jacobians are computed at the center of
c     each control surface (the locations for the integration rule
c     are at the center of each control surface).
c
c formal parameters - input:
c     nelem         int   number of elements in the workset
c     npe           int   number of nodes per element
c     nint          int   number of sub control surfaces
c     deriv         real  shape function derivatives evaluated at the
c                         integration stations
c     cordel        real  element local coordinates
c
c formal parameters - output:
c     gradop        real  element gradient operator at each integration
c                         station
c     det_j         real  determinate of the jacobian at each integration
c                         station
c     err           real  positive volume check (0. = no error, 1. = error)
c     nerr          int   element number which fails positive volume check
c
c***********************************************************************
c
      implicit none
      
      integer nelem, npe, nint, nerr
      double precision deriv, cordel, gradop, det_j, err
c
      dimension deriv(3,npe,nint),
     *     cordel(3,npe,nelem),
     *     gradop(3,npe,nelem,nint),
     *     det_j(nelem,nint),
     *     err(nelem)
c
      integer ke, ki, kn
      double precision dx_ds1, dx_ds2, dx_ds3
      double precision dy_ds1, dy_ds2, dy_ds3
      double precision dz_ds1, dz_ds2, dz_ds3
      double precision ds1_dx, ds1_dy, ds1_dz
      double precision ds2_dx, ds2_dy, ds2_dz
      double precision ds3_dx, ds3_dy, ds3_dz

      double precision test, denom, sum, realmin
      realmin = 2.2250738585072014d-308
c
      do ke = 1,nelem
         err(ke) = 0.d0
      end do
c
      do ki = 1,nint
         do ke = 1,nelem
            dx_ds1 = 0.d0
            dx_ds2 = 0.d0
            dx_ds3 = 0.d0
            dy_ds1 = 0.d0
            dy_ds2 = 0.d0
            dy_ds3 = 0.d0
            dz_ds1 = 0.d0
            dz_ds2 = 0.d0
            dz_ds3 = 0.d0
c 
c calculate the jacobian at the integration station -
            do kn = 1,npe              
c
               dx_ds1 = dx_ds1+deriv(1,kn,ki)*cordel(1,kn,ke)
               dx_ds2 = dx_ds2+deriv(2,kn,ki)*cordel(1,kn,ke)
               dx_ds3 = dx_ds3+deriv(3,kn,ki)*cordel(1,kn,ke)
c                                                           
               dy_ds1 = dy_ds1+deriv(1,kn,ki)*cordel(2,kn,ke)
               dy_ds2 = dy_ds2+deriv(2,kn,ki)*cordel(2,kn,ke)
               dy_ds3 = dy_ds3+deriv(3,kn,ki)*cordel(2,kn,ke)
c                                             
               dz_ds1 = dz_ds1+deriv(1,kn,ki)*cordel(3,kn,ke)
               dz_ds2 = dz_ds2+deriv(2,kn,ki)*cordel(3,kn,ke)
               dz_ds3 = dz_ds3+deriv(3,kn,ki)*cordel(3,kn,ke)
c
            end do
c
c calculate the determinate of the jacobian at the integration station -
            det_j(ke,ki) = dx_ds1*( dy_ds2*dz_ds3 - dz_ds2*dy_ds3 )
     *                   + dy_ds1*( dz_ds2*dx_ds3 - dx_ds2*dz_ds3 )
     *                   + dz_ds1*( dx_ds2*dy_ds3 - dy_ds2*dx_ds3 )
c
c protect against a negative or small value for the determinate of the 
c jacobian. The value of real_min (set in precision.par) represents 
c the smallest Real value (based upon the precision set for this 
c compilation) which the machine can represent - 
            test = det_j(ke,ki)
            if( test .le. 1.d6*realmin ) then
               test = 1.d0
               err(ke) = 1.d0
            end if
            denom = 1.d0/test
c
c compute the gradient operators at the integration station -
c
            ds1_dx = denom*(dy_ds2*dz_ds3 - dz_ds2*dy_ds3)
            ds2_dx = denom*(dz_ds1*dy_ds3 - dy_ds1*dz_ds3)
            ds3_dx = denom*(dy_ds1*dz_ds2 - dz_ds1*dy_ds2)
c
            ds1_dy = denom*(dz_ds2*dx_ds3 - dx_ds2*dz_ds3)
            ds2_dy = denom*(dx_ds1*dz_ds3 - dz_ds1*dx_ds3)
            ds3_dy = denom*(dz_ds1*dx_ds2 - dx_ds1*dz_ds2)
c
            ds1_dz = denom*(dx_ds2*dy_ds3 - dy_ds2*dx_ds3)
            ds2_dz = denom*(dy_ds1*dx_ds3 - dx_ds1*dy_ds3)
            ds3_dz = denom*(dx_ds1*dy_ds2 - dy_ds1*dx_ds2)
c
            do kn = 1,npe 
c
              gradop(1,kn,ke,ki) = 
     *             deriv(1,kn,ki)*ds1_dx
     *           + deriv(2,kn,ki)*ds2_dx
     *           + deriv(3,kn,ki)*ds3_dx
c       
              gradop(2,kn,ke,ki) = 
     *             deriv(1,kn,ki)*ds1_dy
     *           + deriv(2,kn,ki)*ds2_dy
     *           + deriv(3,kn,ki)*ds3_dy
c       
              gradop(3,kn,ke,ki) = 
     *             deriv(1,kn,ki)*ds1_dz
     *           + deriv(2,kn,ki)*ds2_dz
     *           + deriv(3,kn,ki)*ds3_dz
c
            end do
         end do
      end do
c
c summarize volume error checks - 
      sum = 0.d0
      do ke = 1 , nelem
         sum = sum + err(ke)
      end do
      if( sum .ne. 0.d0 ) then
c flag error -
         do ke = 1 , nelem
            if ( err(ke) .ne. 0.d0 ) nerr = ke
         end do
      end if
c     
      return
      end

      subroutine quad3d_shape_fcn( npts, 
     *     par_coord, shape_fcn )
c       
c***********************************************************************
c***********************************************************************
c
c formal parameters - input:
c     npts          int   number of points to evaluate (usually 
c                         the number of Gauss Points)
c     par_coord     real  parametric coordinates of the points to be
c                         evaluated (typically, the gauss pts)
c
c formal parameters - output:
c     shape_fcn     real  shape functions evaluated at the evaluation
c                         points
c
c***********************************************************************
c
      implicit none
c
      integer npts
      double precision par_coord, shape_fcn

      dimension par_coord(2,npts)
      dimension shape_fcn(4,npts)
c
      integer j
      double precision s1, s2, one4th, half
c
      one4th = 1.d0/4.d0
      half = 1.d0/2.d0
c
      do j = 1,npts
c
         s1 = par_coord(1,j)
         s2 = par_coord(2,j)
c
         shape_fcn(1,j) = one4th + half*(-s1 - s2) + s1*s2
         shape_fcn(2,j) = one4th + half*( s1 - s2) - s1*s2
         shape_fcn(3,j) = one4th + half*( s1 + s2) + s1*s2
         shape_fcn(4,j) = one4th + half*(-s1 + s2) - s1*s2
c
      end do
c
      return
      end
